"""background-level inference for phase 1 experiment D.

WHY THIS MODULE EXISTS.  `run_experiment_d._pairSummary` reports cell-level
descriptive fractions -- `frac_worse_than_null` over 3450 double cells.  those
3450 cells are NOT 3450 independent observations.  they are 46 parameter draws,
each contributing 25 cells that share one parameter vector, one fold, one
baseline steady state and one pair of single-perturbation scores.  a cell-level
p-value would overstate the sample size by roughly 75x.  the independent unit is
the BACKGROUND, and every estimate here is formed at that level.

SIGN CONVENTIONS, STATED ONCE.  the three nulls do not live on the same scale
and "worse than null" therefore does not have the same sign for all three:

  additive        excess = burden_12 - (burden_1 + burden_2 - burden_0)
                  burden scale, worse = MORE burden = excess > 0
  multiplicative  excess = log(burden_12) - log(burden_1*burden_2/burden_0)
                  log-burden scale, worse = excess > 0
  bliss           excess = survival_12 - survival_1*survival_2/survival_0
                  survival scale, worse = LESS survival = excess < 0

so for bliss a negative median excess and a `frac_worse_than_null` above 0.5 are
the SAME statement, not a contradiction.  `NULLS` below carries the direction
explicitly so no downstream caller has to remember it.

TIES ARE NOT "WORSE".  the comparison is strict on both sides.  this matters for
bliss, where `survival = clip(1 - burden/H, 0, 1)` puts a hard floor at zero:
whenever either single perturbation is already lethal the bliss prediction is
exactly 0 and an equally lethal double gives excess exactly 0.  such a cell is
structurally incapable of showing bliss synergy.  `bliss_informative` is the
subset where the null is not floored, and it is reported alongside, never
instead of, the full set.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# design constants -- these are properties of configs/phase1/experiment_d.json
# and are asserted against the data rather than assumed
# --------------------------------------------------------------------------

#: burden multipliers; 1.0 is the unperturbed level
BURDEN_LEVELS = (1.0, 1.2, 1.4, 1.7, 2.0, 2.5)
#: capacity multipliers; 1.0 is the unperturbed level
CAPACITY_LEVELS = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)
#: full factorial cells per (background, pair)
CELLS_PER_PAIR = len(BURDEN_LEVELS) * len(CAPACITY_LEVELS)          # 36
#: cells with neither factor at 1.0 -- the genuinely double perturbations
DOUBLE_CELLS_PER_PAIR = (len(BURDEN_LEVELS) - 1) * (len(CAPACITY_LEVELS) - 1)  # 25
#: the three perturbation pairs, in the order run_experiment_d.PAIRS declares them
PAIRS = ("influx_x_total_capacity", "influx_x_chaperone_only",
         "nascent_x_total_capacity")

#: bootstrap seed.  fixed so every rerun of the closure reproduces the intervals.
BOOTSTRAP_SEED = 20260801
#: bootstrap replicates
N_BOOT = 10000


@dataclass(frozen=True)
class Null:
    name: str
    excess_col: str
    #: "greater" -> worse means excess > 0 ; "less" -> worse means excess < 0
    worse: str
    scale: str
    definition: str

    def worseMask(self, values: np.ndarray) -> np.ndarray:
        """strict comparison; ties are not counted as worse."""
        v = np.asarray(values, dtype=float)
        return v > 0.0 if self.worse == "greater" else v < 0.0

    @property
    def worseSign(self) -> float:
        return 1.0 if self.worse == "greater" else -1.0


NULLS: Tuple[Null, ...] = (
    Null("additive", "excess_additive", "greater", "burden",
         "burden_12 - (burden_1 + burden_2 - burden_0); additivity of the "
         "increment over baseline"),
    Null("multiplicative", "log_excess_multiplicative", "greater", "log burden",
         "log(burden_12) - log(burden_1*burden_2/burden_0); multiplicativity of "
         "the fold change over baseline"),
    Null("bliss", "excess_bliss", "less", "survival",
         "survival_12 - survival_1*survival_2/survival_0 with "
         "survival = clip(1 - burden/H, 0, 1); bliss independence"),
)

NULL_BY_NAME: Dict[str, Null] = {n.name: n for n in NULLS}


# --------------------------------------------------------------------------
# cell subsets
# --------------------------------------------------------------------------

def doubleCells(df: pd.DataFrame) -> pd.DataFrame:
    """the genuinely double perturbations, exactly as `_pairSummary` selects them."""
    return df[~df["single_perturbation"].astype(bool)]


def subsetMask(d: pd.DataFrame, subset: str) -> np.ndarray:
    """named, prespecified cell subsets.

    `all`                   every double cell (the primary analysis set)
    `both_singles_damaging` both single perturbations raise burden above baseline.
                            outside this set "supra-additive" does not mean
                            synergy: if one single is protective the
                            multiplicative prediction drops BELOW the additive
                            one and is easier to exceed for a purely arithmetic
                            reason.
    `bliss_informative`     both singles leave nonzero survival, so the bliss
                            prediction is not pinned to the zero floor.
    `uncensored`            the double burden did not hit the censor.  the censor
                            is reached ONLY by trajectories that escaped, whose
                            true burden is unbounded and is recorded as the
                            ceiling instead -- so excluding them removes the most
                            strongly supra-additive cells.  this subset
                            conditions on the DOUBLE outcome and is therefore
                            conservative by construction, not neutral.
    `clean`                 both of the above.
    `converged`             the double integration reached a steady state
                            (final relative rate < 1e-6) rather than stopping at
                            t_end still moving.

    `both_singles_damaging` and `bliss_informative` condition only on the SINGLE
    perturbation scores, never on the double outcome, so neither can select for
    or against an interaction.  they do select on model behaviour, so the
    backgrounds surviving them are not a random subsample and the surviving
    count is always reported.
    """
    if subset == "all":
        return np.ones(len(d), dtype=bool)
    if subset == "clean":
        return subsetMask(d, "both_singles_damaging") & subsetMask(d, "uncensored")
    if subset == "both_singles_damaging":
        return ((d["burden_1"].to_numpy() > d["burden_0"].to_numpy())
                & (d["burden_2"].to_numpy() > d["burden_0"].to_numpy()))
    if subset == "bliss_informative":
        return ((d["survival_1"].to_numpy() > 0.0)
                & (d["survival_2"].to_numpy() > 0.0))
    if subset == "uncensored":
        return ~d["censored_12"].to_numpy().astype(bool)
    if subset == "converged":
        return d["status_12"].to_numpy() == "converged"
    raise ValueError(f"unknown subset '{subset}'")


# --------------------------------------------------------------------------
# background-level summaries
# --------------------------------------------------------------------------

def perBackground(d: pd.DataFrame, null: Null) -> pd.DataFrame:
    """one row per background: its own fraction worse and its own median excess.

    a background with no cells in the subset is dropped, and the caller is told
    how many were dropped -- silently averaging over a shrinking denominator is
    exactly how a subset analysis stops being comparable to the primary one.
    """
    rows = []
    for bg, g in d.groupby("background", sort=True):
        v = g[null.excess_col].replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
        if not len(v):
            continue
        worse = null.worseMask(v)
        rows.append(dict(background=int(bg), n_cells=int(len(v)),
                         n_worse=int(worse.sum()),
                         frac_worse=float(worse.mean()),
                         median_excess=float(np.median(v)),
                         majority_worse=bool(worse.mean() > 0.5)))
    return pd.DataFrame(rows).set_index("background") if rows else pd.DataFrame(
        columns=["n_cells", "n_worse", "frac_worse", "median_excess",
                 "majority_worse"]).rename_axis("background")


def perBackgroundCollapse(d: pd.DataFrame) -> pd.DataFrame:
    """synthetic collapse per background: both singles viable, the double not."""
    rows = []
    for bg, g in d.groupby("background", sort=True):
        v = g["synthetic_collapse"].to_numpy().astype(bool)
        rows.append(dict(background=int(bg), n_cells=int(len(v)),
                         n_collapse=int(v.sum()),
                         frac_collapse=float(v.mean()),
                         any_collapse=bool(v.any())))
    return pd.DataFrame(rows).set_index("background") if rows else pd.DataFrame(
        columns=["n_cells", "n_collapse", "frac_collapse",
                 "any_collapse"]).rename_axis("background")


# --------------------------------------------------------------------------
# grouped (cluster) bootstrap
# --------------------------------------------------------------------------

def bootstrapIndices(n: int, n_boot: int = N_BOOT,
                     seed: int = BOOTSTRAP_SEED) -> np.ndarray:
    """the (n_boot, n) matrix of resampled group positions.

    shared by the reference and the fast path so that, for the same seed, the
    two take literally the same resamples and can be asserted equal rather than
    merely close.
    """
    return np.random.default_rng(seed).integers(0, int(n), size=(int(n_boot), int(n)))


def _percentileCI(draws: np.ndarray, alpha: float) -> Tuple[float, float]:
    lo, hi = np.nanpercentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi)


def scalarBootstrap(values: Sequence[float], stat: Callable[[np.ndarray], float],
                    n_boot: int = N_BOOT, seed: int = BOOTSTRAP_SEED,
                    alpha: float = 0.05) -> Dict[str, float]:
    """grouped bootstrap for a statistic that depends only on per-group scalars.

    every statistic used in this closure -- median of per-background medians,
    median of per-background fractions, fraction of backgrounds with a majority,
    fraction of backgrounds with any collapse -- reduces each background to one
    number first.  resampling that vector is identical to resampling the whole
    backgrounds, and `testFastPathMatchesReference` pins the identity.
    """
    v = np.asarray(values, dtype=float)
    n = len(v)
    if n == 0:
        return dict(point=float("nan"), lo=float("nan"), hi=float("nan"),
                    n_groups=0, n_boot=int(n_boot), seed=int(seed))
    idx = bootstrapIndices(n, n_boot, seed)
    draws = np.array([stat(v[row]) for row in idx], dtype=float)
    lo, hi = _percentileCI(draws, alpha)
    return dict(point=float(stat(v)), lo=lo, hi=hi, n_groups=int(n),
                n_boot=int(n_boot), seed=int(seed))


def pooledRateBootstrap(counts: Sequence[float], sizes: Sequence[float],
                        n_boot: int = N_BOOT, seed: int = BOOTSTRAP_SEED,
                        alpha: float = 0.05) -> Dict[str, float]:
    """cell-level rate sum(k)/sum(n), with backgrounds as the resampling unit.

    the point estimate is the ordinary cell-level fraction -- the same number
    `_pairSummary` reports.  only the INTERVAL differs, and that is the whole
    point: it is roughly sqrt(cells per background) wider than an interval that
    pretended the cells were independent.
    """
    k = np.asarray(counts, dtype=float)
    m = np.asarray(sizes, dtype=float)
    n = len(k)
    if n == 0 or m.sum() == 0:
        return dict(point=float("nan"), lo=float("nan"), hi=float("nan"),
                    n_groups=0, n_boot=int(n_boot), seed=int(seed))
    idx = bootstrapIndices(n, n_boot, seed)
    draws = k[idx].sum(axis=1) / m[idx].sum(axis=1)
    lo, hi = _percentileCI(draws, alpha)
    return dict(point=float(k.sum() / m.sum()), lo=lo, hi=hi, n_groups=int(n),
                n_boot=int(n_boot), seed=int(seed))


def clusterBootstrap(groups: Sequence[np.ndarray],
                     stat: Callable[[List[np.ndarray]], float],
                     n_boot: int = N_BOOT, seed: int = BOOTSTRAP_SEED,
                     alpha: float = 0.05) -> Dict[str, float]:
    """percentile CI resampling WHOLE BACKGROUNDS with replacement.

    `groups[i]` holds every cell value belonging to background i.  resampling at
    the group level is what keeps the within-background correlation in the
    interval; resampling cells would treat 25 correlated cells as 25 independent
    ones and produce an interval roughly sqrt(25) times too narrow.

    the point estimate is the statistic on the observed groups, not the
    bootstrap mean -- the bootstrap is used for the interval only.
    """
    groups = [np.asarray(g, dtype=float) for g in groups]
    n = len(groups)
    point = float(stat(groups)) if n else float("nan")
    if n == 0:
        return dict(point=float("nan"), lo=float("nan"), hi=float("nan"),
                    n_groups=0, n_boot=int(n_boot), seed=int(seed))
    idx = bootstrapIndices(n, n_boot, seed)
    draws = np.array([stat([groups[i] for i in row]) for row in idx], dtype=float)
    lo, hi = _percentileCI(draws, alpha)
    return dict(point=point, lo=lo, hi=hi, n_groups=int(n),
                n_boot=int(n_boot), seed=int(seed))


def statMedianOfGroupMedians(groups: List[np.ndarray]) -> float:
    return float(np.median([np.median(g) for g in groups]))


def statMedianOfGroupFractions(sign: float) -> Callable[[List[np.ndarray]], float]:
    """median across backgrounds of the within-background fraction worse."""
    def _f(groups: List[np.ndarray]) -> float:
        return float(np.median([float(((g > 0) if sign > 0 else (g < 0)).mean())
                                for g in groups]))
    return _f


def statFractionOfGroupsWithMajority(sign: float) -> Callable[[List[np.ndarray]], float]:
    """fraction of backgrounds whose own majority of cells is worse than null."""
    def _f(groups: List[np.ndarray]) -> float:
        return float(np.mean([float(((g > 0) if sign > 0 else (g < 0)).mean()) > 0.5
                              for g in groups]))
    return _f


def statPooledMean(groups: List[np.ndarray]) -> float:
    """cell-level rate, but with the interval formed by resampling backgrounds."""
    cat = np.concatenate(groups) if groups else np.array([])
    return float(cat.mean()) if len(cat) else float("nan")


def statFractionOfGroupsAny(groups: List[np.ndarray]) -> float:
    return float(np.mean([bool(g.any()) for g in groups]))


#: scalar-vector statistics, for use with `scalarBootstrap`
statMedian = lambda v: float(np.median(v))                          # noqa: E731
statMean = lambda v: float(np.mean(v))                              # noqa: E731
statFractionAboveHalf = lambda v: float(np.mean(v > 0.5))           # noqa: E731


def signTestBackgrounds(median_excess: Sequence[float], null: Null) -> Dict[str, float]:
    """exact binomial sign test over BACKGROUNDS -- the only legitimate p-value here.

    each background contributes one sign: the direction of its own median excess.
    ties (median exactly zero) are discarded, which is the conservative
    convention.  a cell-level p-value over 3450 cells would be wrong by roughly
    the square root of the cluster size and is never computed anywhere in this
    module.
    """
    from scipy.stats import binomtest
    v = np.asarray(median_excess, dtype=float)
    worse = int(null.worseMask(v).sum())
    better = int(null.worseMask(-v).sum())
    n = worse + better
    if n == 0:
        return dict(n_backgrounds=int(len(v)), n_worse=0, n_better=0, n_tied=int(len(v)),
                    p_two_sided=float("nan"))
    return dict(n_backgrounds=int(len(v)), n_worse=worse, n_better=better,
                n_tied=int(len(v) - n),
                p_two_sided=float(binomtest(worse, n, 0.5).pvalue))


def groupArrays(d: pd.DataFrame, column: str) -> List[np.ndarray]:
    """cell values grouped by background, in ascending background order."""
    return [g[column].to_numpy(dtype=float)
            for _, g in d.groupby("background", sort=True)]


# --------------------------------------------------------------------------
# sensitivity to backgrounds that produced no interaction rows
# --------------------------------------------------------------------------

def unresolvedBounds(k: int, n_usable: int, n_unresolved: int, n_unusable: int,
                     n_requested: int) -> Dict[str, float]:
    """transparent lower/upper bounds; nothing is imputed.

    three denominators, three different questions:

    `conditional`   k / n_usable.  the estimate for the population the experiment
                    can actually speak about: backgrounds with a viable baseline
                    below threshold.  this is the scientifically meaningful one.
    `usable_bounds` the two unresolved backgrounds might have been usable.  the
                    denominator is then n_usable + n_unresolved and the numerator
                    is k (both negative) or k + n_unresolved (both positive).  if
                    instead both were unusable the estimate stays at k/n_usable,
                    so the reported interval is the union.
    `requested_bounds`
                    k over all n_requested draws, i.e. counting the 12
                    structurally unusable backgrounds as backgrounds in which the
                    property was NOT demonstrated.  this is the most conservative
                    reading and it answers a different question -- what fraction
                    of an arbitrary draw from the sampled box shows the property.
    """
    if not 0 <= k <= n_usable:
        raise ValueError("k must lie in 0..n_usable")
    cond = k / n_usable if n_usable else float("nan")
    n_ext = n_usable + n_unresolved
    lo_u = min(k / n_ext, cond) if n_ext else float("nan")
    hi_u = max((k + n_unresolved) / n_ext, cond) if n_ext else float("nan")
    lo_r = k / n_requested if n_requested else float("nan")
    hi_r = (k + n_unresolved) / n_requested if n_requested else float("nan")
    return dict(k=int(k), n_usable=int(n_usable), n_unresolved=int(n_unresolved),
                n_unusable=int(n_unusable), n_requested=int(n_requested),
                conditional=float(cond),
                usable_lo=float(lo_u), usable_hi=float(hi_u),
                requested_lo=float(lo_r), requested_hi=float(hi_r))


# --------------------------------------------------------------------------
# the decision rule
# --------------------------------------------------------------------------

#: a claim of non-additivity/supra-additivity requires ALL THREE, in the same
#: direction.  any one of them alone is a descriptive statement, not a result.
DECISION_CRITERIA = (
    "background-level median of per-background median excess lies in the "
    "worse-than-null direction",
    "its grouped bootstrap 95% CI excludes zero",
    "the grouped bootstrap 95% CI for the fraction of backgrounds whose own "
    "majority of double cells is worse than null lies strictly above 0.5",
)


def supportDecision(median_ci: Dict[str, float], majority_ci: Dict[str, float],
                    null: Null) -> Dict[str, object]:
    """coherent-direction decision; returns the verdict and every failing leg."""
    sign = null.worseSign
    direction_ok = bool(sign * median_ci["point"] > 0.0)
    ci_excludes_zero = bool(median_ci["lo"] > 0.0 or median_ci["hi"] < 0.0)
    ci_direction_ok = bool(median_ci["lo"] > 0.0 if sign > 0 else median_ci["hi"] < 0.0)
    majority_ok = bool(majority_ci["lo"] > 0.5)
    legs = {
        "median_direction": direction_ok,
        "median_ci_excludes_zero": ci_excludes_zero and ci_direction_ok,
        "majority_ci_above_half": majority_ok,
    }
    return dict(supported=bool(all(legs.values())), legs=legs,
                failing=[k for k, v in legs.items() if not v])
