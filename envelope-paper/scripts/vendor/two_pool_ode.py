"""
two-pool ODE analysis for P1 (Proteostasis Law paper).

State:
  P(t) = fraction of proteome misfolded but still monomeric (0..1)
  A(t) = fraction of proteome in aggregated form (0..1)

Equations:
  M = P · Prot_tot                        # misfolded monomer concentration
  C_free(P) = C_tot / (1 + M/K_d)          # free chaperone (competition)
  v_fold(P) = k_obs_max · C_free/(C_free + K_d)
  v_agg(P)  = k_agg · M                    # second-order; units /s
  Phi(P)    = 1 + v_agg/(v_fold + k_deg)   # chaperone-titration amplification
  J_in(P)   = J_bare · Phi(P)
  R(P)      = (k_deg + v_fold) · P
  drain(P)  = k_nuc · P · M/Prot_tot       # fraction/s into A; k_nuc = k_agg
  A_sat(A)  = A / (A + A_half)

  dP/dt = J_in − R − drain · (1 − A_sat)
  dA/dt = drain · (1 − A_sat) − k_clear · A

Quasi-steady elimination of A (dA/dt = 0):
  k_clear · A = drain · A_half / (A + A_half)
  → A² + A·A_half − drain·A_half/k_clear = 0
  → A_qs(P) = 0.5·(−A_half + √(A_half² + 4·drain·A_half/k_clear))

Saddle-node:
  The *mathematical* saddle-node finds argmax of J_curve_two(P) over
  P ∈ (0, 1).  The *operational* saddle-node restricts the domain to
  P ≤ P_death, where P_death is the smallest P at which A_qs(P) ≥ A_max.

  A_max is the fraction of aggregated proteome at which the cell dies
  (inclusion-body toxicity threshold ≈ 0.25; see LITERATURE_ANCHORS.md).
  Two mechanisms can set J_crit:
    - monomer_runaway     : interior peak of J_curve_two
    - aggregation_death   : peak hits P_death before an interior peak exists

The script runs Parts A–F from the task spec and dumps JSON + summary md.
"""

from __future__ import annotations

import json
import math
import sys
import time
import warnings
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.optimize import brentq, fsolve, minimize_scalar

warnings.filterwarnings("ignore", category=RuntimeWarning)

HERE = Path(__file__).resolve().parent


# ──────────────────────────────────────────────────────────────
# empirical E. coli length distribution (for Part D bootstrap)
# ──────────────────────────────────────────────────────────────
# loaded lazily so the script still imports if the TSV is missing
# or pandas is not available for non-MC interactive use

_EMPIRICAL_LENGTHS = None


def _load_empirical_lengths():
    global _EMPIRICAL_LENGTHS
    if _EMPIRICAL_LENGTHS is None:
        import pandas as pd
        df = pd.read_csv(
            '/storage/kiran-stuff/proteostasis-P1/ecoli_proteome_with_genes.tsv',
            sep='\t',
        )
        # column name may be "Length" (capital L)
        if 'Length' in df.columns:
            length_col = 'Length'
        else:
            length_col = [c for c in df.columns if c.lower() == 'length'][0]
        _EMPIRICAL_LENGTHS = df[length_col].astype(int).values
    return _EMPIRICAL_LENGTHS


# ──────────────────────────────────────────────────────────────
# parameters
# ──────────────────────────────────────────────────────────────

# fraction of proteome aggregated at which cell viability collapses
# (inclusion-body toxicity literature; see LITERATURE_ANCHORS.md)
A_MAX_DEATH = 0.25


@dataclass
class Params:
    # rates in /s, concentrations in M (SI for the rate constants)
    k_deg: float = 3e-4
    k_obs_max: float = 1e-2
    C_tot_uM: float = 50.0          # µM
    K_d_uM: float = 1.0             # µM
    k_agg_M_s: float = 1e3          # M⁻¹ s⁻¹
    Prot_tot_uM: float = 300.0      # µM
    T_gen_s: float = 3600.0         # s
    N_prot: float = 300.0
    p_baseline: float = 0.3
    S_avg: float = 0.30
    k_clear: float = 4e-4
    A_half: float = 0.2
    # cell-death aggregation threshold — operational cap on the valid
    # domain of the monomeric saddle-node search
    A_max: float = A_MAX_DEATH

    # derived conveniences in µM for P-scaled concentrations
    def C_tot(self) -> float: return self.C_tot_uM * 1e-6
    def K_d(self) -> float:   return self.K_d_uM * 1e-6
    def Prot_tot(self) -> float: return self.Prot_tot_uM * 1e-6


BASELINE = Params()


# ──────────────────────────────────────────────────────────────
# core model pieces (vectorized over P)
# ──────────────────────────────────────────────────────────────

def c_free(P, p: Params):
    M = P * p.Prot_tot_uM                 # µM
    return p.C_tot_uM / (1.0 + M / p.K_d_uM)

def v_fold(P, p: Params):
    cf = c_free(P, p)
    return p.k_obs_max * cf / (cf + p.K_d_uM)

def v_agg(P, p: Params):
    # k_agg [M⁻¹ s⁻¹] · M [M] → /s
    M = P * p.Prot_tot_uM * 1e-6
    return p.k_agg_M_s * M

def phi(P, p: Params):
    return 1.0 + v_agg(P, p) / (v_fold(P, p) + p.k_deg)

def R_clearance(P, p: Params):
    return (p.k_deg + v_fold(P, p)) * P

def drain_rate(P, p: Params):
    # k_nuc = k_agg by assumption; drain = k_agg · P · M / Prot_tot
    # with M = P · Prot_tot [in M], so drain = k_agg · P² · Prot_tot [/s]
    Prot_tot_M = p.Prot_tot_uM * 1e-6
    return p.k_agg_M_s * P * P * Prot_tot_M   # fraction/s

def A_qs(P, p: Params):
    # raw (uncapped) quasi-steady A.  The cap A_max is used ONLY to
    # restrict the saddle-node domain and in reporting; the J-curve
    # functional form itself uses this uncapped A_qs.
    d = drain_rate(P, p)
    Ah = p.A_half
    disc = Ah * Ah + 4.0 * d * Ah / p.k_clear
    return 0.5 * (-Ah + np.sqrt(disc))


# ──────────────────────────────────────────────────────────────
# single-pool reduction (A ≡ 0, no drain term)
# ──────────────────────────────────────────────────────────────

def fixed_point_residual_single(P, p: Params, J_bare: float):
    return J_bare * phi(P, p) - R_clearance(P, p)

def fixed_point_residual_two(P, p: Params, J_bare: float):
    A = A_qs(P, p)
    # effective loss from P into A at steady state equals k_clear·A (by balance)
    return J_bare * phi(P, p) - R_clearance(P, p) - p.k_clear * A


# ──────────────────────────────────────────────────────────────
# saddle-node via inverse-curve maximum
# ──────────────────────────────────────────────────────────────

def J_curve_single(P, p: Params):
    """J_bare value making P a fixed point of the single-pool model."""
    return R_clearance(P, p) / phi(P, p)

def J_curve_two(P, p: Params):
    """J_bare value making P a fixed point of the two-pool (with A quasi-steady)."""
    return (R_clearance(P, p) + p.k_clear * A_qs(P, p)) / phi(P, p)


def saddle_node(J_curve: Callable, p: Params,
                P_lo: float = 1e-6, P_hi: float = 0.999):
    """Mathematical saddle-node: argmax and max of J_curve(P) on (P_lo, P_hi)."""
    P_grid = np.geomspace(P_lo, P_hi, 4000)
    vals = J_curve(P_grid, p)
    i = int(np.nanargmax(vals))
    lo = P_grid[max(i - 2, 0)]
    hi = P_grid[min(i + 2, len(P_grid) - 1)]
    try:
        res = minimize_scalar(lambda P: -J_curve(P, p),
                              bounds=(max(lo, 1e-12), min(hi, 0.9999)),
                              method="bounded",
                              options={"xatol": 1e-10})
        P_star = float(res.x)
        J_star = float(-res.fun)
    except Exception:
        P_star = float(P_grid[i])
        J_star = float(vals[i])
    return P_star, J_star


def saddle_node_operational(J_curve: Callable, A_qs_fn: Callable, p: Params,
                            P_lo: float = 1e-6, P_hi: float = 0.999):
    """
    Operational saddle-node = min(monomer_runaway, aggregation_death).

    Restricts the J_curve domain to P ≤ P_death, where P_death is the
    smallest P with A_qs(P) ≥ p.A_max.  If the J-curve peak on that
    restricted domain lies at the right edge (and A truly crosses A_max),
    the bound is set by aggregation death; otherwise by monomer runaway.
    """
    P_grid = np.geomspace(P_lo, P_hi, 4000)
    A_vals = np.array([A_qs_fn(Pi, p) for Pi in P_grid])

    idx_death = np.where(A_vals >= p.A_max)[0]
    P_death = float(P_grid[idx_death[0]]) if len(idx_death) else float(P_hi)

    P_feasible = P_grid[P_grid <= P_death]
    if len(P_feasible) < 2:
        # A_max exceeded before the first grid point; degenerate — return
        # the first feasible point as the operational bound
        P_feasible = P_grid[:2]
    J_vals = J_curve(P_feasible, p)
    i_peak = int(np.nanargmax(J_vals))
    P_peak = P_feasible[i_peak]
    J_peak = J_vals[i_peak]

    # mechanism: aggregation_death if peak is at the right (P_death) edge
    # AND A actually crossed A_max (otherwise P_death == P_hi and the peak
    # landing at the edge just means monomer runaway out beyond the grid).
    at_edge = (i_peak == len(P_feasible) - 1) and (P_death < P_hi)
    mechanism = "aggregation_death" if at_edge else "monomer_runaway"

    # refine with minimize_scalar around the peak
    lo_idx = max(i_peak - 2, 0)
    hi_idx = min(i_peak + 2, len(P_feasible) - 1)
    try:
        res = minimize_scalar(
            lambda P: -J_curve(P, p),
            bounds=(max(P_feasible[lo_idx], 1e-12),
                    min(P_feasible[hi_idx], 0.9999)),
            method="bounded",
            options={"xatol": 1e-10},
        )
        P_star = float(res.x)
        J_star = float(-res.fun)
    except Exception:
        P_star = float(P_peak)
        J_star = float(J_peak)

    return P_star, J_star, mechanism, P_death


def f_codon_from_J(J_bare: float, p: Params) -> float:
    return J_bare * p.T_gen_s / (p.N_prot * (1.0 - p.S_avg) * p.p_baseline)


# ──────────────────────────────────────────────────────────────
# steady-state finder for given J_bare (low-P stable branch)
# ──────────────────────────────────────────────────────────────

def steady_state(J_bare: float, p: Params, model: str = "two"):
    """Return (P*, A*) on the low-P stable branch, or (nan, nan) if past fold."""
    if model == "single":
        f = lambda P: fixed_point_residual_single(P, p, J_bare)
    else:
        f = lambda P: fixed_point_residual_two(P, p, J_bare)

    # bracket: low-P where f > 0 (at P → 0, f → J_bare > 0); increase until sign change
    P_lo, P_hi = 1e-10, 0.5
    try:
        f_lo = f(P_lo); f_hi = f(P_hi)
        if f_lo * f_hi > 0:
            # scan for sign change
            P_grid = np.geomspace(P_lo, 0.95, 500)
            fvals = np.array([f(Pi) for Pi in P_grid])
            sign_change = np.where(np.diff(np.sign(fvals)) != 0)[0]
            if len(sign_change) == 0:
                return float("nan"), float("nan")
            k = sign_change[0]           # first (low-P) root is the stable one
            P_lo, P_hi = P_grid[k], P_grid[k + 1]
        P_star = brentq(f, P_lo, P_hi, xtol=1e-12, maxiter=200)
    except Exception:
        return float("nan"), float("nan")

    A_star = float(A_qs(P_star, p)) if model == "two" else 0.0
    return float(P_star), A_star


# ──────────────────────────────────────────────────────────────
# Part A — steady-state curves
# ──────────────────────────────────────────────────────────────

def part_A(p: Params):
    J_grid = np.geomspace(1e-7, 1e-2, 80)
    P_two, A_two = [], []
    P_one = []
    for J in J_grid:
        Ps, As = steady_state(J, p, "two")
        P_two.append(Ps); A_two.append(As)
        Pss, _ = steady_state(J, p, "single")
        P_one.append(Pss)

    P_dag_two, J_crit_two, mech_two, P_death_two = saddle_node_operational(
        J_curve_two, A_qs, p)
    P_dag_one, J_crit_one = saddle_node(J_curve_single, p)

    return {
        "J_grid": J_grid.tolist(),
        "P_two": P_two,
        "A_two": A_two,
        "P_single": P_one,
        "P_dagger_two": P_dag_two,
        "J_bare_crit_two": J_crit_two,
        "A_dagger_two": float(A_qs(P_dag_two, p)),
        "f_codon_crit_two": f_codon_from_J(J_crit_two, p),
        "mechanism_two": mech_two,
        "P_death_two": P_death_two,
        "P_dagger_single": P_dag_one,
        "J_bare_crit_single": J_crit_one,
        "f_codon_crit_single": f_codon_from_J(J_crit_one, p),
        "mechanism_single": "monomer_runaway",
        "P_death_single": None,
    }


# ──────────────────────────────────────────────────────────────
# Part B — compare single-pool vs two-pool (baseline)
# ──────────────────────────────────────────────────────────────

def part_B(p: Params):
    P1, J1 = saddle_node(J_curve_single, p)
    P2, J2, mech2, Pd2 = saddle_node_operational(J_curve_two, A_qs, p)
    return {
        "single": {
            "P_dagger": P1,
            "J_bare_crit": J1,
            "f_codon_crit": f_codon_from_J(J1, p),
            "mechanism": "monomer_runaway",
            "P_death": None,
        },
        "two_pool": {
            "P_dagger": P2,
            "J_bare_crit": J2,
            "A_dagger": float(A_qs(P2, p)),
            "f_codon_crit": f_codon_from_J(J2, p),
            "mechanism": mech2,
            "P_death": Pd2,
        },
        "shift": {
            "P_dagger_ratio_two_over_single": P2 / P1 if P1 > 0 else None,
            "J_crit_ratio_two_over_single": J2 / J1 if J1 > 0 else None,
            "f_codon_ratio_two_over_single": (
                f_codon_from_J(J2, p) / f_codon_from_J(J1, p)
                if J1 > 0 else None),
        },
    }


# ──────────────────────────────────────────────────────────────
# Part C — k_clear sensitivity
# ──────────────────────────────────────────────────────────────

def part_C(p: Params):
    rows = []
    for kc in (1e-5, 3e-5, 1e-4, 3e-4, 1e-3):
        pc = Params(**{**asdict(p), "k_clear": kc})
        P, J, mech, Pd = saddle_node_operational(J_curve_two, A_qs, pc)
        rows.append({
            "k_clear": kc,
            "P_dagger": P,
            "A_dagger": float(A_qs(P, pc)),
            "J_bare_crit": J,
            "f_codon_crit": f_codon_from_J(J, pc),
            "mechanism": mech,
            "P_death": Pd,
        })
    return rows


# ──────────────────────────────────────────────────────────────
# Part D — Monte Carlo over parameter envelope
# ──────────────────────────────────────────────────────────────

def sample_params(rng: np.random.Generator) -> Params:
    def logu(lo, hi): return math.exp(rng.uniform(math.log(lo), math.log(hi)))
    def u(lo, hi):    return rng.uniform(lo, hi)
    # N_prot drawn as empirical bootstrap from UniProt E. coli length
    # distribution (matches paired_mc.py to keep the two MCs consistent)
    _lens = _load_empirical_lengths()
    N_prot = int(rng.choice(_lens))
    return Params(
        k_deg=logu(1e-4, 1e-3),
        k_obs_max=logu(3e-3, 8.4e-2),
        C_tot_uM=u(30.0, 80.0),
        K_d_uM=logu(0.06, 2.0),
        k_agg_M_s=logu(3e2, 3e3),
        Prot_tot_uM=u(250.0, 350.0),
        T_gen_s=u(30 * 60, 180 * 60),
        N_prot=N_prot,
        p_baseline=u(0.1, 0.5),
        S_avg=u(0.25, 0.35),
        k_clear=logu(1e-5, 1e-3),
        A_max=u(0.15, 0.35),
    )


def part_D(n: int = 5000, seed: int = 17):
    rng = np.random.default_rng(seed)
    P_dag, A_dag, f_crit, mech_list = [], [], [], []
    n_success = 0
    n_fail = 0
    t0 = time.time()
    for i in range(n):
        p = sample_params(rng)
        try:
            Pd, Jd, mech, P_death = saddle_node_operational(J_curve_two, A_qs, p)
            if not (math.isfinite(Pd) and math.isfinite(Jd)) or Pd <= 0 or Jd <= 0:
                n_fail += 1
                continue
            Ad = float(A_qs(Pd, p))
            fc = f_codon_from_J(Jd, p)
            if not (math.isfinite(fc) and fc > 0):
                n_fail += 1
                continue
            P_dag.append(Pd); A_dag.append(Ad); f_crit.append(fc)
            mech_list.append(mech)
            n_success += 1
        except Exception:
            n_fail += 1
    dt = time.time() - t0

    def pct(arr, qs):
        a = np.array(arr)
        if len(a) == 0:
            return [float("nan") for _ in qs]
        return [float(np.nanpercentile(a, q)) for q in qs]

    P_arr = np.array(P_dag)
    f_arr = np.array(f_crit)
    mech_arr = np.array(mech_list)
    frac_arith_bound = float(np.mean((f_arr >= 1e-4) & (f_arr <= 1e-3)))
    frac_order = float(np.mean((f_arr >= 1e-4) & (f_arr <= 1e-2)))
    frac_plausible_P = float(np.mean((P_arr >= 0.03) & (P_arr <= 0.30)))

    # mechanism breakdown
    n_agg = int(np.sum(mech_arr == "aggregation_death"))
    n_mon = int(np.sum(mech_arr == "monomer_runaway"))
    frac_aggregation_death = (n_agg / n_success) if n_success else 0.0
    frac_monomer_runaway = (n_mon / n_success) if n_success else 0.0

    f_agg = f_arr[mech_arr == "aggregation_death"]
    f_mon = f_arr[mech_arr == "monomer_runaway"]

    return {
        "n_requested": n,
        "n_success": n_success,
        "n_fail": n_fail,
        "seconds": dt,
        "P_dagger": {
            "median": float(np.nanmedian(P_dag)) if P_dag else float("nan"),
            "p16_p84": pct(P_dag, [16, 84]),
            "p2p5_p97p5": pct(P_dag, [2.5, 97.5]),
        },
        "A_dagger": {
            "median": float(np.nanmedian(A_dag)) if A_dag else float("nan"),
            "p16_p84": pct(A_dag, [16, 84]),
            "p2p5_p97p5": pct(A_dag, [2.5, 97.5]),
        },
        "f_codon_crit": {
            "median": float(np.nanmedian(f_crit)) if f_crit else float("nan"),
            "p16_p84": pct(f_crit, [16, 84]),
            "p2p5_p97p5": pct(f_crit, [2.5, 97.5]),
        },
        "frac_f_crit_1e-4_1e-3": frac_arith_bound,
        "frac_f_crit_1e-4_1e-2": frac_order,
        "frac_P_dagger_0.03_0.30": frac_plausible_P,
        "mechanism": {
            "frac_aggregation_death": frac_aggregation_death,
            "frac_monomer_runaway": frac_monomer_runaway,
            "n_aggregation_death": n_agg,
            "n_monomer_runaway": n_mon,
        },
        "f_codon_crit_by_mechanism": {
            "aggregation_death": {
                "n": int(len(f_agg)),
                "median": float(np.nanmedian(f_agg)) if len(f_agg) else float("nan"),
                "p16_p84": pct(f_agg.tolist(), [16, 84]) if len(f_agg) else [float("nan")] * 2,
                "p2p5_p97p5": pct(f_agg.tolist(), [2.5, 97.5]) if len(f_agg) else [float("nan")] * 2,
            },
            "monomer_runaway": {
                "n": int(len(f_mon)),
                "median": float(np.nanmedian(f_mon)) if len(f_mon) else float("nan"),
                "p16_p84": pct(f_mon.tolist(), [16, 84]) if len(f_mon) else [float("nan")] * 2,
                "p2p5_p97p5": pct(f_mon.tolist(), [2.5, 97.5]) if len(f_mon) else [float("nan")] * 2,
            },
        },
    }


# ──────────────────────────────────────────────────────────────
# Part E — cross-check at observed E. coli rate
# ──────────────────────────────────────────────────────────────

def part_E(p: Params, f_codon: float = 1e-4):
    J = f_codon * p.N_prot * (1.0 - p.S_avg) * p.p_baseline / p.T_gen_s
    Ps, As = steady_state(J, p, "two")
    Pd, Jd, mech, P_death = saddle_node_operational(J_curve_two, A_qs, p)

    def safe_ratio(num, den):
        try:
            if den and math.isfinite(num) and math.isfinite(den) and den > 0:
                return float(num / den)
        except Exception:
            pass
        return None

    P_below_dagger = (
        bool(Ps < Pd) if math.isfinite(Ps) and math.isfinite(Pd) else None
    )
    A_below_max = (
        bool(As < p.A_max) if math.isfinite(As) else None
    )

    return {
        "f_codon": f_codon,
        "J_bare": J,
        "P_star": Ps,
        "A_star": As,
        "P_dagger": Pd,
        "A_max": p.A_max,
        "P_death": P_death,
        "mechanism_at_bound": mech,
        "P_below_dagger": P_below_dagger,
        "A_below_max": A_below_max,
        "P_headroom_ratio": safe_ratio(Pd, Ps),   # P_dagger / P_star
        "A_headroom_ratio": safe_ratio(p.A_max, As),  # A_max / A_star
        "note": "low-P stable branch at observed E. coli mistranslation rate",
    }


# ──────────────────────────────────────────────────────────────
# Part F — A_max sensitivity sweep
# ──────────────────────────────────────────────────────────────

def part_F(p_base: Params):
    rows = []
    for a_max in (0.15, 0.20, 0.25, 0.30, 0.35):
        p = Params(**{**asdict(p_base), "A_max": a_max})
        P, J, mech, Pd = saddle_node_operational(J_curve_two, A_qs, p)
        rows.append({
            "A_max": a_max,
            "P_dagger": P,
            "A_dagger": float(A_qs(P, p)),
            "J_bare_crit": J,
            "f_codon_crit": f_codon_from_J(J, p),
            "mechanism": mech,
            "P_death": Pd,
        })
    return rows


# ──────────────────────────────────────────────────────────────
# optional plotting
# ──────────────────────────────────────────────────────────────

def make_plots(partA: dict, out_dir: Path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[plot] matplotlib unavailable: {e}")
        return

    J = np.array(partA["J_grid"])
    P2 = np.array(partA["P_two"], dtype=float)
    A2 = np.array(partA["A_two"], dtype=float)
    P1 = np.array(partA["P_single"], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    ax.loglog(J, P2, label="two-pool P*")
    ax.loglog(J, P1, "--", label="single-pool P*")
    ax.axvline(partA["J_bare_crit_two"], color="C0", ls=":",
               label=f"J_crit two = {partA['J_bare_crit_two']:.2e} ({partA['mechanism_two']})")
    ax.axvline(partA["J_bare_crit_single"], color="C1", ls=":",
               label=f"J_crit single = {partA['J_bare_crit_single']:.2e}")
    ax.set_xlabel("J_bare (/s)")
    ax.set_ylabel("P*")
    ax.set_title("steady-state misfolded monomer fraction")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    ax = axes[1]
    ax.loglog(J, A2, color="C2", label="two-pool A*")
    ax.axhline(partA["A_dagger_two"], color="C2", ls=":",
               label=f"A_dagger = {partA['A_dagger_two']:.2e}")
    ax.set_xlabel("J_bare (/s)")
    ax.set_ylabel("A*")
    ax.set_title("steady-state aggregated fraction")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    out_path = out_dir / "two_pool_partA.png"
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[plot] saved {out_path}")


# ──────────────────────────────────────────────────────────────
# summary markdown
# ──────────────────────────────────────────────────────────────

def write_summary(results: dict, out_path: Path):
    B = results["B_compare"]
    C = results["C_kclear_sensitivity"]
    D = results["D_monte_carlo"]
    E = results["E_ecoli_crosscheck"]
    F = results["F_amax_sensitivity"]

    arith_lo, arith_hi = 1e-4, 1e-3
    f_crit_two = B["two_pool"]["f_codon_crit"]
    f_crit_one = B["single"]["f_codon_crit"]
    gap_lo = arith_lo / f_crit_two if f_crit_two else float("inf")
    gap_hi = arith_hi / f_crit_two if f_crit_two else float("inf")

    lines = []
    lines.append("# two-pool ODE — summary (v2, mass-balance cap)\n")
    lines.append("Generated by `two_pool_ode.py` against the literature-anchored")
    lines.append("parameter envelope in `LITERATURE_ANCHORS.md`.")
    lines.append("")
    lines.append("Operational saddle-node = min(monomer_runaway, aggregation_death).")
    lines.append("A_max (cell-death aggregation fraction) restricts the valid domain")
    lines.append("of the J-curve to P ≤ P_death, where A_qs(P_death) = A_max.\n")

    lines.append("## Part B — single-pool vs two-pool (baseline)\n")
    lines.append(f"- single-pool: P_dagger = {B['single']['P_dagger']:.3e}, "
                 f"J_crit = {B['single']['J_bare_crit']:.3e} /s, "
                 f"f_codon_crit = {B['single']['f_codon_crit']:.3e}, "
                 f"mechanism = {B['single']['mechanism']}")
    lines.append(f"- two-pool:    P_dagger = {B['two_pool']['P_dagger']:.3e}, "
                 f"A_dagger = {B['two_pool']['A_dagger']:.3e}, "
                 f"J_crit = {B['two_pool']['J_bare_crit']:.3e} /s, "
                 f"f_codon_crit = {B['two_pool']['f_codon_crit']:.3e}, "
                 f"mechanism = {B['two_pool']['mechanism']}, "
                 f"P_death = {B['two_pool']['P_death']:.3e}")
    lines.append(f"- P_dagger shift (two/single): "
                 f"{B['shift']['P_dagger_ratio_two_over_single']:.3f}")
    lines.append(f"- f_codon_crit shift (two/single): "
                 f"{B['shift']['f_codon_ratio_two_over_single']:.3f}\n")

    lines.append("## Part C — k_clear sensitivity\n")
    lines.append("| k_clear (/s) | P_dagger | A_dagger | f_codon_crit | mechanism |")
    lines.append("|---:|---:|---:|---:|:---|")
    for row in C:
        lines.append(f"| {row['k_clear']:.0e} | {row['P_dagger']:.3e} | "
                     f"{row['A_dagger']:.3e} | {row['f_codon_crit']:.3e} | "
                     f"{row['mechanism']} |")
    lines.append("")

    lines.append("## Part D — Monte Carlo envelope (two-pool, with A_max sampled)\n")
    lines.append(f"- samples: {D['n_success']} succeeded / "
                 f"{D['n_requested']} requested "
                 f"({100.0*D['n_success']/D['n_requested']:.1f}%), "
                 f"wall time {D['seconds']:.1f} s")
    lines.append(f"- P_dagger median = {D['P_dagger']['median']:.3e}, "
                 f"[16, 84] = {D['P_dagger']['p16_p84']}, "
                 f"[2.5, 97.5] = {D['P_dagger']['p2p5_p97p5']}")
    lines.append(f"- A_dagger median = {D['A_dagger']['median']:.3e}, "
                 f"[16, 84] = {D['A_dagger']['p16_p84']}, "
                 f"[2.5, 97.5] = {D['A_dagger']['p2p5_p97p5']}")
    lines.append(f"- f_codon_crit median = {D['f_codon_crit']['median']:.3e}, "
                 f"[16, 84] = {D['f_codon_crit']['p16_p84']}, "
                 f"[2.5, 97.5] = {D['f_codon_crit']['p2p5_p97p5']}")
    lines.append(f"- P(f_codon_crit ∈ [1e-4, 1e-3]) = "
                 f"{D['frac_f_crit_1e-4_1e-3']:.3f}")
    lines.append(f"- P(f_codon_crit ∈ [1e-4, 1e-2]) = "
                 f"{D['frac_f_crit_1e-4_1e-2']:.3f}")
    lines.append(f"- P(P_dagger ∈ [0.03, 0.30]) = "
                 f"{D['frac_P_dagger_0.03_0.30']:.3f}\n")

    lines.append("### Mechanism breakdown (Part D)\n")
    lines.append(f"- fraction bound by aggregation_death: "
                 f"{D['mechanism']['frac_aggregation_death']:.3f} "
                 f"(n = {D['mechanism']['n_aggregation_death']})")
    lines.append(f"- fraction bound by monomer_runaway : "
                 f"{D['mechanism']['frac_monomer_runaway']:.3f} "
                 f"(n = {D['mechanism']['n_monomer_runaway']})")
    for key, label in (("aggregation_death", "aggregation_death"),
                       ("monomer_runaway", "monomer_runaway")):
        sub = D["f_codon_crit_by_mechanism"][key]
        lines.append(f"- f_codon_crit | {label} (n={sub['n']}): "
                     f"median = {sub['median']:.3e}, "
                     f"[16, 84] = {sub['p16_p84']}, "
                     f"[2.5, 97.5] = {sub['p2p5_p97p5']}")
    lines.append("")

    lines.append("## Part E — cross-check at observed E. coli rate\n")
    lines.append(f"- f_codon = {E['f_codon']:.1e} → J_bare = {E['J_bare']:.3e} /s")
    lines.append(f"- P* = {E['P_star']:.3e}, A* = {E['A_star']:.3e}")
    lines.append(f"- P_dagger = {E['P_dagger']:.3e} ({E['mechanism_at_bound']}); "
                 f"A_max = {E['A_max']:.3e}")
    lines.append(f"- P* < P_dagger: {E['P_below_dagger']}; "
                 f"A* < A_max: {E['A_below_max']}")
    if E["P_headroom_ratio"] is not None:
        lines.append(f"- headroom: P_dagger/P* = ×{E['P_headroom_ratio']:.2e}")
    if E["A_headroom_ratio"] is not None:
        lines.append(f"- headroom: A_max/A*   = ×{E['A_headroom_ratio']:.2e}")
    lines.append("")

    lines.append("## Part F — A_max sensitivity (baseline params)\n")
    lines.append("| A_max | P_dagger | A_dagger | f_codon_crit | mechanism |")
    lines.append("|---:|---:|---:|---:|:---|")
    for row in F:
        lines.append(f"| {row['A_max']:.2f} | {row['P_dagger']:.3e} | "
                     f"{row['A_dagger']:.3e} | {row['f_codon_crit']:.3e} | "
                     f"{row['mechanism']} |")
    lines.append("")

    lines.append("## ODE bound vs arithmetic bound (gap diagnosis)\n")
    lines.append(f"- arithmetic bound window: [{arith_lo:.0e}, {arith_hi:.0e}] /codon")
    lines.append(f"- two-pool ODE operational bound (baseline): {f_crit_two:.3e} /codon")
    lines.append(f"- dominant mechanism at baseline: **{B['two_pool']['mechanism']}**")
    lines.append(f"- gap (arithmetic lower / ODE): ×{gap_lo:.3f}")
    lines.append(f"- gap (arithmetic upper / ODE): ×{gap_hi:.3f}")
    lines.append(f"- inverse gap (ODE / arithmetic lower): ×{1.0/gap_lo:.2f}")
    lines.append(f"- inverse gap (ODE / arithmetic upper): ×{1.0/gap_hi:.2f}")
    lines.append(f"- single-pool baseline for comparison: {f_crit_one:.3e} /codon\n")

    out_path.write_text("\n".join(lines))
    print(f"[summary] wrote {out_path}")


# ──────────────────────────────────────────────────────────────
# main
# ──────────────────────────────────────────────────────────────

def main():
    p = BASELINE
    results = {}

    print("[A] steady-state sweep over J_bare")
    results["A_sweep"] = part_A(p)

    print("[B] single-pool vs two-pool (baseline)")
    results["B_compare"] = part_B(p)

    print("[C] k_clear sensitivity")
    results["C_kclear_sensitivity"] = part_C(p)

    print("[D] Monte Carlo envelope (5000 samples)")
    results["D_monte_carlo"] = part_D(n=5000, seed=17)

    print("[E] cross-check at f_codon = 1e-4")
    results["E_ecoli_crosscheck"] = part_E(p, f_codon=1e-4)

    print("[F] A_max sensitivity sweep")
    results["F_amax_sensitivity"] = part_F(p)

    results["baseline_params"] = asdict(p)

    out_json = HERE / "two_pool_results.json"
    out_json.write_text(json.dumps(results, indent=2, default=float))
    print(f"[json] wrote {out_json}")

    make_plots(results["A_sweep"], HERE)
    write_summary(results, HERE / "two_pool_summary.md")


if __name__ == "__main__":
    main()
