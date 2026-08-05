"""Corollary 2's converse is stated with a condition that is false above n = 2.

Task B4. The manuscript says the converse of Corollary 2 needs "(G2)-(G4) plus
full rank". (G2) is `tr J != 0`, and in two dimensions that does two jobs at
once: with a zero eigenvalue present, the other eigenvalue IS `tr J`, so
`tr J != 0` makes the zero simple AND puts the rest of the spectrum off the
imaginary axis. Neither survives to `n >= 3`.

Writing the characteristic polynomial as

    p(lam) = lam^n + c_{n-1} lam^{n-1} + ... + c_1 lam + c_0,
    c_k = (-1)^{n-k} E_{n-k},   E_m = sum of the m x m principal minors,

a zero eigenvalue is `c_0 = 0`, and it is ALGEBRAICALLY SIMPLE exactly when
`c_1 != 0`. For n = 3 this is `c_1 = E_2`, the sum of the three 2x2 principal
minors, which has nothing to do with the trace: a matrix can have `tr J != 0` and
a double zero eigenvalue. Partial hyperbolicity is the separate requirement that
no other eigenvalue sit on the imaginary axis, which for n = 3 means the
remaining pair has nonzero real part -- a Hopf-type degeneracy is not excluded by
the trace either.

For n = 2, `c_1 = -tr J` and the remaining eigenvalue is `tr J`, so both reduce to
the published condition. That is why the error was invisible.

This script evaluates the corrected conditions at the fold states of the two
three-state systems Section 4.1 verifies -- the sigma-32 regulated chaperone pool
and the reactive/sequestered aggregate split -- and reports the margins.
"""

from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402
import regulation as REG  # noqa: E402
import sequestration as SEQ  # noqa: E402
import dilution as D  # noqa: E402

OUT = REPO_ROOT / "data" / "computed" / "corollary2_spectrum.json"


def spectrumConditions(J: np.ndarray) -> dict:
    """c_0, c_1, the spectrum, and the two conditions the converse needs."""
    n = J.shape[0]
    eig = np.linalg.eigvals(J)

    def elementary(m: int) -> float:
        """E_m = sum of the m x m principal minors."""
        return float(sum(np.linalg.det(J[np.ix_(idx, idx)])
                         for idx in combinations(range(n), m)))

    c0 = ((-1) ** n) * float(np.linalg.det(J))
    c1 = ((-1) ** (n - 1)) * elementary(n - 1)

    # the zero eigenvalue, and everything else
    order = np.argsort(np.abs(eig))
    eig = eig[order]
    rest = eig[1:]
    return {
        "n": n,
        "c0": c0,
        "c1": c1,
        "abs_c1": abs(c1),
        "tr_J": float(np.trace(J)),
        "det_J": float(np.linalg.det(J)),
        "lambda_zero_abs": float(abs(eig[0])),
        # partial hyperbolicity: how far the REST of the spectrum sits from the
        # imaginary axis. this is the quantity tr J does not control.
        "min_abs_re_rest": float(np.min(np.abs(rest.real))) if rest.size else np.nan,
        "max_abs_im_rest": float(np.max(np.abs(rest.imag))) if rest.size else np.nan,
        "eigs": [[float(z.real), float(z.imag)] for z in eig],
    }


# ---------------------------------------------------------------------------
# system 1: sigma-32 regulated chaperone pool
# ---------------------------------------------------------------------------


def regulatedFolds(k: int = 30, seed: int = 41, sigma0: float = 0.6) -> pd.DataFrame:
    """fold states of the sigma-32 system.

    THE `sigma0 = 0` ARM IS EXCLUDED, and the reason is not fastidiousness. With
    synthesis switched off, `dc/dt` is identically zero, so the third row of the
    3x3 Jacobian is identically zero, so `det J = 0` EVERYWHERE. `foldSolve3`'s
    third residual is then satisfied by construction and the solve converges to
    an arbitrary point of `{G = 0}` rather than to a fold. Its 24 "fold states"
    are not fold states, and a `c_1` margin computed over them measures nothing.

    The same structural zero makes the unregulated identity check of Section 4.1
    vacuous: `grad C = 0`, so `-det[grad R; grad G; grad C]` has a zero row and
    both sides of the identity are exactly zero. "median relative error
    0.000e+00" is `0 = 0`, not a verification. Compare the caution Section 5
    already states about the self-damage check.
    """
    run = FT.phase1RunDir()
    c_ = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c_ = c_[c_["C1_fold_exists"] == True]  # noqa: E712
    draws = c_.sample(n=k, random_state=seed)

    rows = []
    for idx, (_, r) in enumerate(draws.iterrows()):
        try:
            p = FT.paramsFromSampleRow(r)
        except (M.ModelError, ValueError, KeyError):
            continue
        for label, reg in (("regulated", REG.Regulator(sigma0=sigma0,
                                                       kappa_r=0.1, delta=0.5)),
                           ("regulated_strong", REG.Regulator(sigma0=1.5,
                                                              kappa_r=0.05,
                                                              delta=0.8))):
            try:
                out = REG.foldSolve3(p, reg, seed=(0.3, 0.1, p.c_tot))
            except (M.ModelError, OverflowError, np.linalg.LinAlgError, ValueError):
                continue
            if out is None:
                continue
            _, u, a, c = out
            try:
                J = REG.jacobian3(u, a, c, p, reg)
            except (M.ModelError, np.linalg.LinAlgError):
                continue
            if not np.all(np.isfinite(J)):
                continue
            row = spectrumConditions(np.asarray(J, float))
            row.update({"draw": idx, "regime": label})
            rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# system 2: reactive / sequestered aggregate split
# ---------------------------------------------------------------------------


def sequesteredFoldSolve(p: M.Params, seq: SEQ.Sequestration, g,
                         seed=(0.3, 0.1, 0.05)):
    """{da_r = 0, da_s = 0, det J = 0} for the reactive/sequestered system.

    the same structure as `regulation.foldSolve3`, one system over: three
    equations in three unknowns, none of which contains `j`.
    """
    from scipy.optimize import root

    def residual(z):
        u, ar, a_s = (float(np.exp(np.clip(v, -34, 6))) for v in z)
        try:
            _, dar, das = SEQ.rhs3s(u, ar, a_s, p, seq, g)
            d = float(np.linalg.det(SEQ.jacobian3s(u, ar, a_s, p, seq, g)))
        except (M.ModelError, OverflowError, np.linalg.LinAlgError, TypeError):
            return [1e6, 1e6, 1e6]
        if not all(np.isfinite(v) for v in (dar, das, d)):
            return [1e6, 1e6, 1e6]
        return [dar, das, d]

    s = root(residual, [np.log(v) for v in seed], method="hybr",
             options={"xtol": 1e-12})
    if not s.success:
        return None
    u, ar, a_s = (float(np.exp(np.clip(v, -34, 6))) for v in s.x)
    if not all(np.isfinite(v) and v > 0 for v in (u, ar, a_s)):
        return None
    try:
        _, dar, das = SEQ.rhs3s(u, ar, a_s, p, seq, g)
    except (M.ModelError, OverflowError, TypeError):
        return None
    if max(abs(dar), abs(das)) > 1e-7:
        return None
    return u, ar, a_s


def sequesteredFolds() -> pd.DataFrame:
    """fold states of the reactive/sequestered system over Section 4.1's grid.

    The `k_seq = 0` control is EXCLUDED: it is the two-state model with a third
    equation that is identically zero, so its Jacobian has a structural zero row
    and `c_1` vanishes for a reason that has nothing to do with the converse.
    Counting it would have manufactured the very failure this task looks for.
    """
    p = M.Params().validate()
    rows, n_attempt, n_solved = [], 0, 0
    for seq in SEQ.SEQ_GRID:
        if seq.k_seq == 0.0 and seq.k_rel == 0.0:
            continue
        for g in (None, D.Growth(0.05, 0.5), D.Growth(0.1, 0.5)):
            for seed in ((0.3, 0.1, 0.05), (0.8, 0.3, 0.15), (0.12, 0.03, 0.01)):
                n_attempt += 1
                out = sequesteredFoldSolve(p, seq, g, seed=seed)
                if out is None:
                    continue
                u, ar, a_s = out
                try:
                    J = np.asarray(SEQ.jacobian3s(u, ar, a_s, p, seq, g), float)
                except (M.ModelError, np.linalg.LinAlgError, TypeError):
                    continue
                if not np.all(np.isfinite(J)):
                    continue
                n_solved += 1
                row = spectrumConditions(J)
                row.update({"k_seq": seq.k_seq, "k_rel": seq.k_rel, "q": seq.q,
                            "u": u, "a_r": ar, "a_s": a_s})
                rows.append(row)
    D_ = pd.DataFrame(rows)
    D_.attrs["n_attempt"] = n_attempt
    D_.attrs["n_solved"] = n_solved
    return D_


def _margins(D_: pd.DataFrame, label: str) -> dict:
    if D_.empty:
        return {"n": 0}
    return {
        "n": int(len(D_)),
        "abs_c1_min": float(D_["abs_c1"].min()),
        "abs_c1_median": float(D_["abs_c1"].median()),
        "abs_tr_J_min": float(D_["tr_J"].abs().min()),
        "n_c1_below_1e-8": int((D_["abs_c1"] < 1e-8).sum()),
        "min_abs_re_rest": float(np.nanmin(D_["min_abs_re_rest"])),
        "n_rest_on_imaginary_axis": int(
            (D_["min_abs_re_rest"] < 1e-8).sum()),
        # the point of the task: does tr J track c_1 at all?
        "corr_abs_c1_abs_tr": float(
            np.corrcoef(D_["abs_c1"], D_["tr_J"].abs())[0, 1])
        if len(D_) > 2 else float("nan"),
    }


def run() -> dict:
    R = regulatedFolds()
    S = sequesteredFolds()
    out = {"regulated_chaperone": {
               lab: _margins(R[R["regime"] == lab], lab)
               for lab in sorted(R["regime"].unique())} if not R.empty else {},
           "sequestered_aggregate": _margins(S, "folds")}
    # report the counts nobody asserts against: how many solves were attempted
    # and how many converged. a margin over survivors is not a margin.
    out["sequestered_aggregate"]["n_solve_attempts"] = int(S.attrs.get("n_attempt", 0))
    out["sequestered_aggregate"]["n_solve_converged"] = int(S.attrs.get("n_solved", 0))
    if not S.empty:
        S.to_csv(REPO_ROOT / "data" / "computed" / "corollary2_sequestered.tsv",
                 sep="\t", index=False)
    if not R.empty:
        R.to_csv(REPO_ROOT / "data" / "computed" / "corollary2_regulated.tsv",
                 sep="\t", index=False)
    return out


def main() -> int:
    o = run()
    OUT.write_text(json.dumps(o, indent=2) + "\n")
    print(json.dumps(o, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
