"""Transversality under load-dependent capacity, in closed form.

Task B5. When capacity is frozen, `F_j = (1, 0)` and the left null vector
normalised to `w_1 = 1` gives `w . F_j = 1`: transversality is automatic, which
is task B3. When capacity depends on the influx that argument fails, and the
earlier handoff said only "use the full `F_j`". It evaluates.

With `F_1 = j - R - G` and `F_2 = G`, both `R` and `G` now depending on `j`
through the capacity,

    F_j = ( 1 - R_j - G_j ,  G_j )

and with `w = (1, 1+lambda)` for `grad R = lambda grad G` at the fold,

    w . F_j = (1 - R_j - G_j) + (1 + lambda) G_j
            = 1 - R_j + lambda G_j.

So transversality is the scalar condition

    R_j - lambda G_j  !=  1.                                            (*)

This is exactly the statement that B3's automatic transversality is lost
precisely when capacity feeds back on the influx, and it is checkable: the two
derivatives are one-dimensional differences in `j` at a solved fold state, and
`lambda` comes from the gradient ratio.

TWO QUANTITIES, AND ONLY ONE IS A MARGIN. `|w . F_j|` in the normalisation
`w_1 = 1` proves NONVANISHING, but a left null vector rescales freely, so its
size means nothing on its own; it is reported as `coefficient`. The
scale-invariant quantity is the cosine `|w . F_j| / (||w|| ||F_j||)`, reported as
`margin`. They differ by more than presentation: the coefficient never drops
below 1.0 while the cosine reaches 0.0070.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402
import capacity_self_damage as SD  # noqa: E402

OUT = REPO_ROOT / "data" / "computed" / "transversality_selfdamage.json"
TSV = REPO_ROOT / "data" / "computed" / "transversality_selfdamage.tsv"


def transversalityAt(u: float, a: float, p: M.Params, sd: SD.SelfDamage,
                     h_rel: float = 1e-5) -> dict | None:
    """the two sides of (*) at one solved self-damaged fold state.

    `R_j` and `G_j` are partial derivatives in the INFLUX at fixed state, which
    is why they can be nonzero at all: the capacity factor carries `j`. In
    burden mode the factor does not contain `j`, so both should come out at the
    differencing floor and (*) should reduce to `1 != 0`. Reporting both modes
    is what makes that a check rather than an assumption.
    """
    try:
        Ru, Ra = SD._grad(SD.removalRsd, u, a, p, sd)
        Gu, Ga = SD._grad(SD.aggregateGsd, u, a, p, sd)
    except (M.ModelError, OverflowError):
        return None
    n_G2 = Gu * Gu + Ga * Ga
    if not np.isfinite(n_G2) or n_G2 <= 0.0:
        return None
    lam = float((Ru * Gu + Ra * Ga) / n_G2)

    h = h_rel * max(abs(p.j), 1e-8)
    pp, pm = p.with_(j=p.j + h), p.with_(j=max(p.j - h, 0.0))
    step = pp.j - pm.j
    if step <= 0.0:
        return None
    try:
        R_j = (SD.removalRsd(u, a, pp, sd) - SD.removalRsd(u, a, pm, sd)) / step
        G_j = (SD.aggregateGsd(u, a, pp, sd) - SD.aggregateGsd(u, a, pm, sd)) / step
    except (M.ModelError, OverflowError):
        return None

    w_dot_Fj = 1.0 - R_j + lam * G_j
    # the quantity above is taken in the normalisation w_1 = 1, which proves
    # NONVANISHING but is not a margin: a left null vector rescales freely. the
    # scale-invariant quantity is the cosine between w and F_j.
    n_w = float(np.hypot(1.0, 1.0 + lam))
    n_Fj = float(np.hypot(1.0 - R_j - G_j, G_j))
    return {
        "u": u, "a": a, "j": p.j, "eps": sd.eps, "mode": sd.mode,
        "lambda": lam, "R_j": float(R_j), "G_j": float(G_j),
        "w_dot_Fj": float(w_dot_Fj),
        "coefficient": float(abs(w_dot_Fj)),
        "margin": float(abs(w_dot_Fj) / max(n_w * n_Fj, 1e-300)),
        "norm_w": n_w, "norm_Fj": n_Fj,
        "grad_G": float(np.hypot(Gu, Ga)),
    }


def run(k: int = 20, seed: int = 11) -> dict:
    nets = SD._networks(k, seed)
    rows = []
    n_attempt = n_solved = 0
    for mode in (SD.INFLUX, SD.BURDEN):
        for eps in SD.EPS_LADDER:
            sd = SD.SelfDamage(eps=eps, mode=mode).validate()
            for item in nets:
                p0 = item[0] if isinstance(item, tuple) else item
                n_attempt += 1
                try:
                    out = SD.foldSolveSd(p0, sd)
                except (M.ModelError, OverflowError, ValueError,
                        np.linalg.LinAlgError):
                    out = None
                if out is None:
                    continue
                j_c, u, a = out
                n_solved += 1
                r = transversalityAt(u, a, p0.with_(j=j_c), sd)
                if r is not None:
                    rows.append(r)

    D = pd.DataFrame(rows)
    D.to_csv(TSV, sep="\t", index=False)

    out = {"n_solve_attempts": n_attempt, "n_solve_converged": n_solved,
           "n_evaluable": int(len(D))}
    for mode in (SD.INFLUX, SD.BURDEN):
        sub = D[D["mode"] == mode]
        if sub.empty:
            out[mode] = {"n": 0}
            continue
        out[mode] = {
            "n": int(len(sub)),
            "coefficient_min": float(sub["coefficient"].min()),
            "coefficient_max": float(sub["coefficient"].max()),
            "margin_min": float(sub["margin"].min()),
            "margin_median": float(sub["margin"].median()),
            "abs_R_j_max": float(sub["R_j"].abs().max()),
            "abs_G_j_max": float(sub["G_j"].abs().max()),
            "n_margin_below_1e-3": int((sub["margin"] < 1e-3).sum()),
            "per_eps_coefficient_min": {
                f"{e:g}": float(sub[sub["eps"] == e]["coefficient"].min())
                for e in sorted(sub["eps"].unique())
                if (sub["eps"] == e).any()},
            "per_eps_margin_min": {
                f"{e:g}": float(sub[sub["eps"] == e]["margin"].min())
                for e in sorted(sub["eps"].unique())
                if (sub["eps"] == e).any()},
        }
    if not D.empty:
        out["margin_min_overall"] = float(D["margin"].min())
    return out


def main() -> int:
    o = run()
    OUT.write_text(json.dumps(o, indent=2) + "\n")
    print(json.dumps(o, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
