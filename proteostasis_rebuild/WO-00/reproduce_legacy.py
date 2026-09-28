#!/usr/bin/env python3
"""
WO-00 gate G0.4: re-run the legacy model read-only and compare with its stored
outputs. this establishes WHAT the legacy computes. it says nothing about
whether the legacy is right -- that is WO-01..WO-09.

the vendored module is imported with bytecode writing disabled so that no
__pycache__ file is created or touched inside the legacy tree.
"""
import json
import sys

sys.dont_write_bytecode = True
from pathlib import Path  # noqa: E402

ENV = Path("/storage/kiran-stuff/proteostasis_law/envelope-paper")
sys.path.insert(0, str(ENV / "scripts" / "vendor"))
import two_pool_ode as m  # noqa: E402

TOL = 1e-6


def J_legacy(f, p):
    """the mapping used by legacy scripts 09/11/12 (includes (1-S))."""
    return f * p.N_prot * (1.0 - p.S_avg) * p.p_baseline / p.T_gen_s


def headroom_at(f):
    p = m.Params()
    P_dag, J_crit, mech, P_death = m.saddle_node_operational(m.J_curve_two, m.A_qs, p)
    P_star, A_star = m.steady_state(J_legacy(f, p), p)
    return {"P_dagger": P_dag, "J_crit": J_crit, "mechanism": mech,
            "P_death": P_death, "P_star": P_star, "A_star": A_star,
            "headroom_P": P_dag / P_star, "headroom_A": p.A_max / A_star,
            "f_codon_crit": m.f_codon_from_J(J_crit, p)}


def arithmetic(N, P, S, pm):
    """exact per-protein arithmetic threshold, written independently."""
    return (1.0 - P ** (1.0 / N)) / ((1.0 - S) * pm)


def run():
    comp = ENV / "data" / "computed"
    burden = json.loads((comp / "translation_burden.json").read_text())
    head = json.loads((comp / "headroom_sensitivity_summary.json").read_text())
    f_uw = burden["usage_weighted_mean_mu_per_codon"]

    r = headroom_at(f_uw)
    base = headroom_at(1e-4)
    arith = json.loads((ENV / "data" / "raw" / "arithmetic_results.json").read_text())
    tp = json.loads(Path("/storage/kiran-stuff/proteostasis-P1/two_pool_results.json")
                    .read_text())["B_compare"]["two_pool"]
    n300 = [x for x in arith["B_length_sweep"] if x["N"] == 300][0]
    checks = [
        ("headroom_P at usage-weighted mu", r["headroom_P"],
         head["internally_consistent_headroom_P"]),
        ("headroom_A at usage-weighted mu", r["headroom_A"],
         head["internally_consistent_headroom_A"]),
        ("headroom_P at window bottom 1e-4", base["headroom_P"],
         head["previously_reported_headroom_P"]),
        ("baseline two-pool f_codon_crit", r["f_codon_crit"], tp["f_codon_crit"]),
        ("baseline two-pool P_dagger", r["P_dagger"], tp["P_dagger"]),
        ("arith exact, stated params", arithmetic(300, 0.7, 0.3, 0.3),
         n300["f_codon_crit"]),
        ("arith with (1-S)p_m forced to 1", arithmetic(300, 0.7, 0.0, 1.0),
         arith["A_reproduction"]["f_codon_crit_exact"]),
    ]
    out = {"mechanism_at_baseline": r["mechanism"],
           "legacy_A_reproduction_factor": arith["A_reproduction"]["factor"],
           "rows": []}
    ok = True
    for name, got, ref in checks:
        rel = abs(got - ref) / abs(ref)
        passed = rel < TOL
        ok &= passed
        out["rows"].append({"check": name, "recomputed": got, "stored": ref,
                            "rel_err": rel, "tol": TOL, "pass": passed})
    out["all_pass"] = ok
    out["detail_at_usage_weighted_mu"] = r
    return out


if __name__ == "__main__":
    res = run()
    for row in res["rows"]:
        print(f"{'PASS' if row['pass'] else 'FAIL'}  {row['check']:<38} "
              f"recomputed={row['recomputed']:.6e} stored={row['stored']:.6e} "
              f"rel={row['rel_err']:.1e}")
    print("mechanism at baseline:", res["mechanism_at_baseline"])
    Path(__file__).with_name("legacy_reproduction.json").write_text(
        json.dumps(res, indent=2, default=float))
    sys.exit(0 if res["all_pass"] else 1)
