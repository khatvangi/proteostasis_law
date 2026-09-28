#!/usr/bin/env python3
"""
WO-01 gates G1.2 and G1.4: dimensional and bookkeeping audit of the legacy
two-pool model (proteostasis-P1/two_pool_ode.py lines 4-25, 127-159, 261-262).

the legacy expressions are transcribed here as sympy, one per legacy function,
with the unit of every numeric variable taken from its name/comment in the
source (e.g. Prot_tot_uM is uM, k_agg_M_s is 1/(M s)). the literal 1e-6 the
legacy uses to turn uM into M is written as the conversion symbol UM_TO_M so the
checker can see the scale change.

two separate questions are asked of every term:
  1. dimension: does it have the unit its equation requires?
  2. bookkeeping: is it a flux with a donor pool and a receiver pool, or a
     declared external source/sink with a named physical process?
a term can pass (1) and fail (2). that is the point of asking both.
"""
import json
import math
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import units as u  # noqa: E402

sym = sp.symbols(
    "P A J_bare k_deg k_obs_max C_tot_uM K_d_uM k_agg_M_s Prot_tot_uM "
    "T_gen_s N_prot p_baseline S_avg k_clear A_half UM_TO_M f_codon")
(P, A, J_bare, k_deg, k_obs_max, C_tot_uM, K_d_uM, k_agg_M_s, Prot_tot_uM,
 T_gen_s, N_prot, p_baseline, S_avg, k_clear, A_half, UM_TO_M, f_codon) = sym

# P and A are "fractions of proteome" in the legacy: dimensionless, and so is
# J_bare's numerator. P_T (the denominator) is never a state.
LEGACY_UNITS = {
    "P": u.ONE, "A": u.ONE, "J_bare": u.PER_S,
    "k_deg": u.PER_S, "k_obs_max": u.PER_S, "k_clear": u.PER_S,
    "C_tot_uM": u.UM, "K_d_uM": u.UM, "Prot_tot_uM": u.UM,
    "k_agg_M_s": u.PER_M_PER_S, "UM_TO_M": u.M_PER_UM,
    "T_gen_s": u.S, "N_prot": u.CODON_PER_PROTEIN,
    "p_baseline": u.ONE, "S_avg": u.ONE, "A_half": u.ONE,
    "f_codon": u.U(codon=-1, protein=1),   # errors (in proteins affected) per codon
}

# legacy pieces, line numbers from two_pool_ode.py
M_ = P * Prot_tot_uM                                      # :128
c_free = C_tot_uM / (1 + M_ / K_d_uM)                     # :129
v_fold = k_obs_max * c_free / (c_free + K_d_uM)           # :133
v_agg = k_agg_M_s * P * Prot_tot_uM * UM_TO_M             # :137-138
phi = 1 + v_agg / (v_fold + k_deg)                        # :141
R = (k_deg + v_fold) * P                                  # :144
drain = k_agg_M_s * P * P * Prot_tot_uM * UM_TO_M         # :150
A_sat = A / (A + A_half)                                  # :17
A_qs = sp.Rational(1, 2) * (-A_half + sp.sqrt(A_half**2 + 4 * drain * A_half / k_clear))  # :158-159
dPdt = J_bare * phi - R - drain * (1 - A_sat)             # :19
dAdt = drain * (1 - A_sat) - k_clear * A                  # :20
f_from_J = J_bare * T_gen_s / (N_prot * (1 - S_avg) * p_baseline)  # :262

PIECES = [
    ("M (misfolded conc.)", M_, u.UM, "two_pool_ode.py:128"),
    ("C_free", c_free, u.UM, ":129"),
    ("v_fold", v_fold, u.PER_S, ":133"),
    ("v_agg", v_agg, u.PER_S, ":137-138"),
    ("Phi", phi, u.ONE, ":141"),
    ("R (clearance flux)", R, u.PER_S, ":144"),
    ("drain", drain, u.PER_S, ":150"),
    ("A_qs", A_qs, u.ONE, ":158-159"),
    ("dP/dt", dPdt, u.PER_S, ":19"),
    ("dA/dt", dAdt, u.PER_S, ":20"),
    ("f_codon from J", f_from_J, u.U(codon=-1, protein=1), ":262"),
]

# bookkeeping. each flux: which pool loses it, which gains it. "external"
# means outside the modelled pools; an external source must name a process.
FLUXES = [
    ("J_bare", "external:synthesis of error-bearing chains", "P", True,
     "source is protein synthesis, but synthesis itself is not modelled; see ln2 finding"),
    ("J_bare*(Phi-1)", None, "P", False,
     "extra inflow proportional to v_agg; no pool is debited and no synthesis "
     "process is named. mass enters P from nowhere"),
    ("k_deg*P", "P", "external:degradation", True, "degradation sink"),
    ("v_fold*P", "P", "untracked native pool", True,
     "refolded protein leaves to a native pool that is not a state; harmless "
     "only because nothing depends on the native pool"),
    ("drain*(1-A_sat)", "P", "A", True, "paired transfer P -> A"),
    ("k_clear*A", "A", "external:clearance", True,
     "disaggregation (Mogk) would return A to P or N; modelled as total removal"),
    ("dilution mu*P, mu*A", None, None, False,
     "ABSENT. in growth at T_gen = 3600 s, mu = ln2/T_gen = 1.93e-4 /s, the same "
     "order as k_deg (3e-4) and k_clear (4e-4)"),
]


def run():
    rows, ok = [], True
    for name, expr, want, where in PIECES:
        try:
            got = u.unit_of(expr, LEGACY_UNITS)
            passed = got.same(want)
            rows.append({"piece": name, "where": where, "unit": str(got),
                         "expected": str(want), "dimension_ok": passed})
            ok &= passed
        except u.DimensionError as e:
            rows.append({"piece": name, "where": where, "unit": None,
                         "expected": str(want), "dimension_ok": False, "error": str(e)})
            ok = False

    # semantic check of the J mapping. f_codon*N*(1-S)*p is misfolded proteins
    # per synthesized protein. converting to a per-second fraction of the
    # proteome needs the synthesis rate per existing protein, which in balanced
    # exponential growth is mu = ln2/T_gen (no turnover), not 1/T_gen.
    synth_legacy = 1.0 / 3600.0
    synth_balanced = math.log(2) / 3600.0
    ln2 = {"legacy_synthesis_rate_per_protein": synth_legacy,
           "balanced_growth_rate_mu": synth_balanced,
           "legacy_over_balanced": synth_legacy / synth_balanced}

    # sensitivity of the dimension checker: a deliberately broken expression
    # (uM added to M) must be rejected, or the checker proves nothing.
    try:
        u.unit_of(Prot_tot_uM + Prot_tot_uM * UM_TO_M, LEGACY_UNITS)
        negative_control = False
    except u.DimensionError:
        negative_control = True

    return {"pieces": rows, "all_dimension_ok": ok,
            "fluxes": [dict(zip(("flux", "donor", "receiver", "conservative", "note"), f))
                       for f in FLUXES],
            "synthesis_rate_semantics": ln2,
            "negative_control_rejected": negative_control}


if __name__ == "__main__":
    res = run()
    for r in res["pieces"]:
        print(f"{'OK  ' if r['dimension_ok'] else 'FAIL'} {r['piece']:<22} {r['where']:<22} unit={r['unit']}")
    print("negative control (uM + M) rejected:", res["negative_control_rejected"])
    print("legacy 1/T_gen vs ln2/T_gen:", round(res["synthesis_rate_semantics"]["legacy_over_balanced"], 4))
    for f in res["fluxes"]:
        print(f"{'cons ' if f['conservative'] else 'NONC '} {f['flux']:<22} {f['note'][:70]}")
    Path(__file__).with_name("legacy_units_audit.json").write_text(json.dumps(res, indent=2))
    sys.exit(0 if res["all_dimension_ok"] and res["negative_control_rejected"] else 1)
