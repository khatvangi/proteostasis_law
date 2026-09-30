"""
writes ../EQUATION_TERM_LEDGER.tsv (the reader-test instrument, architecture
section 11) for every term introduced in stages 0-4.

the ledger is text written for a reader, so it is kept here as data rather than
generated from the flux table; check_conservation.c5_ledger then verifies that
every flux and dilution term has exactly one row and that each row's donor and
receiver are the pools the stoichiometry actually debits and credits.

parameter_status_WO07 annotates E. coli parameter status from
proteostasis_rebuild/WO-07/bundles.tsv (EXP = exponential bundle, STAT =
stationary bundle). it is an annotation only: no value enters the derivation.
"""
import csv
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "EQUATION_TERM_LEDGER.tsv"
COLS = ["term_id", "kind", "expression", "stage_introduced", "active_in", "process",
        "molecular_actor", "donor_pool", "receiver_pool", "rate_law_justification", "units",
        "evidence_class", "parameters", "parameter_status_WO07"]
DIL = "boundary: dilution (concentration, not molecules destroyed)"

# WO-07 parameter status strings (bundles.tsv, read 2026-09-30)
ST = {
    "s_P": "EXP MATCHED 0.478 uM/s (balanced-growth closure mu*P_T, neglects turnover); STAT UNMEASURED (mu*P_T fails at negative net growth)",
    "eps": "UNMEASURED in both bundles (synthesis_error_rate_per_codon and p_misfold UNMEASURED); the standing substitution frequency f is a different quantity (EXP: MISMATCHED PHASE_MIXED candidates only; STAT: MATCHED 1.82e-3/codon standing, not per synthesis)",
    "mu": "EXP MATCHED 1.611e-4 1/s (0.58 1/h); STAT UNMEASURED (candidate net growth -2.78e-6 1/s MISMATCHED; negative, so no dilution sink)",
    "k_mis": "UNMEASURED in both bundles (unfolding_k_mis)",
    "k_a": "UNMEASURED in both bundles (aggregation_k_a)",
    "k_dA": "UNMEASURED in both bundles (aggregate_degradation_k_dA)",
    "k_r": "no WO-07 row: stage-1 reference rate; its exact stage-2 limit is phi*k_on*C_T",
    "k_d": "UNMEASURED in both bundles (misfolded_degradation_k_d); exact stage-2 limit k_onD*D_T",
    "k_dis": "UNMEASURED in both bundles (ClpB_disaggregation_rate_k_dis); exact stage-2 limit k_onE*E_T",
    "s_C": "enters as C_T = s_C/mu. WO-07 forbids summed pools, so C is ONE machine: DnaK EXP MATCHED 10.88 uM monomer (GroEL14 EXP MATCHED 0.950 uM if C is read as GroEL); STAT UNMEASURED (candidates MISMATCHED: exponential normalization)",
    "k_on": "no WO-07 row for k_on; DnaK_affinity_K UNMEASURED in both bundles (in-vitro candidates 0.06-107 uM MISMATCHED IN_VITRO; an equilibrium K is not the driven-cycle K_M)",
    "k_off": "no WO-07 row for k_off (see k_on)",
    "k_cat": "DnaK_cycle_rate_k_cat UNMEASURED in both bundles (in-vitro 0.003-0.084 1/s MISMATCHED IN_VITRO); GroEL_cycle_rate UNMEASURED",
    "phi": "UNMEASURED in both bundles (cycle_partition_phi)",
    "prot": "no WO-07 row and no saved abundance (WO-06 schmidt_pools.json has no Lon, ClpXP or FtsH)",
    "s_E": "enters as E_T = s_E/mu. ClpB6 EXP MATCHED 0.01166 uM hexamer (full assembly assumed); STAT MISMATCHED 0.106 uM. An independent E pool is a counterfactual (ClpB needs DnaK; stage 5)",
    "k_catE": "UNMEASURED in both bundles (ClpB_disaggregation_rate_k_dis)",
    "k_onE": "no WO-07 row",
    "k_onA": "no WO-07 row; DnaK binding to aggregates in vivo not yet record-checked (WO-06 procedure)",
}
MA = "mass action, 1:1: rate proportional to the product of the two free concentrations"


def R(tid, kind, expr, st, act, proc, actor, don, rec, just, units, ev, pars, status):
    return dict(zip(COLS, [tid, kind, expr, str(st), act, proc, actor, don, rec, just, units, ev, pars, status]))


ALL = "S1 S2 S2C S4 S4C"
FIN = "S2 S2C S4 S4C"
ROWS = [
    R("syn_N", "flux", "(1-eps) s_P", 1, ALL, "synthesis of a chain that folds without help", "ribosome",
      "boundary: synthesis", "N", "constant flux: in balanced growth synthesis does not depend on the proteostasis state (declared; relaxed only by a growth-coupling stage)",
      "uM/s", "established in kind; constancy is a declared assumption", "s_P, eps", ST["s_P"] + " | eps: " + ST["eps"]),
    R("syn_U", "flux", "eps s_P", 1, ALL, "synthesis of a chain that is non-native at release (error-bearing or misfolded)", "ribosome (mistranslation; cotranslational misfolding)",
      "boundary: synthesis", "U", "fixed fraction eps of the same synthesis flux; no donor other than synthesis",
      "uM/s", "established in kind; eps not measured (P1-C25)", "s_P, eps", ST["eps"]),
    R("unfold", "flux", "k_mis N", 1, ALL, "spontaneous unfolding of native protein", "none (thermal)",
      "N", "U", "first-order: each native molecule unfolds independently at a constant hazard",
      "uM/s (k_mis 1/s)", "established in kind; value default 0 in WO-02", "k_mis", ST["k_mis"]),
    R("agg", "flux", "k_a U^2", 1, ALL, "aggregation of free non-native monomer", "none (encounter of two non-native monomers)",
      "U", "A", "second-order mass action: two free non-native monomers must meet; counted in monomer-equivalents; nucleation and templated growth not written (templated growth is stage 8)",
      "uM/s (k_a 1/(uM s))", "modelling choice (WO-02 self-review)", "k_a", ST["k_a"]),
    R("degA", "flux", "k_dA A", 1, ALL, "degradation of aggregated client", "unspecified (no named protease)",
      "A", "boundary: degradation", "first-order per aggregate monomer-equivalent; weak and flagged: no machine pool is written for it",
      "uM/s (k_dA 1/s)", "weak; flagged", "k_dA", ST["k_dA"]),
    R("resc", "flux", "k_r U", 1, "S1", "chaperone rescue with machinery in excess", "chaperone not limiting (DnaK/GroEL)",
      "U", "N", "first-order reference (unlimited capacity). exact limit of the stage-2 cycle when k_cat -> inf at fixed k_on, C_T: k_r = phi k_on C_T",
      "uM/s (k_r 1/s)", "reference, not a claim about cells", "k_r", ST["k_r"]),
    R("degU", "flux", "k_d U", 1, "S1 S2C S4C", "proteolysis of non-native monomer with protease in excess", "Lon, ClpXP (not limiting)",
      "U", "boundary: degradation", "first-order reference; exact limit of the stage-2 protease cycle when k_catD -> inf: k_d = k_onD D_T",
      "uM/s (k_d 1/s)", "reference; retained in S2C/S4C only to reproduce WO-02/WO-04", "k_d", ST["k_d"]),
    R("dis", "flux", "k_dis A", 1, "S1 S2C S4C", "disaggregation with disaggregase in excess", "ClpB (with DnaK in cells; not limiting)",
      "A", "U", "first-order reference; exact limit of the stage-2 disaggregase cycle when k_catE -> inf: k_dis = k_onE E_T",
      "uM/s (k_dis 1/s)", "reference", "k_dis", ST["k_dis"]),
    R("syn_C", "flux", "s_C", 2, FIN, "chaperone synthesis (constant: no induction before stage 7)", "ribosome (hsp genes)",
      "boundary: synthesis", "C", "constant flux; with dilution it fixes C_T = s_C/mu at steady state",
      "uM/s", "constancy is a declared assumption (stage 7 relaxes it)", "s_C", ST["s_C"]),
    R("bind", "flux", "k_on C U", 2, FIN, "chaperone captures a free non-native client", "DnaK/DnaJ; GroEL",
      "C + U", "B", MA, "uM/s (k_on 1/(uM s))", "kinetic form justified (P1-C05); values unverified", "k_on", ST["k_on"]),
    R("rel", "flux", "k_off B", 2, FIN, "unproductive dissociation of the complex", "DnaK/DnaJ; GroEL",
      "B", "C + U", "first-order dissociation of one complex", "uM/s (k_off 1/s)", "kinetic form standard for a 1:1 complex; value unverified", "k_off", ST["k_off"]),
    R("cycN", "flux", "phi k_cat B", 2, FIN, "completed cycle releases a native client", "DnaK/DnaJ/GrpE; GroEL/ES",
      "B", "C + N", "first-order cycle completion (ATP-driven) times the probability phi that the released client is native",
      "uM/s (k_cat 1/s)", "cycle form: WO-03 G3.2; phi unmeasured", "k_cat, phi", ST["k_cat"] + " | phi: " + ST["phi"]),
    R("cycU", "flux", "(1-phi) k_cat B", 2, FIN, "completed cycle releases a still non-native client", "DnaK/DnaJ/GrpE; GroEL/ES",
      "B", "C + U", "first-order cycle completion times 1 - phi", "uM/s", "as cycN", "k_cat, phi", ST["k_cat"] + " | phi: " + ST["phi"]),
    R("syn_D", "flux", "s_D", 2, "S2 S4", "protease synthesis", "ribosome (lon, clpXP genes)",
      "boundary: synthesis", "D", "constant flux; D_T = s_D/mu at steady state", "uM/s", "declared assumption", "s_D", ST["prot"]),
    R("bindD", "flux", "k_onD D U", 2, "S2 S4", "protease recognises a free non-native client", "Lon; ClpXP (ClpX recognition)",
      "D + U", "DU", MA, "uM/s (k_onD 1/(uM s))", "kinetic form by the same construction as the chaperone; values unmeasured", "k_onD", ST["prot"]),
    R("relD", "flux", "k_offD DU", 2, "S2 S4", "release before commitment to degradation", "Lon; ClpXP",
      "DU", "D + U", "first-order dissociation", "uM/s (k_offD 1/s)", "as bindD", "k_offD", ST["prot"]),
    R("catD", "flux", "k_catD DU", 2, "S2 S4", "processive degradation; protease freed", "Lon; ClpXP (ATP-driven)",
      "DU", "D + boundary: degradation", "first-order completion of one degradation cycle; the client monomer-equivalent leaves the client pool as peptides",
      "uM/s (k_catD 1/s)", "as bindD", "k_catD", ST["prot"]),
    R("syn_E", "flux", "s_E", 2, "S2 S4", "disaggregase synthesis", "ribosome (clpB gene)",
      "boundary: synthesis", "E", "constant flux; E_T = s_E/mu at steady state", "uM/s", "declared assumption", "s_E", ST["s_E"]),
    R("bindE", "flux", "k_onE E A", 2, "S2 S4", "disaggregase engages an aggregate site", "ClpB (counterfactual: DnaK-independent)",
      "A + E", "EA", MA + "; one site per aggregate monomer-equivalent (declared: exposed sites proportional to mass)",
      "uM/s (k_onE 1/(uM s))", "counterfactual pool (architecture Stage 2); site stoichiometry unmeasured", "k_onE", ST["k_onE"]),
    R("relE", "flux", "k_offE EA", 2, "S2 S4", "disaggregase dissociates without extraction", "ClpB",
      "EA", "A + E", "first-order dissociation", "uM/s (k_offE 1/s)", "as bindE", "k_offE", ST["k_onE"]),
    R("catE", "flux", "k_catE EA", 2, "S2 S4", "extraction of one monomer from the aggregate; disaggregase freed", "ClpB (ATP-driven threading)",
      "EA", "E + U", "first-order completion of one extraction cycle; the extracted monomer is non-native",
      "uM/s (k_catE 1/s)", "as bindE", "k_catE", ST["k_catE"]),
    R("bindA", "flux", "k_onA C A", 4, "S4 S4C", "chaperone binds an exposed site on an aggregate (feedback A)", "DnaK (with DnaJ)",
      "A + C", "CA", MA + "; one chaperone per aggregate monomer-equivalent site (declared: sites proportional to aggregate mass; for compact aggregates sites scale sub-linearly)",
      "uM/s (k_onA 1/(uM s))", "biology to be record-checked; not cited", "k_onA", ST["k_onA"]),
    R("relA", "flux", "k_offA CA", 4, "S4 S4C", "chaperone dissociates from the aggregate site", "DnaK",
      "CA", "A + C", "first-order dissociation; no catalytic exit at stage 4 (k_dcat = 0; that is stage 5)",
      "uM/s (k_offA 1/s)", "as bindA", "k_offA", ST["k_onA"]),
]
for pool, meaning, stg in [("N", "native client", 1), ("U", "free non-native client", 1), ("A", "aggregated client", 1),
                           ("B", "chaperone-client complex", 2), ("C", "free chaperone", 2),
                           ("DU", "protease-client complex", 2), ("D", "free protease", 2),
                           ("EA", "disaggregase-aggregate complex", 2), ("E", "free disaggregase", 2),
                           ("CA", "chaperone-aggregate complex", 4)]:
    ROWS.append(R("dil_" + pool, "dilution", f"mu {pool}", stg, "every stage containing " + pool,
                  f"dilution of {meaning} by growth", "cell growth and division", pool, DIL,
                  "exact under balanced exponential growth: every intracellular concentration carries -mu x; a complex dilutes both its client and its machine",
                  "uM/s (mu 1/s)", "exact kinematics; mu measured in EXP", "mu", ST["mu"]))
ND = "not a flux (derived quantity)"
DER = [
    ("C_T", "C + B + CA = s_C/mu", 2, "chaperone total, fixed by synthesis and dilution; at mu = 0 a measured constant", ST["s_C"]),
    ("D_T", "D + DU = s_D/mu", 2, "protease total", ST["prot"]),
    ("E_T", "E + EA = s_E/mu", 2, "disaggregase total", ST["s_E"]),
    ("K_M", "(k_off + k_cat + mu)/k_on", 2, "driven-cycle occupancy constant, exact at steady state (not K_d)", ST["k_on"]),
    ("K_MD", "(k_offD + k_catD + mu)/k_onD", 2, "protease occupancy constant", ST["prot"]),
    ("K_ME", "(k_offE + k_catE + mu)/k_onE", 2, "disaggregase occupancy constant", ST["k_onE"]),
    ("K_Ae", "(k_offA + mu)/k_onA", 4, "occupancy constant of chaperone on aggregate sites; dilution enters because the complex is diluted", ST["k_onA"]),
    ("S_0", "s_P (k_mis + eps mu)/(k_mis + mu)", 3, "net non-native inflow after eliminating N; equals s_P - mu N at the (1-eps) part", ST["s_P"]),
    ("kappa_B", "phi k_cat mu/(k_mis + mu) + mu", 3, "removal weight of chaperone-bound client: rescued protein counts as removed only if it is diluted before it unfolds again", ST["k_cat"]),
    ("G_eps", "s_P mu/(k_mis + mu)", 3, "dG/deps; constant in U", ST["s_P"]),
    ("G", "S_0 - kappa_B B - (k_catD+mu) DU - mu U - (k_dA+mu) A - mu EA - mu CA", 3, "organising function = client balance dP_T/dt on the manifold where every other pool is at steady state; CA term from stage 4", "derived"),
    ("F_A", "k_a U^2 - (k_dA+mu) A - (k_catE+mu) EA - mu CA", 3, "aggregate balance on the manifold; defines A(U)", "derived"),
    ("V_D", "k_catD D_T", 3, "protease capacity (maximum degradation flux)", ST["prot"]),
    ("V_E", "k_catE E_T", 3, "disaggregase capacity (maximum extraction flux)", ST["k_catE"]),
    ("s_crit", "V_D U_c/(K_MD + U_c), U_c = sqrt(V_E/k_a)", 3, "mu = 0, k_dA = 0 existence ceiling on s_P", "derived; no bundle value (s_P UNMEASURED in STAT)"),
]
for tid, expr, st, just, status in DER:
    ROWS.append(R(tid, "derived", expr, st, "S2 S4" if st >= 2 else ALL, just, "-", ND, ND, just,
                  "uM or uM/s or 1/s as written", "derived", "-", status))


def main():
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(ROWS)
    print(f"wrote {OUT} ({len(ROWS)} rows)")


if __name__ == "__main__":
    main()
