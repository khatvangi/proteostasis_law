"""WO-07 G7.3: the effect of each mismatch on model outputs, as a range.

each MISMATCHED row in bundles.tsv names effect ids; each id below varies ONE
mismatched quantity over a stated span, holds everything else at the bundle's
MATCHED value or at an explicitly ILLUSTRATIVE / declared reference value, and
reports the range of every output it touches. because the varied input is
mismatched by construction, every effect range is a SENSITIVITY (or
LEGACY_CONDITIONAL), never a prediction.

outputs
  O1 phi_sub      fraction of chains in the STANDING (copy-weighted) proteome
                  carrying >= 1 MS-detectable substitution:
                  sum_i w_i (1 - exp(-f L_i)) / sum_i w_i, w = copies. in balanced
                  growth (no turnover) this is also the fraction of new chains.
  O4 chains_per_X P_T / pool of ONE machine X.
  O2 sigma_DnaK   s_P phi_sub / (k_cat DnaK_T): substituted chains made per DnaK
                  cycle; the WO-02 saturation index with p_misfold = 1, one cycle
                  per client, standing frequency used as per-synthesis rate.
  O2s sP_break    k_cat DnaK_T / phi_sub: synthesis rate at which substituted
                  chains alone would equal DnaK in-vitro cycle capacity (used where
                  s_P is UNMEASURED).
  O3 headroom_P   legacy two-pool headroom P_dagger/P_star, WO-05 corrected
                  mapping. LEGACY_CONDITIONAL (phantom inflow, imposed A_max gate).
  O5 dilution     max(mu, 0) / (max(mu, 0) + k_d), k_d ILLUSTRATIVE. negative net
                  growth is not a dilution sink.
  O6 DnaK_bound   exact finite-pool bound fraction of DnaK (WO-03) at client
                  load M = eps P_T, eps ILLUSTRATIVE.

a REFERENCE output (not an effect) is DERIVED_OUTPUT only if every input it
depends on is MATCHED in that bundle and nothing illustrative is held.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for d in ("WO-03", "WO-05", "WO-06"):
    sys.path.insert(0, str(ROOT / d))
import derived as dv  # noqa: E402
import build_bundles as bb  # noqa: E402

LN2 = math.log(2.0)
ILLUSTRATIVE = {"k_d": 3e-4, "eps": 0.04, "K_ref_uM": 1.0, "norm_factor": [0.5, 2.0]}
SYN_PROXY = "standing substitution frequency used as per-synthesis error rate (UNMEASURED)"
HELD_ILL = {"O2_sigma_DnaK": ["p_misfold=1 (bound)", "one DnaK cycle per client", SYN_PROXY],
            "O2s_sP_breakeven": ["p_misfold=1 (bound)", "one DnaK cycle per client", SYN_PROXY],
            "O5_dilution_share": ["k_d=3e-4 /s (WO-02 SCENARIO, ILLUSTRATIVE)"],
            "O6_DnaK_bound_fraction": ["eps=0.04 (WO-02 SCENARIO, ILLUSTRATIVE)",
                                       "M = eps*P_T treated as free client"],
            "O3_legacy_headroom_P": ["legacy k_deg, k_agg, p_baseline, S, k_clear, A_half, "
                                     "A_max (all ASSUMED/UNVERIFIED, WO-06)", SYN_PROXY]}
# inputs each output depends on (context keys)
DEPS = {"O1_phi_sub": {"f", "Lw"}, "O2_sigma_DnaK": {"mu", "P_T", "f", "Lw", "k_cat", "DnaK"},
        "O2s_sP_breakeven": {"k_cat", "DnaK", "f", "Lw"}, "O3_legacy_headroom_P": {"*"},
        "O5_dilution_share": {"mu"}, "O6_DnaK_bound_fraction": {"DnaK", "P_T", "K"}}
# context keys that hold a MATCHED bundle value (everything else is a reference point)
MATCHED_KEYS = {"EXPONENTIAL": {"mu", "P_T", "pools", "DnaK", "Lw", "N_mean"},
                "STATIONARY": {"f"}}


# ---------------------------------------------------------------- outputs
def phi_sub(f, Lw):
    L_, w = Lw
    return float(np.sum(w * -np.expm1(-f * L_)) / np.sum(w))


def sigma(c):
    return c["mu"] * c["P_T"] * phi_sub(c["f"], c["Lw"]) / (c["k_cat"] * c["DnaK"])


def sP_break(c):
    return c["k_cat"] * c["DnaK"] / phi_sub(c["f"], c["Lw"])


def dnak_bound(c):
    import chaperone as ch
    return float(ch.bound_exact(c["DnaK"], ILLUSTRATIVE["eps"] * c["P_T"], c["K"]) / c["DnaK"])


def legacy_headroom(c):
    import error_semantics as es
    import run_wo05 as r5
    p = r5.m.Params(Prot_tot_uM=c["P_T"], C_tot_uM=c["DnaK"], T_gen_s=LN2 / c["mu"],
                    N_prot=c["N_mean"], K_d_uM=c["K"], k_obs_max=c["k_cat"])
    J = es.flux(es.ErrorRate(c["f"], es.MS), p.N_prot, p.p_baseline, p.T_gen_s, p.S_avg,
                balanced_growth=True)
    return float(r5.headroom(J, p)["headroom_P"])


def dilution(c):
    d = max(c["mu"], 0.0)
    return d / (d + ILLUSTRATIVE["k_d"])


OUT = {"O1_phi_sub": lambda c: phi_sub(c["f"], c["Lw"]),
       "O2_sigma_DnaK": sigma, "O2s_sP_breakeven": sP_break,
       "O3_legacy_headroom_P": legacy_headroom, "O5_dilution_share": dilution,
       "O6_DnaK_bound_fraction": dnak_bound}


def chains_per(machine):
    return lambda c: c["P_T"] / c["pools"][machine]


def effect_kind(name):
    return "LEGACY_CONDITIONAL" if name.startswith("O3") else "SENSITIVITY"


def reference_kind(name, bundle):
    if name.startswith("O3"):
        return "LEGACY_CONDITIONAL"
    if name.startswith("O4"):
        deps = {"P_T", "pools"}
    else:
        deps = DEPS[name]
    ok = deps <= MATCHED_KEYS[bundle] and not HELD_ILL.get(name)
    return "DERIVED_OUTPUT" if ok else "SENSITIVITY"


# ---------------------------------------------------------------- contexts
def contexts():
    per, der, _ = bb.load_sources()
    gl, st = per[bb.EXP_COND], per[bb.STAT_COND]
    pools = lambda cond: {m: per[cond][f"{g}_uM_oligomer"] for _, m, g, _ in bb.MACHINES}
    Le = dv.weighted_lengths(bb.EXP_COND)[:2]
    Ls = dv.weighted_lengths(bb.STAT_COND)[:2]
    Lg = dv.genome_lengths()
    et = der["etel_translated_weights_exp_glucose"]
    exp = {"mu": gl["growth_per_h"] / 3600, "P_T": gl["total_protein_mM_wholecell"] * 1e3,
           "pools": pools(bb.EXP_COND), "Lw": Le,
           "N_mean": der["N_synthesis_weighted_exp_glucose"]["weighted_mean_codons"],
           # reference points INSIDE mismatch spans (no matched value exists):
           "f": et["zero_inclusive_mean"], "k_cat": 0.04, "K": ILLUSTRATIVE["K_ref_uM"]}
    stat = {"mu": st["growth_per_h"] / 3600, "P_T": st["total_protein_mM_wholecell"] * 1e3,
            "pools": pools(bb.STAT_COND), "Lw": Ls, "f": 1.82e-3, "k_cat": 0.04,
            "K": ILLUSTRATIVE["K_ref_uM"]}
    for c in (exp, stat):
        c["DnaK"] = c["pools"]["DnaK"]
    refs = {"EXPONENTIAL": {"f": "eTEL zero-inclusive, translated weights (MIXED phase)",
                            "k_cat": "0.04 /s in vitro T->R", "K": "1 uM ILLUSTRATIVE"},
            "STATIONARY": {"Lw": "Schmidt copy-weighted 1 d (MISMATCHED)",
                           "mu": "Schmidt 1 d, -0.01 /h (MISMATCHED)",
                           "k_cat": "0.04 /s in vitro", "K": "1 uM ILLUSTRATIVE",
                           "P_T/pools": "Schmidt 1 d (MISMATCHED)"}}
    return exp, stat, {"Le": Le, "Ls": Ls, "Lg": Lg}, per, der, refs


def with_(c, **kw):
    d = dict(c)
    d.update(kw)
    if "pools" in kw:
        d["DnaK"] = kw["pools"]["DnaK"]
    return d


def evaluate(ctx, key, values, outputs):
    """outputs over the span; values are the context entries to substitute
    (key '_ctx' means each value is a whole context)."""
    res = {}
    for name, fn in outputs.items():
        ys = np.array([fn(v) if key == "_ctx" else fn(with_(ctx, **{key: v})) for v in values],
                      float)
        lo, hi = float(np.nanmin(ys)), float(np.nanmax(ys))
        res[name] = {"kind": effect_kind(name), "at_reference": float(fn(ctx)), "lo": lo,
                     "hi": hi, "fold_span": hi / lo if lo > 0 else None,
                     "illustrative_held": HELD_ILL.get(name, [])}
    return res


def grid(lo, hi, ref=None, n=9):
    v = list(np.geomspace(lo, hi, n)) if lo > 0 else list(np.linspace(lo, hi, n))
    return sorted(set(v + ([ref] if ref is not None else [])))


# ---------------------------------------------------------------- effects
def run():
    exp, stat, Ls_, per, der, refs = contexts()
    et, eg = der["etel_translated_weights_exp_glucose"], der["etel_genomic_weights"]
    ft = sorted(et.values())
    E = {}

    def add(eid, bundle, params, span, basis, ctx, key, values, outputs, span_kind="DATA",
            **extra):
        E[eid] = {"bundle": bundle, "params": params, "span": span, "span_basis": basis,
                  "span_kind": span_kind, "outputs": evaluate(ctx, key, values, outputs),
                  **extra}

    pick = lambda *ks: {k: OUT[k] for k in ks}
    exp_core = pick("O2_sigma_DnaK", "O3_legacy_headroom_P")
    o4 = {f"O4_chains_per_{m}": chains_per(m) for m in ("DnaK", "GroEL", "ClpB")}

    # ---- EXPONENTIAL
    bnid = 4000.0
    add("X_PT_BNID", "EXPONENTIAL", ["total_protein_P_T"], [exp["P_T"], bnid],
        "matched Schmidt glucose vs BNID 104726 (B/r, 40 min doubling, assumed volume)",
        exp, "P_T", [exp["P_T"], bnid], {**o4, **exp_core})
    milo = [x * 1e3 for x in json.loads((ROOT / "WO-06/audit_derived.json").read_text())
            ["prot_tot"]["milo2013_mM_from_2to4e6_per_um3"]]
    add("X_PT_MILO", "EXPONENTIAL", ["total_protein_P_T"], [exp["P_T"], milo[1]],
        "matched value through the Milo 2013 generic range", exp, "P_T",
        grid(milo[0], milo[1], exp["P_T"]), {**o4, **exp_core})
    add("X_N_GENOME", "EXPONENTIAL", ["copy_weighted_length_N"],
        ["copy-weighted (glucose)", "genome-unweighted"],
        "two weightings of the length distribution", exp, "_ctx",
        [exp, with_(exp, Lw=Ls_["Lg"], N_mean=der["N_genome_unweighted"]["weighted_mean_codons"])],
        pick("O1_phi_sub", "O2_sigma_DnaK", "O3_legacy_headroom_P"))
    # codon weights: every aggregate under both weightings; the headline range uses
    # the aggregate that moves most, and the per-aggregate O1 ratios are kept
    ratios = {k: phi_sub(eg[k], exp["Lw"]) / phi_sub(et[k], exp["Lw"]) for k in et}
    kmax = max(ratios, key=ratios.get)
    add("X_CODON_GENOMIC", "EXPONENTIAL", ["codon_usage_weights"], [et[kmax], eg[kmax]],
        f"translated vs genomic codon weights for the aggregate that moves most ({kmax}); "
        "O1 ratio genomic/translated for every aggregate in 'O1_ratio_by_aggregate'",
        with_(exp, f=et[kmax]), "f", [et[kmax], eg[kmax]], pick("O1_phi_sub", "O2_sigma_DnaK"),
        O1_ratio_by_aggregate=ratios,
        f_ratio_by_aggregate={k: eg[k] / et[k] for k in et})
    sp_ = json.loads((ROOT / "WO-05/wo05_results.json").read_text())["G5.5"]["per_dataset_spread"]
    add("X_F_ETEL", "EXPONENTIAL", ["substitution_freq_standing_f"], [ft[0], ft[-1]],
        "min..max of the four eTEL aggregates (translated weights); all MIXED phase. every "
        "value is an MS LOWER bound, so the range is OPEN ABOVE", exp, "f", grid(ft[0], ft[-1]),
        pick("O1_phi_sub", "O2_sigma_DnaK", "O3_legacy_headroom_P"), open_above=True,
        secondary_per_dataset_span={
            "span": [sp_["min"], sp_["max"]], "basis": "80 per-dataset usage-weighted rates "
            "(genomic weights, WO-05)",
            "outputs": evaluate(exp, "f", grid(sp_["min"], sp_["max"]),
                                pick("O1_phi_sub", "O2_sigma_DnaK"))})
    add("X_KCAT", "EXPONENTIAL", ["DnaK_cycle_rate_k_cat"], [0.003, 1.0],
        "every in-vitro DnaK rate in saved records (L02 binding 3e-3..8.4e-2; L02b T->R 0.04, "
        "R->T 1.0); NOT a bound on the in-vivo rate", exp, "k_cat", grid(0.003, 1.0, 0.04),
        pick("O2_sigma_DnaK", "O3_legacy_headroom_P"))
    add("X_K", "EXPONENTIAL", ["DnaK_affinity_K"], [0.06, 107.0],
        "R-state 0.06..2 through T-state 2.2..107 uM (in vitro, peptides)", exp, "K",
        grid(0.06, 107.0, 1.0), pick("O6_DnaK_bound_fraction", "O3_legacy_headroom_P"))

    # ---- STATIONARY
    st_core = pick("O2s_sP_breakeven")
    add("S_MU", "STATIONARY", ["growth_rate_mu"], [-0.013 / 3600, -0.007 / 3600],
        "Schmidt stationary net growth -0.01 +- 0.003 /h: negative throughout, so the "
        "dilution sink max(mu, 0) is zero over the whole span", stat, "mu",
        grid(-0.013 / 3600, -0.007 / 3600), pick("O5_dilution_share"))
    pool_params = [p for p, *_ in bb.MACHINES]
    for param, m, g_, _ in bb.MACHINES:
        s1, s3 = per[bb.STAT_COND][f"{g_}_uM_oligomer"], per[bb.STAT_COND3][f"{g_}_uM_oligomer"]
        ratio = per["LB"][f"{g_}_uM_oligomer"] / per[bb.EXP_COND][f"{g_}_uM_oligomer"]
        r = max(ratio, 1 / ratio)
        lo, hi = min(s1, s3) / r, max(s1, s3) * r
        outs = {f"O4_chains_per_{m}": chains_per(m), **(st_core if m == "DnaK" else {})}
        add(f"S_MEDIUM_{m}", "STATIONARY", [param], [lo, hi],
            f"Schmidt 1 d/3 d ({s1:.4g}, {s3:.4g}) widened by the exponential LB:glucose "
            f"factor {r:.3g} for this machine. that factor confounds medium with growth rate "
            "(1.9 vs 0.58 /h); the medium effect in stationary phase is unmeasured, so the "
            "scale is BORROWED from exponential growth and labelled as such", stat, "pools",
            [dict(stat["pools"], **{m: v}) for v in grid(lo, hi, s1)], outs,
            span_kind="SENSITIVITY")
    t1, t3 = (per[c]["total_protein_mM_wholecell"] * 1e3 for c in (bb.STAT_COND, bb.STAT_COND3))
    rp = per["LB"]["total_protein_mM_wholecell"] / per[bb.EXP_COND]["total_protein_mM_wholecell"]
    add("S_MEDIUM_P_T", "STATIONARY", ["total_protein_P_T"], [min(t1, t3) / rp, max(t1, t3) * rp],
        f"1 d/3 d widened by the exponential LB:glucose factor {rp:.3g}; that factor is an "
        "artefact of Schmidt's constant-concentration normalization (composition only), "
        "SENSITIVITY with a borrowed scale", stat, "P_T",
        grid(min(t1, t3) / rp, max(t1, t3) * rp, t1), {"O4_chains_per_DnaK": chains_per("DnaK")},
        span_kind="SENSITIVITY")
    lo_c, hi_c = ILLUSTRATIVE["norm_factor"]
    scaled = [with_(stat, P_T=stat["P_T"] * c_, pools={k: v * c_ for k, v in stat["pools"].items()})
              for c_ in grid(lo_c, hi_c, 1.0)]
    add("S_NORM", "STATIONARY", ["total_protein_P_T"] + pool_params, [lo_c, hi_c],
        "common factor on the stationary volumetric normalization (Schmidt used the "
        "glucose-exponential value); factor range ILLUSTRATIVE", stat, "_ctx", scaled,
        {**o4, **st_core}, span_kind="SENSITIVITY")
    for o in E["S_NORM"]["outputs"].values():
        o["illustrative_held"] = o["illustrative_held"] + ["normalization factor 0.5..2"]
    add("S_N", "STATIONARY", ["copy_weighted_length_N"],
        ["copy-weighted 1 d (Schmidt)", "genome-unweighted"], "two MISMATCHED weightings; "
        "the anchor culture's proteome is UNMEASURED", stat, "_ctx",
        [stat, with_(stat, Lw=Ls_["Lg"])], pick("O1_phi_sub", "O2s_sP_breakeven"))
    add("S_KCAT", "STATIONARY", ["DnaK_cycle_rate_k_cat"], [0.003, 1.0],
        "every in-vitro DnaK rate in saved records", stat, "k_cat", grid(0.003, 1.0, 0.04),
        st_core)
    add("S_K", "STATIONARY", ["DnaK_affinity_K"], [0.06, 107.0],
        "R-state through T-state in-vitro K_d", stat, "K", grid(0.06, 107.0, 1.0),
        pick("O6_DnaK_bound_fraction"))

    # ---- reference outputs per bundle, each with its honest kind
    REF = {}
    for name, c, keys in (("EXPONENTIAL", exp, ("O1_phi_sub", "O2_sigma_DnaK", "O2s_sP_breakeven",
                                                "O3_legacy_headroom_P", "O5_dilution_share")),
                          ("STATIONARY", stat, ("O1_phi_sub", "O2s_sP_breakeven",
                                                "O5_dilution_share"))):
        REF[name] = {k: {"value": OUT[k](c), "kind": reference_kind(k, name)} for k in keys}
        REF[name].update({f"O4_chains_per_{m}": {"value": c["P_T"] / v,
                                                 "kind": reference_kind("O4", name)}
                          for m, v in c["pools"].items()})
    REF["EXPONENTIAL"]["s_P_uM_per_s"] = {"value": exp["mu"] * exp["P_T"],
                                          "kind": "DERIVED_OUTPUT"}
    REF["STATIONARY"]["O1_phi_sub_at_f_pm_2SE"] = {
        "value": [phi_sub(1.82e-3 + s * 2 * 5.92e-5, stat["Lw"]) for s in (-1, 1)],
        "kind": "SENSITIVITY"}

    # ---- counterfactuals: what the prohibited borrowings would do (never bundle rows)
    CF = {}
    f_st = 1.82e-3
    CF["CF_STIKELEATHER_INTO_EXP"] = {
        "what": "stationary-harvest f (Stikeleather) placed into the EXPONENTIAL bundle",
        "at_stikeleather": {k: OUT[k](with_(exp, f=f_st))
                            for k in ("O1_phi_sub", "O2_sigma_DnaK", "O3_legacy_headroom_P")},
        "etel_span": {k: [E["X_F_ETEL"]["outputs"][k]["lo"], E["X_F_ETEL"]["outputs"][k]["hi"]]
                      for k in ("O1_phi_sub", "O2_sigma_DnaK", "O3_legacy_headroom_P")},
        "f_ratio_to_etel": [f_st / ft[-1], f_st / ft[0]]}
    CF["CF_EXP_POOLS_INTO_STAT"] = {
        "what": "exponential glucose DnaK pool placed into the STATIONARY bundle",
        "O2s_sP_breakeven": [OUT["O2s_sP_breakeven"](with_(stat, pools=exp["pools"])),
                             OUT["O2s_sP_breakeven"](stat)],
        "O4_chains_per_DnaK": [stat["P_T"] / exp["pools"]["DnaK"], stat["P_T"] / stat["DnaK"]]}
    ad = json.loads((ROOT / "WO-06/audit_derived.json").read_text())["pools"]["Glucose"]
    sums = {"legacy_50": 50.0, "DnaK+GroEL_protomers": ad["dnaK_plus_groL_protomers_uM"],
            "DnaK+GroEL14": ad["dnaK_plus_groEL14_uM"], "DnaK_only": exp["DnaK"]}
    CF["CF_SUMMED_POOL"] = {
        "what": "machines summed into one pool (legacy C_tot) instead of DnaK alone",
        "O2_sigma_DnaK": {k: sigma(with_(exp, DnaK=v)) for k, v in sums.items()},
        "O3_legacy_headroom_P": {k: legacy_headroom(with_(exp, DnaK=v)) for k, v in sums.items()}}
    leg = with_(exp, P_T=300.0, DnaK=50.0, mu=LN2 / 3600, N_mean=300.0, k_cat=1e-2)
    CF["CF_LEGACY_VALUES"] = {
        "what": "legacy P_T 300 and C_tot 50 vs bundle values (all else legacy for O3)",
        "O4_chains_per_pool": {"legacy_300_over_50": 300.0 / 50.0,
                               "bundle_P_T_over_DnaK": exp["P_T"] / exp["DnaK"]},
        "O3_legacy_headroom_P": {
            "as_published_params_zero_inclusive_genomic_f": legacy_headroom(
                with_(leg, f=eg["zero_inclusive_mean"])),
            "only_P_T_to_bundle": legacy_headroom(with_(leg, f=eg["zero_inclusive_mean"],
                                                        P_T=exp["P_T"])),
            "all_matched_values_substituted": legacy_headroom(exp)},
        "control": "as_published value must equal WO-05 zero_inclusive_input_both (63.10)"}
    return {"illustrative": ILLUSTRATIVE, "references": refs, "reference_outputs": REF,
            "effects": E, "counterfactuals": CF,
            "label": "every effect range varies a MISMATCHED input and is a SENSITIVITY or "
                     "LEGACY_CONDITIONAL range, not a prediction"}


if __name__ == "__main__":
    r = run()
    (HERE / "mismatch_effects.json").write_text(json.dumps(r, indent=1, default=float))
    for k, v in r["effects"].items():
        print(k, {o: (round(x["lo"], 4), round(x["hi"], 4)) for o, x in v["outputs"].items()})
    print(json.dumps(r["reference_outputs"], indent=0, default=float))
    print(json.dumps(r["counterfactuals"], indent=0, default=float))
