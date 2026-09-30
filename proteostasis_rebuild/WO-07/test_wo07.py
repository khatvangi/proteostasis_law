"""WO-07 tests. run from proteostasis_rebuild/:

    PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s WO-07 -p 'test_wo07.py' -v

every negative control mutates a copy of the real table and must be caught; the
positive controls show the rules do not simply reject everything.
"""
import copy
import json
import math
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_bundles as bb  # noqa: E402
import check_bundles as cb  # noqa: E402

ROWS = cb.load()
SRC = cb.sources()
EFF = json.loads((HERE / "mismatch_effects.json").read_text())


def find(rows, bundle, param, flag=None, candidate=None):
    return next(r for r in rows if r["bundle"] == bundle and r["param"] == param
                and (flag is None or r["match_flag"] == flag)
                and (candidate is None or r["candidate"] == candidate))


def errs(rows):
    return cb.validate(rows, SRC)


class RealTable(unittest.TestCase):
    def test_passes(self):
        self.assertEqual(errs(ROWS), [])

    def test_regenerates_identically(self):
        self.assertTrue(cb.regenerated_equal(ROWS))

    def test_g71_both_bundles_and_flags(self):
        self.assertEqual({r["bundle"] for r in ROWS}, {bb.EXP, bb.STAT})
        for r in ROWS:
            self.assertIn(r["match_flag"], cb.FLAGS)
            self.assertTrue(r["condition_tag"].startswith("phase="))

    def test_every_param_accounted_in_both(self):
        for b in (bb.EXP, bb.STAT):
            self.assertEqual({r["param"] for r in ROWS if r["bundle"] == b}, set(bb.PARAMS))


class G72OppositePhase(unittest.TestCase):
    """the required test: an opposite-phase value without MISMATCHED fails."""

    def inject(self, src_bundle, dst_bundle, param, flag, value_phase=None, axes=None):
        rows = copy.deepcopy(ROWS)
        donor = copy.deepcopy(next(r for r in rows if r["bundle"] == src_bundle
                                   and r["param"] == param and r["value_phase"] == src_bundle))
        donor["bundle"], donor["match_flag"] = dst_bundle, flag
        donor["candidate"] = "injected"
        if value_phase:
            donor["value_phase"] = value_phase
        if axes:
            donor["mismatch_axes"] = axes
        rows.append(donor)
        return errs(rows)

    def test_stationary_value_in_exp_marked_matched_fails(self):
        e = self.inject(bb.STAT, bb.EXP, "pool_DnaK", "MATCHED")
        self.assertTrue(any("opposite-phase value (STATIONARY) in EXPONENTIAL" in x for x in e))

    def test_relabelled_phase_is_still_caught(self):
        # the lie 'value_phase=EXPONENTIAL' does not help: phase comes from the source
        e = self.inject(bb.STAT, bb.EXP, "pool_DnaK", "MATCHED", value_phase=bb.EXP)
        self.assertTrue(any("but the source records STATIONARY" in x for x in e))
        self.assertTrue(any("opposite-phase value" in x for x in e))

    def test_exponential_value_in_stat_unflagged_fails(self):
        for flag in ("MATCHED", "UNMEASURED"):
            e = self.inject(bb.EXP, bb.STAT, "total_protein_P_T", flag)
            self.assertTrue(any("opposite-phase value (EXPONENTIAL) in STATIONARY" in x
                                for x in e), flag)

    def test_flagged_without_phase_axis_fails(self):
        e = self.inject(bb.STAT, bb.EXP, "pool_ClpB6", "MISMATCHED", axes="STRAIN")
        self.assertTrue(any("lacks PHASE_OPPOSITE" in x for x in e))

    def test_positive_control_flagged_with_effect_passes_phase_rule(self):
        # properly flagged: the phase rules are silent (only the effect-bundle and
        # coverage bookkeeping may complain, because no effect was computed for it)
        e = self.inject(bb.STAT, bb.EXP, "pool_DnaK", "MISMATCHED", axes="PHASE_OPPOSITE")
        self.assertFalse(any("opposite-phase" in x for x in e))

    def test_mixed_and_invitro_must_be_mismatched(self):
        rows = copy.deepcopy(ROWS)
        find(rows, bb.EXP, "substitution_freq_standing_f", candidate="eTEL_zero_inclusive_mean")[
            "match_flag"] = "MATCHED"
        find(rows, bb.STAT, "DnaK_cycle_rate_k_cat", candidate="T_to_R_step")[
            "match_flag"] = "MATCHED"
        e = errs(rows)
        self.assertTrue(any("MIXED value not flagged" in x for x in e))
        self.assertTrue(any("IN_VITRO value not flagged" in x for x in e))


class Separation(unittest.TestCase):
    def test_stikeleather_never_in_exponential(self):
        rows = copy.deepcopy(ROWS)
        r = copy.deepcopy(find(rows, bb.STAT, "substitution_freq_standing_f", "MATCHED"))
        r.update(bundle=bb.EXP, match_flag="MISMATCHED", mismatch_axes="PHASE_OPPOSITE",
                 effect_ids="X_F_ETEL", candidate="injected")
        rows.append(r)
        self.assertTrue(any("Stikeleather stationary error input in EXPONENTIAL" in x
                            for x in errs(rows)))

    def test_stikeleather_is_stationary_error_input(self):
        r = find(ROWS, bb.STAT, "substitution_freq_standing_f", "MATCHED")
        self.assertEqual(r["source_key"], "audit:E03")
        self.assertAlmostEqual(float(r["value"]), 1.82e-3)
        self.assertFalse(any("42406629" in x["source"] for x in ROWS if x["bundle"] == bb.EXP))

    def test_summed_pool_rejected(self):
        rows = copy.deepcopy(ROWS)
        r = copy.deepcopy(find(rows, bb.EXP, "pool_DnaK", "MATCHED"))
        r.update(param="pool_DnaK+GroEL+ClpB", machine="DnaK+GroEL+ClpB")
        rows.append(r)
        self.assertTrue(any("summed chaperone pool" in x for x in errs(rows)))

    def test_pool_row_names_one_machine(self):
        rows = copy.deepcopy(ROWS)
        find(rows, bb.EXP, "pool_GroEL14", "MATCHED")["machine"] = "chaperones"
        self.assertTrue(any("exactly one machine" in x for x in errs(rows)))

    def test_machines_separate_in_both_bundles(self):
        for b in (bb.EXP, bb.STAT):
            ms = {r["machine"] for r in ROWS if r["bundle"] == b and r["param"].startswith("pool_")}
            self.assertTrue(cb.REQUIRED_MACHINES <= ms)


class Bookkeeping(unittest.TestCase):
    def test_unmeasured_cannot_carry_value(self):
        rows = copy.deepcopy(ROWS)
        find(rows, bb.STAT, "synthesis_rate_s_P", "UNMEASURED")["value"] = "0.478"
        self.assertTrue(any("UNMEASURED carries a value" in x for x in errs(rows)))

    def test_mismatched_needs_effect(self):
        rows = copy.deepcopy(ROWS)
        find(rows, bb.EXP, "total_protein_P_T", candidate="BNID_104726")["effect_ids"] = "X_NOPE"
        self.assertTrue(any("effect X_NOPE not computed" in x for x in errs(rows)))

    def test_effect_bundle_must_match(self):
        rows = copy.deepcopy(ROWS)
        find(rows, bb.EXP, "total_protein_P_T", candidate="BNID_104726")["effect_ids"] = "S_NORM"
        self.assertTrue(any("belongs to STATIONARY" in x for x in errs(rows)))

    def test_missing_matched_or_unmeasured_fails(self):
        rows = [r for r in ROWS if not (r["bundle"] == bb.STAT and r["param"] == "pool_DnaK"
                                        and r["match_flag"] == "UNMEASURED")]
        self.assertTrue(any("STATIONARY/pool_DnaK: needs exactly one" in x for x in errs(rows)))

    def test_matched_must_hit_anchor(self):
        rows = copy.deepcopy(ROWS)
        find(rows, bb.EXP, "pool_DnaK", "MATCHED")["medium"] = "LB"
        self.assertTrue(any("!= anchor" in x for x in errs(rows)))

    def test_required_fields(self):
        for col in ("strain", "medium", "temperature", "source", "measurement_type",
                    "condition_tag", "match_flag", "value_phase"):
            rows = copy.deepcopy(ROWS)
            rows[0][col] = ""
            self.assertTrue(any(f"empty {col}" in x for x in errs(rows)), col)

    def test_value_must_equal_source(self):
        rows = copy.deepcopy(ROWS)
        find(rows, bb.EXP, "pool_DnaK", "MATCHED")["value"] = "50"
        self.assertTrue(any("not found in source" in x for x in errs(rows)))
        rows = copy.deepcopy(ROWS)
        find(rows, bb.EXP, "DnaK_cycle_rate_k_cat", candidate="T_to_R_step")["value"] = "0.4"
        self.assertTrue(any("not found in source" in x for x in errs(rows)))


class ReviewerBypasses(unittest.TestCase):
    """each case passed the first validator (independent review, 2026-09-30)."""

    def test_self_reported_conditions_not_trusted(self):
        # yeast k_d typed as E. coli BW25113 VERIFIED
        rows = copy.deepcopy(ROWS)
        r = find(rows, bb.EXP, "misfolded_degradation_k_d", "UNMEASURED")
        r.update(value="0.0003", source="Belle 2006", source_key="audit:L01a", value_phase=bb.EXP,
                 organism="Escherichia coli", strain="BW25113", medium="M9 glucose",
                 temperature="37 C", wo06_status="VERIFIED", match_flag="MATCHED",
                 measurement_type="MEASURED")
        self.assertTrue(any("SOURCE conditions" in x or "the source is" in x for x in errs(rows)))

    def test_bnid_cannot_be_primary(self):
        rows = copy.deepcopy(ROWS)
        r = find(rows, bb.EXP, "total_protein_P_T", "MATCHED")
        r.update(value="4000", source_key="audit:L06b")
        self.assertTrue(any("SOURCE conditions" in x for x in errs(rows)))

    def test_status_retyping_caught(self):
        rows = copy.deepcopy(ROWS)
        r = find(rows, bb.EXP, "codon_usage_weights", candidate="genomic")
        r.update(match_flag="MATCHED", wo06_status="VERIFIED", mismatch_axes="-",
                 effect_ids="-", value_phase="PHASE_INVARIANT")
        rows = [x for x in rows if not (x["bundle"] == bb.EXP and x["param"] ==
                                         "codon_usage_weights" and x["candidate"] == "primary")]
        self.assertTrue(any("MISMATCHED_CONDITION" in x for x in errs(rows)))

    def test_vector_must_match_source_key(self):
        rows = copy.deepcopy(ROWS)
        find(rows, bb.EXP, "codon_usage_weights", "MATCHED")["value"] = \
            "vector:derived.json#N_abundance_weighted_stationary_1d"
        self.assertTrue(any("not found in source" in x for x in errs(rows)))

    def test_parser_hyphenated_phase_is_ambiguous(self):
        self.assertEqual(cb.audit_phase("post-exponential (overnight)"), "AMBIGUOUS")
        self.assertEqual(cb.audit_phase("stationary; genome-wide proteomics"), "STATIONARY")
        self.assertEqual(cb.audit_phase("unmixed stationary"), "STATIONARY")
        self.assertEqual(cb.audit_phase("stationary (overnight)"), "STATIONARY")
        self.assertEqual(cb.audit_phase("MIXED (80 PRIDE datasets)"), "MIXED")

    def test_effect_must_vary_the_row_param(self):
        rows = copy.deepcopy(ROWS)
        find(rows, bb.STAT, "pool_ClpB6", candidate="Schmidt_1d")["effect_ids"] = "S_MU"
        self.assertTrue(any("varies ['growth_rate_mu']" in x for x in errs(rows)))

    def test_mM_scaling_only_for_concentrations(self):
        rows = copy.deepcopy(ROWS)
        find(rows, bb.EXP, "DnaK_cycle_rate_k_cat", candidate="T_to_R_step")["value"] = "40"
        self.assertTrue(any("not found in source" in x for x in errs(rows)))


class G73Effects(unittest.TestCase):
    def test_every_mismatch_has_finite_range(self):
        for r in ROWS:
            if r["match_flag"] != "MISMATCHED":
                continue
            for eid in r["effect_ids"].split(";"):
                e = EFF["effects"][eid]
                self.assertTrue(e["outputs"], eid)
                for o in e["outputs"].values():
                    self.assertTrue(math.isfinite(o["lo"]) and o["lo"] <= o["hi"], eid)

    def test_sensitivity_labels(self):
        # every effect varies a mismatched input: never DERIVED_OUTPUT
        for eid, e in EFF["effects"].items():
            for name, o in e["outputs"].items():
                self.assertIn(o["kind"], ("SENSITIVITY", "LEGACY_CONDITIONAL"), f"{eid}/{name}")

    def test_reference_kinds_honest(self):
        ref = EFF["reference_outputs"]
        self.assertEqual(ref["EXPONENTIAL"]["O4_chains_per_DnaK"]["kind"], "DERIVED_OUTPUT")
        self.assertEqual(ref["EXPONENTIAL"]["O1_phi_sub"]["kind"], "SENSITIVITY")  # f is MIXED
        self.assertEqual(ref["STATIONARY"]["O1_phi_sub"]["kind"], "SENSITIVITY")   # N mismatched
        self.assertEqual(ref["STATIONARY"]["O4_chains_per_DnaK"]["kind"], "SENSITIVITY")

    def test_negative_growth_is_no_dilution(self):
        o = EFF["effects"]["S_MU"]["outputs"]["O5_dilution_share"]
        self.assertEqual((o["lo"], o["hi"]), (0.0, 0.0))
        self.assertEqual(EFF["reference_outputs"]["STATIONARY"]["O5_dilution_share"]["value"], 0.0)

    def test_mislabelled_effect_caught(self):
        import copy as c
        bad = c.deepcopy(SRC)
        bad[3]["effects"]["S_N"]["outputs"]["O1_phi_sub"]["kind"] = "DERIVED_OUTPUT"
        self.assertTrue(any("labelled DERIVED_OUTPUT" in x for x in cb.validate(ROWS, bad)))

    def test_legacy_control_reproduces_wo05(self):
        w5 = json.loads((HERE.parent / "WO-05/wo05_results.json").read_text())
        ref = w5["G5.5"]["effect_on_old_x"]["zero_inclusive_input_both"]
        got = EFF["counterfactuals"]["CF_LEGACY_VALUES"]["O3_legacy_headroom_P"][
            "as_published_params_zero_inclusive_genomic_f"]
        self.assertAlmostEqual(got, ref, places=6)

    def test_derived_genomic_etel_reproduces_wo05(self):
        w5 = json.loads((HERE.parent / "WO-05/wo05_results.json").read_text())["G5.5"]["inputs"]
        der = json.loads((HERE / "derived.json").read_text())["etel_genomic_weights"]
        for k in der:
            self.assertAlmostEqual(der[k], w5[k]["f"], places=15)

    def test_o4_is_normalization_invariant(self):
        o = EFF["effects"]["S_NORM"]["outputs"]["O4_chains_per_DnaK"]
        self.assertAlmostEqual(o["lo"], o["hi"], places=9)

    def test_hand_recompute(self):
        per, _, _ = bb.load_sources()
        gl = per["Glucose"]
        ref = EFF["reference_outputs"]["EXPONENTIAL"]
        self.assertAlmostEqual(ref["O4_chains_per_DnaK"]["value"],
                               gl["total_protein_mM_wholecell"] * 1e3 / gl["dnaK_uM_oligomer"])
        sP = gl["growth_per_h"] / 3600 * gl["total_protein_mM_wholecell"] * 1e3
        self.assertAlmostEqual(ref["s_P_uM_per_s"]["value"], sP)


if __name__ == "__main__":
    unittest.main()
