"""WO-06 tests. run from proteostasis_rebuild/:

    PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s WO-06 -p 'test_wo06.py' -v

no network: everything is checked against records/ saved by fetch_records.py.
"""
import copy
import csv
import hashlib
import json
import re
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_audit as ca  # noqa: E402

BASE = Path("/storage/kiran-stuff")
ROWS = ca.load()


def params_named():
    names = set()
    for r in ROWS:
        for tok in re.split(r"[;,()\s]+", r["param"]):
            if tok:
                names.add(tok)
    return names


class G61G63Fields(unittest.TestCase):
    def test_real_table_passes(self):
        self.assertEqual(ca.validate(ROWS), [])

    def test_missing_organism_fails(self):
        for col in ("organism", "strain", "phase", "medium", "temperature", "unit",
                    "source_record", "measurement_type", "verification_status"):
            bad = copy.deepcopy(ROWS)
            bad[0][col] = ""
            self.assertTrue(any(f"empty {col}" in e for e in ca.validate(bad)), col)

    def test_placeholder_cannot_be_verified(self):
        bad = copy.deepcopy(ROWS)
        r = next(x for x in bad if x["measurement_type"] == "ASSUMED")
        r["verification_status"] = "VERIFIED"
        self.assertTrue(ca.validate(bad))

    def test_no_source_only_on_placeholders(self):
        bad = copy.deepcopy(ROWS)
        r = next(x for x in bad if x["verification_status"] == "MISMATCHED_CONDITION")
        r["organism"] = "NO_SOURCE"
        self.assertTrue(any("NO_SOURCE" in e for e in ca.validate(bad)))

    def test_verified_needs_record_and_match(self):
        bad = copy.deepcopy(ROWS)
        r = next(x for x in bad if x["verification_status"] == "MISMATCHED_CONDITION")
        r["verification_status"] = "VERIFIED"      # org/cond not MATCHED
        self.assertTrue(ca.validate(bad))
        bad = copy.deepcopy(ROWS)
        r = next(x for x in bad if x["row_id"] == "L06b")
        r["source_record"] = "BNID:999999"         # record file absent
        self.assertTrue(any("record file missing" in e for e in ca.validate(bad)))

    def test_g62_wording_enforced(self):
        bad = copy.deepcopy(ROWS)
        next(x for x in bad if x["citation_check"] == "UNMATCHED")["verification_status"] = "ASSUMED"
        self.assertTrue(any("UNMATCHED citation" in e for e in ca.validate(bad)))
        bad = copy.deepcopy(ROWS)
        next(x for x in bad if x["citation_check"] == "MISCITED")["verification_status"] = "MISMATCHED_CONDITION"
        self.assertTrue(any("MISCITED citation" in e for e in ca.validate(bad)))

    def test_mismatch_requires_located_value(self):
        bad = copy.deepcopy(ROWS)
        r = next(x for x in bad if x["verification_status"] == "MISMATCHED_CONDITION")
        r["value_located"] = "NO"
        self.assertTrue(any("MISMATCHED_CONDITION but value not located" in e for e in ca.validate(bad)))

    def test_vocabulary_enforced(self):
        bad = copy.deepcopy(ROWS)
        bad[0]["verification_status"] = "VERIFIED_FROM_MEMORY"
        self.assertTrue(ca.validate(bad))

    def test_status_counts(self):
        st = {r["verification_status"] for r in ROWS}
        self.assertTrue({"VERIFIED", "MISCITED", "UNVERIFIED", "ASSUMED"} <= st)


class G62Bibliography(unittest.TestCase):
    def test_every_pmid_saved_and_matches(self):
        for r in ROWS:
            for p in ca.pmids(r["source_record"]):
                self.assertIsNotNone(ca.pubmed_fields(p), p)
        # the check is real: a wrong volume is caught
        bad = copy.deepcopy(ROWS)
        r = next(x for x in bad if x["row_id"] == "L11")
        r["record_citation"] = "EMBO J|1999|19"
        self.assertTrue(any("record is" in e for e in ca.validate(bad)))

    def test_fetch_list_covers_tsv(self):
        src = (HERE / "fetch_records.py").read_text()
        listed = set(re.findall(r'"(\d{6,8})",\s*#', src))
        used = {p for r in ROWS for p in ca.pmids(r["source_record"])}
        self.assertLessEqual(used, listed)

    def test_controls_fired(self):
        cm = json.loads((HERE / "records" / "citmatch.json").read_text())
        for probe, ctl in ca.CONTROLS.items():
            self.assertTrue(cm[ctl]["result"].isdigit(), ctl)
            self.assertIn(probe, cm)

    def test_known_miscitations(self):
        cm = json.loads((HERE / "records" / "citmatch.json").read_text())
        self.assertEqual(cm["Pierpaoli1997_EMBOJ"]["result"], "NOT_FOUND")
        self.assertEqual(cm["DrummondWilke2009_Cell"]["result"], "NOT_FOUND")
        self.assertEqual(cm["Ciryam2013_PNAS_110_E3453"]["result"], "NOT_FOUND")
        self.assertEqual(cm["Bednarska2013_MolCell_52_617"]["result"], "24239291")
        x = (HERE / "records" / "pubmed_24239291.xml").read_text()
        self.assertIn("sliding clamp", x)                 # unrelated to aggregation
        x = (HERE / "records" / "pubmed_25466257.xml").read_text()
        self.assertIn("Yeasts S. cerevisiae and S. pombe", x)

    def test_quoted_values_in_saved_text(self):
        """values the table says were located are present in the saved record."""
        ab = (HERE / "records" / "pubmed_9843444.xml").read_text()
        self.assertIn("Kd = 0.06-2 microM", ab)
        self.assertIn("0.003-0.084 s-1", ab)
        self.assertIn("25 degreesC", ab)
        ab = (HERE / "records" / "pubmed_11135201.xml").read_text()
        self.assertIn("first order kinetics", ab)
        t = re.sub(r"\s+", " ", re.sub("<[^>]+>", " ", (HERE / "records" / "pmc_PMC13335486.xml").read_text()))
        self.assertIn("1.82 × 10 −3 per codon, SE = 5.92 × 10 −5", t)
        self.assertIn("stationary phase", t)
        t = re.sub(r"\s+", " ", re.sub("<[^>]+>", " ", (HERE / "records" / "pmc_PMC2764353.xml").read_text()))
        self.assertIn("remains essentially unknown", t)
        bn = (HERE / "records" / "bionumbers_104726.html").read_text(errors="ignore")
        self.assertIn("4 mM", re.sub(r"\s+", " ", re.sub("<[^>]+>", " ", bn)))


class Derived(unittest.TestCase):
    def test_audit_checks_reproduce(self):
        import audit_checks
        saved = json.loads((HERE / "audit_derived.json").read_text())
        self.assertEqual(json.loads(json.dumps(audit_checks.run())), saved)

    def test_schmidt_reproduce(self):
        import schmidt_pools
        saved = json.loads((HERE / "schmidt_pools.json").read_text())
        self.assertEqual(json.loads(json.dumps(schmidt_pools.run()))["summary"], saved["summary"])

    def test_key_findings(self):
        d = json.loads((HERE / "audit_derived.json").read_text())
        self.assertGreater(d["prot_tot"]["fold_low_vs_schmidt_glucose"], 9)
        self.assertFalse(d["k_deg"]["baseline_inside_cited_consensus"])
        self.assertEqual(d["lengths"]["median"], 271.0)
        self.assertLess(d["S_code"]["all_positions_unweighted"], 0.26)
        self.assertGreater(d["S_code"]["third_position_unweighted"], 0.69)
        s = json.loads((HERE / "schmidt_pools.json").read_text())["summary"]
        self.assertLess(s["clpB_uM_oligomer"]["exp_max"], 0.1)
        self.assertLess(s["dnaK_uM_protomer"]["exp_max"], 30)


class G61Coverage(unittest.TestCase):
    """every parameter name in the legacy and rebuild code appears in the table."""

    def test_legacy_params(self):
        src = (BASE / "proteostasis-P1/two_pool_ode.py").read_text()
        body = src[src.index("class Params"):src.index("BASELINE = Params()")]
        names = set(re.findall(r"^\s+(\w+): float =", body, re.M))
        samp = src[src.index("def sample_params"):src.index("def part_D")]
        names |= set(re.findall(r"^\s+(\w+)=", samp, re.M))
        self.assertGreaterEqual(len(names), 13)
        self.assertEqual(names - params_named(), set())

    def test_arithmetic_and_paired(self):
        src = (BASE / "proteostasis-P1/arithmetic_stress_test.py").read_text()
        body = src[src.index("class Baseline"):src.index("BASE = Baseline()")]
        names = set(re.findall(r"^\s+(\w+): float =", body, re.M))
        self.assertEqual(names - params_named(), set())

    def test_rebuild_params(self):
        sys.path.insert(0, str(HERE.parent / "WO-02"))
        import model as m2
        names = set(m2.PARAM_NAMES) | set(m2.SCENARIO)
        src = (HERE.parent / "WO-04/bifurcation.py").read_text()
        samp = src[src.index("def sample"):src.index("def fold_in_eps")]
        names |= set(re.findall(r'"(\w+)": (?:lu\(|rng\.|[\d.])', samp))   # dict keys only
        src = (HERE.parent / "WO-03/run_wo03.py").read_text()
        names |= set(re.findall(r"(k_\w+)=\d", src[src.index("base = dict"):src.index("K_dT =")]))
        names |= set(re.findall(r'"(\w+)": [\d.]', src[src.index("ext = {"):src.index("# reduction test")]))
        self.assertEqual(names - params_named(), set())


class Preservation(unittest.TestCase):
    def test_legacy_hashes_unchanged(self):
        with open(HERE.parent / "WO-00" / "legacy_hashes.tsv") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                h = hashlib.sha256((BASE / r["path"]).read_bytes()).hexdigest()
                self.assertEqual(h, r["sha256"], r["path"])

    def test_extra_legacy_inputs_unchanged(self):
        # legacy file read by WO-06 but not in the WO-00 manifest
        p = BASE / "proteostasis_law/envelope-paper/data/raw/ecoli_k12_cds.fna"
        self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),
                         "9bd477a9ecd7b84d00fe6c5a90f4f1a7cf6066c43f8e943dae9eb45fef78fbb5")


if __name__ == "__main__":
    unittest.main()
