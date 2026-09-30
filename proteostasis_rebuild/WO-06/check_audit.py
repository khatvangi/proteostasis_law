#!/usr/bin/env python3
"""WO-06 automated check (G6.1-G6.3) for parameter_audit.tsv.

  python check_audit.py            # exit 0 iff no errors; prints every error

validate(rows) returns a list of error strings and is what test_wo06 calls, on
the real table and on deliberately broken copies (negative controls).

rules:
  R1 required fields non-empty (incl. organism, strain, phase, medium,
     temperature, unit, source, measurement type, verification status)
  R2 closed vocabularies
  R3 NO_SOURCE in a condition field only for ASSUMED / ILLUSTRATIVE rows
  R4 ASSUMED / ILLUSTRATIVE measurement type => same verification status
     (UNKNOWN = a citation exists but cannot be matched, so what it measured is unknown)
     (a placeholder can never be VERIFIED)
  R5 VERIFIED => record MATCHED or LOCAL_DATA, value located, org/cond MATCHED,
     and the record resolves to a saved file
  R5b MISMATCHED_CONDITION => value located
  R6 every PMID cited has a saved PubMed record whose journal|year|volume equals
     record_citation (the bibliographic check, from the record, not memory)
  R7 MISCITED status => miscitation_evidence; citation_check MISCITED/UNMATCHED
     with a citmatch key => the saved ecitmatch result is NOT_FOUND or a PMID
     listed in source_record as what it really resolves to; its CONTROL fired
  R8 citation_check MATCHED with a citmatch key => ecitmatch returned the row's PMID
  R9 citation_check UNMATCHED => status UNVERIFIED; citation_check MISCITED => status MISCITED
"""
import csv
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REC = HERE / "records"
BASE = Path("/storage/kiran-stuff")
TSV = HERE / "parameter_audit.tsv"

REQUIRED = ["row_id", "param", "model_scope", "used_at", "value_used", "unit", "model_role",
            "cited_as", "source_record", "organism", "strain", "phase", "medium", "temperature",
            "measured_in_source", "measurement_type", "citation_check", "value_located",
            "org_cond_match", "verification_status", "what_source_supports"]
CONDITION = ["organism", "strain", "phase", "medium", "temperature"]
VOCAB = {
    "measurement_type": {"MEASURED", "DERIVED_FROM_MEASUREMENT", "ESTIMATED", "MODEL_FIT",
                         "THEORETICAL", "ASSUMED", "ILLUSTRATIVE", "UNKNOWN"},
    "citation_check": {"MATCHED", "MISCITED", "AMBIGUOUS", "UNMATCHED", "NO_CITATION", "LOCAL_DATA"},
    "value_located": {"YES", "NO", "NA"},
    "verification_status": {"VERIFIED", "MISMATCHED_CONDITION", "NOT_SUPPORTED", "MISCITED",
                            "UNVERIFIED", "ASSUMED", "ILLUSTRATIVE"},
}
OCM = {"MATCHED", "ORGANISM_MISMATCH", "CONDITION_MISMATCH", "QUANTITY_MISMATCH", "NA"}
PLACEHOLDER = {"ASSUMED", "ILLUSTRATIVE"}
# which control must fire for each probe expected to fail
CONTROLS = {
    "Ciryam2013_PNAS_110_E3453": "CONTROL_Ciryam2013_CellRep_5",
    "Bednarska2013_MolCell_52_617": "CONTROL_Bednarska2013_Microbiology_159",
    "Stirling2018_CellRep_25_2242": "CONTROL_CellRep_2018_25_page",
    "DrummondWilke2009_Cell": "CONTROL_DrummondWilke2008_Cell",
    "Yamanaka2017_CurrBiol": "CONTROL_CurrBiol_2017_author",
    "Pierpaoli1997_EMBOJ": "CONTROL_Pierpaoli1997_JMolBiol",
}


def load(path=TSV):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def pubmed_fields(pmid):
    """journal ISO abbreviation, year, volume from the saved efetch record."""
    f = REC / f"pubmed_{pmid}.xml"
    if not f.exists():
        return None
    x = f.read_text().split("<ReferenceList")[0]
    j = re.search(r"<ISOAbbreviation>(.*?)</ISOAbbreviation>", x)
    v = re.search(r"<Volume>(.*?)</Volume>", x)
    y = re.search(r"<PubDate>.*?<Year>(\d{4})</Year>", x, re.S) or \
        re.search(r"<PubDate>.*?<MedlineDate>(\d{4})", x, re.S)
    return (j.group(1) if j else "", y.group(1) if y else "", v.group(1) if v else "")


def pmids(source_record):
    return re.findall(r"PMID:(\d+)", source_record)


def record_exists(r):
    s = r["source_record"]
    if s.startswith("PMID:"):
        return all((REC / f"pubmed_{p}.xml").exists() for p in pmids(s))
    if s.startswith("LOCAL:"):
        return (BASE / s[len("LOCAL:"):]).exists()
    if s.startswith("BNID:"):
        return (REC / f"bionumbers_{s.split(':')[1].split()[0]}.html").exists()
    return False


def validate(rows, citmatch=None):
    if citmatch is None:
        citmatch = json.loads((REC / "citmatch.json").read_text())
    err = []
    ids = [r.get("row_id", "") for r in rows]
    if len(set(ids)) != len(ids):
        err.append("duplicate row_id")
    for r in rows:
        rid = r.get("row_id") or "?"
        # R1
        for c in REQUIRED:
            if not (r.get(c) or "").strip():
                err.append(f"{rid}: empty {c}")
        # R2
        for c, voc in VOCAB.items():
            if r.get(c) not in voc:
                err.append(f"{rid}: {c}={r.get(c)!r} not in vocabulary")
        if not set((r.get("org_cond_match") or "").split("+")) <= OCM:
            err.append(f"{rid}: org_cond_match={r.get('org_cond_match')!r} not in vocabulary")
        st, mt = r.get("verification_status"), r.get("measurement_type")
        # R3
        if any(r.get(c) == "NO_SOURCE" for c in CONDITION) and st not in PLACEHOLDER:
            err.append(f"{rid}: NO_SOURCE condition on a {st} row")
        # R4
        if mt in PLACEHOLDER and st not in PLACEHOLDER:
            err.append(f"{rid}: measurement_type {mt} but status {st}")
        # R5
        if st == "VERIFIED":
            if r.get("citation_check") not in {"MATCHED", "LOCAL_DATA"}:
                err.append(f"{rid}: VERIFIED without a matched record")
            if r.get("value_located") != "YES":
                err.append(f"{rid}: VERIFIED but value not located")
            if r.get("org_cond_match") != "MATCHED":
                err.append(f"{rid}: VERIFIED with org/cond {r.get('org_cond_match')}")
            if not record_exists(r):
                err.append(f"{rid}: VERIFIED but record file missing for {r.get('source_record')}")
        # R5b a condition mismatch is only claimed for a value actually seen
        if st == "MISMATCHED_CONDITION" and r.get("value_located") != "YES":
            err.append(f"{rid}: MISMATCHED_CONDITION but value not located (use UNVERIFIED)")
        # R6
        cits = (r.get("record_citation") or "").split(";")
        ps = pmids(r.get("source_record", ""))
        if ps and len(cits) != len(ps):
            err.append(f"{rid}: {len(ps)} PMIDs but {len(cits)} record_citation entries")
        for p, c in zip(ps, cits):
            got = pubmed_fields(p)
            if got is None:
                err.append(f"{rid}: no saved record for PMID {p}")
            elif "|".join(got) != c.strip():
                err.append(f"{rid}: PMID {p} record is {'|'.join(got)!r}, row says {c.strip()!r}")
        # R7 / R8
        # R9 (G6.2 wording): unmatched -> UNVERIFIED, mismatched -> MISCITED
        if r.get("citation_check") == "UNMATCHED" and st != "UNVERIFIED":
            err.append(f"{rid}: UNMATCHED citation but status {st}")
        if r.get("citation_check") == "MISCITED" and st != "MISCITED":
            err.append(f"{rid}: MISCITED citation but status {st}")
        key = r.get("citmatch_key", "-")
        if st == "MISCITED" and not (r.get("miscitation_evidence") or "").strip("- "):
            err.append(f"{rid}: MISCITED without evidence")
        if key not in ("-", ""):
            if key not in citmatch:
                err.append(f"{rid}: citmatch key {key} not in citmatch.json")
                continue
            res = citmatch[key]["result"]
            if r.get("citation_check") in {"MISCITED", "UNMATCHED"}:
                if not (res == "NOT_FOUND" or res in ps):
                    err.append(f"{rid}: citmatch {key} -> {res}, not NOT_FOUND nor a listed PMID")
                ctl = CONTROLS.get(key)
                if ctl is None or not citmatch.get(ctl, {}).get("result", "").isdigit():
                    err.append(f"{rid}: no fired control for probe {key}")
            if r.get("citation_check") == "MATCHED" and res not in ps:
                err.append(f"{rid}: citmatch {key} -> {res} but row cites {ps}")
    return err


def main():
    err = validate(load())
    for e in err:
        print("FAIL", e)
    print(f"{len(load())} rows, {len(err)} errors")
    return 1 if err else 0


if __name__ == "__main__":
    sys.exit(main())
