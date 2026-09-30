"""WO-07 validator for bundles.tsv (G7.1-G7.3). exit 1 on any error.

the phase of every value is re-derived from its SOURCE (the WO-06 audit row,
Schmidt's per-condition phase, or derived.json), never taken from the row's own
value_phase label. so a stationary value relabelled 'EXPONENTIAL' is caught.
"""
import csv
import json
import math
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import build_bundles as bb  # noqa: E402

FLAGS = {"MATCHED", "MISMATCHED", "UNMEASURED"}
PHASES = {"EXPONENTIAL", "STATIONARY", "MIXED", "IN_VITRO", "PHASE_INVARIANT", "AMBIGUOUS", "NONE"}
OPPOSITE = {bb.EXP: bb.STAT, bb.STAT: bb.EXP}
REQUIRED = [c for c in bb.COLUMNS if c != "machine"]
STIKELEATHER = ("PMID:42406629", "audit:E03")
MACHINE_NAMES = {m for _, m, _, _ in bb.MACHINES}
REQUIRED_MACHINES = {"DnaK", "GroEL", "ClpB"}
OUTPUT_KINDS = {"DERIVED_OUTPUT", "SENSITIVITY", "LEGACY_CONDITIONAL"}


def load(path=HERE / "bundles.tsv"):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sources():
    audit = {r["row_id"]: r for r in csv.DictReader(
        open(ROOT / "WO-06/parameter_audit.tsv", newline=""), delimiter="\t")}
    per, der, _ = bb.load_sources()
    eff = json.loads((HERE / "mismatch_effects.json").read_text())
    return audit, per, der, eff


def audit_phase(text):
    """keyword parser on WO-06's phase field. whole words only; a hyphenated
    'post-exponential' or any two-phase text is AMBIGUOUS (forces MISMATCHED)."""
    t = text.lower()
    word = lambda w: re.search(rf"(?<![-\w]){w}(?![-\w])", t) is not None
    st, ex = word("stationary"), word("exponential") or word("balanced")
    if word("mixed"):
        return "MIXED"
    if "in vitro" in t:
        return "IN_VITRO"
    if "genome" in t and not (st or ex):
        return "PHASE_INVARIANT"
    if st and not ex:
        return "STATIONARY"
    if ex and not st:
        return "EXPONENTIAL"
    return "AMBIGUOUS"


def schmidt_phase(per, cond):
    return bb.STAT if per[cond]["phase"] == "stationary" else bb.EXP


def source_phase(key, audit, per, der):
    """phase of the value as the SOURCE records it."""
    if key == "none":
        return "NONE"
    kind, _, rest = key.partition(":")
    if kind == "audit":
        return audit_phase(audit[rest]["phase"])
    if kind == "schmidt":
        return schmidt_phase(per, rest.split(":")[0])
    if kind == "derived":
        top = rest.split(".")[0]
        if top.startswith("etel_"):
            return "MIXED"
        if top == "N_genome_unweighted":
            return "PHASE_INVARIANT"
        return schmidt_phase(per, der[top]["schmidt_condition"])
    raise KeyError(key)


def source_conditions(key, audit, per, der):
    """organism/strain/medium/temperature/status as the SOURCE records them;
    MATCHED rows are judged on these, not on what the row itself claims."""
    kind, _, rest = key.partition(":")
    if kind == "derived":
        top = rest.split(".")[0]
        if "schmidt_condition" not in der.get(top, {}):
            return None
        kind, rest = "schmidt", der[top]["schmidt_condition"] + ":-"
    if kind == "schmidt":
        cond = rest.split(":")[0]
        ph = schmidt_phase(per, cond)
        medium = bb.ANCHOR[bb.EXP]["medium"] if cond == bb.EXP_COND else (
            bb.SCHMIDT_STAT_MEDIUM if ph == bb.STAT else cond)
        return {"organism": "Escherichia coli", "strain": "BW25113", "medium": medium,
                "temperature": "42 C" if "42" in cond else "37 C", "status": "VERIFIED",
                "org_cond": "MATCHED"}
    if kind == "audit":
        a = audit[rest]
        return {"organism": a["organism"], "strain": a["strain"], "medium": a["medium"],
                "temperature": a["temperature"], "status": a["verification_status"],
                "org_cond": a["org_cond_match"]}
    return None


def matches_anchor(src, b):
    """strain: the anchor's strain token must appear in the source strain field."""
    an = bb.ANCHOR[b]
    return (src["organism"].startswith("Escherichia coli")
            and re.search(rf"(?<!\w){re.escape(an['strain'])}(?!\w)", src["strain"]) is not None
            and src["medium"] == an["medium"] and src["temperature"] == an["temperature"])


def numbers(s):
    return [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?(?:e-?\d+)?", s)]


def close(a, b, rel):
    return abs(a - b) <= rel * max(abs(a), abs(b), 1e-300)


def source_value_ok(r, audit, per, der):
    """the row's value equals what its source file holds."""
    v, key = r["value"], r["source_key"]
    if v.startswith("vector:"):
        # the vector must be the one the source_key names
        return (v == "vector:derived.json#" + key.partition(":")[2] if key.startswith("derived:")
                else key == "audit:E05" and v == "vector:global_codon_usage_ecoli.tsv")
    vals = numbers(v)
    kind, _, rest = key.partition(":")
    if kind == "schmidt":
        cond, field = rest.split(":")
        fields = field.split("*")
        x = math.prod(per[cond][f] for f in fields)
        scale = {"1/s": 1 / 3600, "uM chains": 1e3, "uM chains/s": 1e3 / 3600}.get(r["unit"], 1.0)
        return len(vals) == 1 and close(vals[0], x * scale, 1e-5)
    if kind == "derived":
        top, _, sub = rest.partition(".")
        x = der[top][sub] if sub else der[top]["weighted_mean_codons"]
        return len(vals) == 1 and close(vals[0], x, 1e-5)
    if kind == "audit":
        a = audit[rest]
        txt = " ".join(a[c] for c in ("value_used", "range_used", "measured_in_source",
                                      "what_source_supports"))
        # abs: a range written '2.2-107' tokenizes as 2.2 and -107
        pool = [abs(x) for x in numbers(txt)]
        if r["unit"].startswith("uM chains"):
            pool += [x * 1e3 for x in pool]        # mM -> uM, as written in the audit text
        return all(any(close(abs(x), y, 5e-3) for y in pool) for x in vals)
    return False


def check_effects(r, eff):
    errs = []
    for eid in r["effect_ids"].split(";"):
        e = eff["effects"].get(eid)
        if e is None:
            errs.append(f"effect {eid} not computed")
            continue
        if e["bundle"] != r["bundle"]:
            errs.append(f"effect {eid} belongs to {e['bundle']}")
        if r["param"] not in e["params"]:
            errs.append(f"effect {eid} varies {e['params']}, not {r['param']}")
        if not e["outputs"]:
            errs.append(f"effect {eid} has no output")
        for name, o in e["outputs"].items():
            lo, hi = o["lo"], o["hi"]
            if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi):
                errs.append(f"effect {eid}/{name}: range not finite or lo > hi")
            if o["kind"] not in OUTPUT_KINDS:
                errs.append(f"effect {eid}/{name}: unknown kind {o['kind']}")
            # an effect varies a MISMATCHED input: never a derived output / prediction
            if o["kind"] == "DERIVED_OUTPUT":
                errs.append(f"effect {eid}/{name}: mismatch range labelled DERIVED_OUTPUT")
    return errs


def validate(rows, src=None):
    audit, per, der, eff = src or sources()
    errs = []
    for i, r in enumerate(rows):
        tag = f"row {i} {r.get('bundle')}/{r.get('param')}/{r.get('candidate')}"
        for c in REQUIRED:
            if not str(r.get(c, "")).strip():
                errs.append(f"{tag}: empty {c}")
        b, flag, vp = r["bundle"], r["match_flag"], r["value_phase"]
        if b not in (bb.EXP, bb.STAT):
            errs.append(f"{tag}: unknown bundle {b}")
            continue
        if flag not in FLAGS:
            errs.append(f"{tag}: match_flag {flag!r} not in {sorted(FLAGS)}")
        if vp not in PHASES:
            errs.append(f"{tag}: value_phase {vp!r} not in vocabulary")
        try:
            sph = source_phase(r["source_key"], audit, per, der)
        except KeyError:
            errs.append(f"{tag}: source_key {r['source_key']!r} does not resolve")
            continue
        if sph != vp:
            errs.append(f"{tag}: value_phase {vp} but the source records {sph}")
        # G7.2: never an opposite-phase (or non-phase) value without MISMATCHED
        if sph == OPPOSITE[b] and flag != "MISMATCHED":
            errs.append(f"{tag}: opposite-phase value ({sph}) in {b} bundle without MISMATCHED")
        if sph == OPPOSITE[b] and "PHASE_OPPOSITE" not in r["mismatch_axes"]:
            errs.append(f"{tag}: opposite-phase value lacks PHASE_OPPOSITE axis")
        if sph in ("MIXED", "IN_VITRO", "AMBIGUOUS") and flag != "MISMATCHED":
            errs.append(f"{tag}: {sph} value not flagged MISMATCHED")
        if sph == "NONE" and flag != "UNMEASURED":
            errs.append(f"{tag}: no source but flagged {flag}")
        if flag == "MATCHED":
            if sph not in (b, "PHASE_INVARIANT"):
                errs.append(f"{tag}: MATCHED but source phase is {sph}")
            if r["mismatch_axes"] != "-":
                errs.append(f"{tag}: MATCHED with mismatch axes {r['mismatch_axes']}")
            if not r["organism"].startswith("Escherichia coli"):
                errs.append(f"{tag}: MATCHED but organism {r['organism']}")
            if sph == b:
                for k, v in bb.ANCHOR[b].items():
                    if r[k] != v:
                        errs.append(f"{tag}: MATCHED but {k} {r[k]!r} != anchor {v!r}")
            src_c = source_conditions(r["source_key"], audit, per, der)
            if src_c is None:
                errs.append(f"{tag}: MATCHED but the source records no conditions")
            else:
                if sph == b and not matches_anchor(src_c, b):
                    errs.append(f"{tag}: MATCHED but the SOURCE conditions {src_c} miss the anchor")
                if src_c["status"] != "VERIFIED" or src_c["org_cond"] != "MATCHED":
                    errs.append(f"{tag}: MATCHED but the source is {src_c['status']}/"
                                f"{src_c['org_cond']} in WO-06")
            if r["wo06_status"] != "VERIFIED":
                errs.append(f"{tag}: MATCHED needs a WO-06 VERIFIED source")
        if flag == "MISMATCHED":
            if r["mismatch_axes"] in ("-", ""):
                errs.append(f"{tag}: MISMATCHED without axes")
            errs += [f"{tag}: {e}" for e in check_effects(r, eff)]
        if flag == "UNMEASURED":
            if r["value"] != "NA":
                errs.append(f"{tag}: UNMEASURED carries a value {r['value']!r} (borrowed?)")
            if r["source"] != "NONE" or r["source_key"] != "none" or vp != "NONE":
                errs.append(f"{tag}: UNMEASURED must have no source")
        elif not source_value_ok(r, audit, per, der):
            errs.append(f"{tag}: value {r['value']} not found in source {r['source_key']}")
        # stationary error input stays out of the exponential bundle
        rec = audit.get(r["source_key"].partition(":")[2], {}).get("source_record", "") \
            if r["source_key"].startswith("audit:") else ""
        if b == bb.EXP and (r["source_key"] in STIKELEATHER or STIKELEATHER[0] in r["source"]
                            or STIKELEATHER[0] in rec):
            errs.append(f"{tag}: Stikeleather stationary error input in EXPONENTIAL bundle")
        # machine-by-machine pools: never summed
        blob = (r["param"] + " " + r["machine"]).lower()
        if any(s in blob for s in ("+", "c_tot", "combined", "sum")):
            errs.append(f"{tag}: summed chaperone pool")
        if r["param"].startswith("pool_"):
            if r["machine"] not in MACHINE_NAMES:
                errs.append(f"{tag}: pool row must name exactly one machine, got {r['machine']!r}")
            elif r["param"] != next(p for p, m, _, _ in bb.MACHINES if m == r["machine"]):
                errs.append(f"{tag}: machine {r['machine']} does not own param {r['param']}")
    # coverage: each (bundle, param) has exactly one of MATCHED / UNMEASURED
    for b in (bb.EXP, bb.STAT):
        for p in bb.PARAMS:
            sel = [r for r in rows if r["bundle"] == b and r["param"] == p]
            n_m = sum(r["match_flag"] == "MATCHED" for r in sel)
            n_u = sum(r["match_flag"] == "UNMEASURED" for r in sel)
            if n_m + n_u != 1:
                errs.append(f"{b}/{p}: needs exactly one MATCHED or UNMEASURED row "
                            f"(has {n_m} MATCHED, {n_u} UNMEASURED)")
        present = {r["machine"] for r in rows if r["bundle"] == b and r["param"].startswith("pool_")}
        for m in REQUIRED_MACHINES - present:
            errs.append(f"{b}: machine {m} missing")
        extra = {r["param"] for r in rows if r["bundle"] == b} - set(bb.PARAMS)
        for p in extra:
            errs.append(f"{b}/{p}: param not in the declared list")
    return errs


def regenerated_equal(rows):
    fresh = [{k: str(v) for k, v in r.items()} for r in bb.build()]
    return fresh == rows


if __name__ == "__main__":
    rows = load()
    errs = validate(rows)
    if not regenerated_equal(rows):
        errs.append("bundles.tsv differs from build_bundles.build(): regenerate, never hand-edit")
    for e in errs:
        print("ERROR", e)
    n = {f: sum(r["match_flag"] == f for r in rows) for f in sorted(FLAGS)}
    print(f"{len(rows)} rows, {n}, {len(errs)} errors")
    sys.exit(1 if errs else 0)
