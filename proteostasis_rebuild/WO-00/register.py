"""register helper: load, save, and update CLAIM_REGISTER.tsv rows."""
import csv, sys
REG = "/storage/kiran-stuff/proteostasis_law/proteostasis_rebuild/CLAIM_REGISTER.tsv"
def load():
    with open(REG, newline="") as f: return list(csv.DictReader(f, delimiter="\t"))
def save(rows):
    with open(REG, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
def setstatus(rows, cid, status, note):
    for r in rows:
        if r["claim_id"] == cid:
            r["status"] = status; r["note"] = (r["note"] + " | " if r["note"] else "") + note; return
    raise KeyError(cid)
