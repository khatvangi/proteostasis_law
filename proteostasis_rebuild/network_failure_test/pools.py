"""measured EXPONENTIAL (glucose, BW25113) machine pools, used ONLY as relative
pool sizes.

sources
  * chaperones (DnaK, DnaJ2, GrpE2, TF, GroEL14, GroES7, HtpG2, ClpB6): the
    MATCHED EXPONENTIAL rows of WO-07/bundles.tsv (read, never written).
  * degradation and sequestration machines (Lon, ClpX, ClpP, HslU, HslV,
    IbpA, IbpB): not in the WO-07 bundle. they are read here from the SAME
    Schmidt 2016 supplement (Table S6 copies/cell, Table S23 calculated volume,
    glucose column), through WO-06's sha256-checked loader. same source, same
    condition, same volume convention as the WO-07 rows; a test re-derives the
    WO-07 chaperone values through this path to prove the two agree.

what these numbers are NOT: capacities. a capacity is pool x per-unit
throughput for a given client, and no matched in-vivo throughput exists
(WO-07: every kinetic parameter UNMEASURED). the model therefore multiplies each
pool by an explicit, swept multiplier s_m.

oligomer convention: functional units (DnaK monomer, GroEL14, ClpB6, Lon6,
ClpX6, ClpP14, HslU6, HslV12). IbpA/IbpB form polydisperse oligomers with no
fixed stoichiometry, so they are kept in protomers. this convention is itself
arbitrary for capacity purposes, which is one more reason s_m is swept.
"""
import csv
import sys
from pathlib import Path

sys.dont_write_bytecode = True  # never write __pycache__ into WO-06

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BUNDLES = ROOT / "WO-07" / "bundles.tsv"

# machine -> WO-07 bundle param
WO07_PARAMS = {
    "DnaK": "pool_DnaK", "DnaJ": "pool_DnaJ2", "GrpE": "pool_GrpE2", "TF": "pool_TF",
    "GroEL": "pool_GroEL14", "GroES": "pool_GroES7", "HtpG": "pool_HtpG2", "ClpB": "pool_ClpB6",
}
# machine -> (uniprot accession, subunits per functional unit), read from Schmidt S6
SCHMIDT_EXTRA = {
    "Lon": ("P0A9M0", 6), "ClpX": ("P0A6H1", 6), "ClpP": ("P0A6G7", 14),
    "HslU": ("P0A6H5", 6), "HslV": ("P0A7B8", 12),
    "IbpA": ("P0C054", 1), "IbpB": ("P0C058", 1),
}
# accessions of the WO-07 machines, used only for the cross-check test
WO07_ACC = {"DnaK": ("P0A6Y8", 1), "DnaJ": ("P08622", 2), "GrpE": ("P09372", 2),
            "TF": ("P0A850", 1), "GroEL": ("P0A6F5", 14), "GroES": ("P0A6F9", 7),
            "HtpG": ("P0A6Z3", 2), "ClpB": ("P63284", 6)}


def readWO07Pools():
    """EXPONENTIAL MATCHED machine pools from WO-07 (uM functional units)."""
    out = {}
    with open(BUNDLES) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["bundle"] == "EXPONENTIAL" and r["match_flag"] == "MATCHED" \
                    and r["param"] in WO07_PARAMS.values():
                out[r["param"]] = float(r["value"])
    return {m: out[p] for m, p in WO07_PARAMS.items()}


def readSchmidtGlucose(table):
    """uM functional units, glucose condition, via WO-06's checked loader."""
    sys.path.insert(0, str(ROOT / "WO-06"))
    import schmidt_pools as sp  # read-only; sha256 asserted inside load()
    copies, _mass, vol, _gr = sp.load()
    v = vol["Glucose"]
    return {m: sp.uM(float(copies.loc[acc, "Glucose"]), v) / n for m, (acc, n) in table.items()}


def findPools():
    """all 15 machine pools (uM functional units), glucose exponential."""
    p = readWO07Pools()
    p.update(readSchmidtGlucose(SCHMIDT_EXTRA))
    return p


if __name__ == "__main__":
    for m, v in findPools().items():
        print(f"{m:6s} {v:10.4g} uM")
