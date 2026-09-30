"""WO-06: measured chaperone pools and total protein concentration, from a
locally saved primary supplement (Schmidt et al. 2016 Nat Biotechnol 34:104,
PMID 26641532, PMC4888949; Supplementary Tables xlsx, sha256 below).

what is measured vs inferred in this source (read from its Methods, saved as
records/pmc_PMC4888949.xml):
  * copies/cell per protein: MS, 41 proteins by SRM + isotope dilution, the
    rest by summed MS intensity calibrated on those 41 (Table S6, dataset 2,
    strain BW25113, biological triplicates).
  * total protein mass per cell: measured by LC-MS in triplicate for GLUCOSE
    only; for every other condition it was ADJUSTED "assuming that the
    volumetric protein concentration is condition independent". so total
    protein concentration is an independent measurement for glucose only.
  * single-cell volume: CALCULATED from growth rate (Volkmer & Heinemann 2011),
    not measured per sample (Table S23 footnote 1). whole-cell volume includes
    the periplasm, so every concentration below is per whole-cell volume and is
    a lower bound on the cytoplasmic concentration.

outputs schmidt_pools.json. every value is an MS copy number divided by a
calculated volume; nothing here is fitted to the legacy model.
"""
import hashlib
import json
import warnings
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
XLSX = HERE / "records" / "schmidt2016" / "NIHMS65833-supplement-Supplementary_tables.xlsx"
XLSX_SHA256 = "3280a13ff67a73f25440cff6ee73fb99b5ce3ef57854213dbbf6272be241912f"
N_A = 6.02214076e23

# uniprot accession -> (gene, subunits per functional oligomer)
CHAPERONES = {
    "P0A6Y8": ("dnaK", 1),
    "P0A6F5": ("groL", 14),   # GroEL tetradecamer
    "P0A6F9": ("groS", 7),    # GroES heptamer
    "P63284": ("clpB", 6),    # ClpB hexamer
    "P08622": ("dnaJ", 2),
    "P09372": ("grpE", 2),
    "P0A850": ("tig", 1),     # trigger factor
    "P0A6Z3": ("htpG", 2),
}
# conditions reported per strain BW25113 in Table S23 whose names match Table S6
PHASE = {"Stationary phase 1 day": "stationary", "Stationary phase 3 days": "stationary"}


def load():
    assert hashlib.sha256(XLSX.read_bytes()).hexdigest() == XLSX_SHA256, "supplement changed"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s6 = pd.read_excel(XLSX, sheet_name="Table S6", header=None)
        s23 = pd.read_excel(XLSX, sheet_name="Table S23", header=None)
    # table S6: row 2 is the header; columns 7..28 are copies/cell, 29..50 fg/cell
    hdr = [str(h).strip() for h in s6.iloc[2]]
    copies = s6.iloc[3:, 7:29].apply(pd.to_numeric, errors="coerce")
    copies.columns = hdr[7:29]
    mass = s6.iloc[3:, 29:51].apply(pd.to_numeric, errors="coerce")
    mass.columns = hdr[29:51]
    acc = s6.iloc[3:, 0].astype(str).str.strip()
    copies.index = mass.index = acc
    # ClpB has a second row (isoform ClpB-3, 1 peptide, <1% of ClpB); sum rows per accession
    copies = copies.groupby(level=0).sum(min_count=1)
    mass = mass.groupby(level=0).sum(min_count=1)
    # table S23: BW25113 rows only; column 4 = volume (fl), 2 = growth rate (1/h)
    t = s23.iloc[3:29, [0, 1, 2, 4]]
    t.columns = ["cond", "strain", "growth_per_h", "volume_fl"]
    t = t[t.strain.astype(str).str.strip() == "BW25113"].copy()
    t["cond"] = t.cond.astype(str).str.strip().str.rstrip("3").str.strip()
    t["cond"] = t.cond.str.replace("chemostat", "Chemostat")
    vol = dict(zip(t.cond, t.volume_fl.astype(float)))
    gr = dict(zip(t.cond, pd.to_numeric(t.growth_per_h, errors="coerce")))
    return copies, mass, vol, gr


def uM(copies_per_cell, volume_fl):
    return copies_per_cell / (N_A * volume_fl * 1e-15) * 1e6


def run():
    copies, mass, vol, gr = load()
    rows = []
    for cond in copies.columns:
        if cond not in vol:
            raise KeyError(f"no volume for condition {cond!r}")
        v = vol[cond]
        tot = float(copies[cond].sum(skipna=True))
        r = {"condition": cond, "phase": PHASE.get(cond, "exponential/steady"),
             "growth_per_h": float(gr[cond]), "volume_fl_calculated": v,
             "total_copies_per_cell": tot,
             "total_protein_mM_wholecell": uM(tot, v) / 1e3,
             "total_mass_fg_per_cell": float(mass[cond].sum(skipna=True)),
             "total_is_independent_measurement": cond == "Glucose"}
        for a, (g, n) in CHAPERONES.items():
            c = float(copies.loc[a, cond]) if a in copies.index else float("nan")
            r[f"{g}_copies"] = c
            r[f"{g}_uM_protomer"] = uM(c, v)
            r[f"{g}_uM_oligomer"] = uM(c, v) / n
        rows.append(r)
    df = pd.DataFrame(rows)
    exp = df[df.phase != "stationary"]
    summ = {}
    for col in ["total_protein_mM_wholecell"] + [f"{g}_uM_protomer" for g, _ in CHAPERONES.values()] \
            + ["groL_uM_oligomer", "clpB_uM_oligomer"]:
        summ[col] = {"exp_min": float(exp[col].min()), "exp_median": float(exp[col].median()),
                     "exp_max": float(exp[col].max()),
                     "glucose": float(df.loc[df.condition == "Glucose", col].iloc[0]),
                     "LB": float(df.loc[df.condition == "LB", col].iloc[0]),
                     "42C_glucose": float(df.loc[df.condition == "42°C glucose", col].iloc[0]),
                     "stationary_1d": float(df.loc[df.condition == "Stationary phase 1 day", col].iloc[0]),
                     "stationary_3d": float(df.loc[df.condition == "Stationary phase 3 days", col].iloc[0])}
    out = {"source": "Schmidt et al. 2016 Nat Biotechnol 34:104-110, PMID 26641532, "
                     "Supplementary Tables S6 (dataset 2, BW25113) and S23",
           "xlsx_sha256": XLSX_SHA256, "n_proteins": int(copies.shape[0]),
           "caveats": ["volume calculated (Volkmer & Heinemann 2011), whole-cell incl. periplasm",
                       "total protein mass/cell measured for glucose only; other conditions "
                       "adjusted assuming condition-independent volumetric protein concentration",
                       "only ~55% of genes quantified; unquantified proteins make totals a lower bound"],
           "summary": summ, "per_condition": rows}
    (HERE / "schmidt_pools.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    o = run()
    for k, v in o["summary"].items():
        print(f"{k:28s} " + "  ".join(f"{a}={b:.4g}" for a, b in v.items()))
