"""WO-06: writes parameter_audit.tsv.

one row per (parameter, source) pair. the condition columns (organism, strain,
phase, medium, temperature) describe the SOURCE's measurement, not what the
model intends. rows with no source carry the literal NO_SOURCE there, and the
test only allows that for ASSUMED / ILLUSTRATIVE rows.

vocabularies (enforced by test_wo06.py):
  measurement_type   MEASURED | DERIVED_FROM_MEASUREMENT | ESTIMATED | MODEL_FIT |
                     THEORETICAL | ASSUMED | ILLUSTRATIVE | UNKNOWN (unmatchable citation)
  citation_check     bibliographic only: does the reference AS WRITTEN resolve to
                     one record? MATCHED | MISCITED (details wrong, or they resolve
                     to a different paper) | AMBIGUOUS | UNMATCHED (no record) |
                     NO_CITATION | LOCAL_DATA
  value_located      YES (seen in saved record/full text/data) | NO | NA
  org_cond_match     source vs model use (E. coli, balanced exponential growth,
                     in vivo, the model's quantity): MATCHED | ORGANISM_MISMATCH |
                     CONDITION_MISMATCH | QUANTITY_MISMATCH | NA; combined with '+'
  verification_status
                     VERIFIED            record MATCHED, value located, org/cond MATCHED
                     MISMATCHED_CONDITION record matched, value located, but another
                                         organism / condition / quantity
                     NOT_SUPPORTED       record matched, but it does not contain or
                                         contradicts the value attributed to it
                     MISCITED            citation misattributes: wrong bibliographic
                                         details, or a property the paper lacks
                     UNVERIFIED          no record, or value not locatable in any
                                         retrievable text
                     ASSUMED | ILLUSTRATIVE  model choice / placeholder, no source

every computed number quoted here is in audit_derived.json or schmidt_pools.json.
"""
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
D = json.loads((HERE / "audit_derived.json").read_text())
SP = json.loads((HERE / "schmidt_pools.json").read_text())["summary"]

COLS = ["row_id", "param", "model_scope", "used_at", "value_used", "range_used", "unit",
        "model_role", "cited_as", "source_record", "record_citation", "source_location",
        "organism", "strain", "phase", "medium", "temperature", "measured_in_source",
        "measurement_type", "citation_check", "citmatch_key", "value_located",
        "org_cond_match", "verification_status", "miscitation_evidence",
        "what_source_supports", "leverage", "note"]

NS = "NO_SOURCE"
TP = "proteostasis-P1/two_pool_ode.py"
AR = "proteostasis-P1/arithmetic_stress_test.py"
PM = "proteostasis-P1/paired_mc.py"
LA = "proteostasis-P1/LITERATURE_ANCHORS.md"
EPS = "proteostasis_law/envelope-paper/scripts"
RB = "proteostasis_law/proteostasis_rebuild"

g = lambda k, f=".3g": format(k, f)
sch = lambda col, key: g(SP[col][key])


def nosrc(**kw):
    """defaults for a model choice with no cited source."""
    base = dict(cited_as="none", source_record="NONE", record_citation="-",
                source_location="-", organism=NS, strain=NS, phase=NS, medium=NS,
                temperature=NS, measured_in_source="nothing (no source)",
                measurement_type="ASSUMED", citation_check="NO_CITATION", citmatch_key="-",
                value_located="NA", org_cond_match="NA", verification_status="ASSUMED",
                miscitation_evidence="-", what_source_supports="-")
    base.update(kw)
    return base


def schmidt(**kw):
    base = dict(cited_as="(WO-06 reference value; not used by legacy)",
                source_record="PMID:26641532",
                record_citation="Nat Biotechnol|2016|34",
                source_location="Supplementary Tables S6 (combined table, both datasets) and S23; records/schmidt2016/*.xlsx sha256 3280a13f...; computed by schmidt_pools.py",
                organism="Escherichia coli", strain="BW25113",
                measurement_type="DERIVED_FROM_MEASUREMENT", citation_check="MATCHED",
                citmatch_key="-", value_located="YES", miscitation_evidence="-",
                model_scope="WO-06 reference (for WO-07)", used_at="-", leverage="-")
    base.update(kw)
    return base


rows = []
add = lambda **kw: rows.append({"note": "-", **kw})   # note is the only optional column
L, K, P, T = D["lengths"], D["k_deg"], D["prot_tot"], D["S_code"]
pg, pl = D["pools"]["Glucose"], D["pools"]["LB"]

# ---------------------------------------------------------------- k_deg
kdr = (f"legacy range 1e-4..1e-3 /s = t1/2 {K['half_life_min_at_range'][0]:.1f}..{K['half_life_min_at_range'][1]:.1f} min; "
       f"baseline 3e-4 /s = t1/2 {K['half_life_min_at_baseline']:.1f} min. cited '1-10 hr' = "
       f"{K['k_equiv_to_cited_consensus_per_s'][0]:.2e}..{K['k_equiv_to_cited_consensus_per_s'][1]:.2e} /s; "
       f"baseline inside cited window: {K['baseline_inside_cited_consensus']}")
common_kdeg = dict(param="k_deg", model_scope="legacy two-pool; paired_mc",
                   used_at=f"{TP}:98 (baseline), {TP}:400 (MC), {PM}:95",
                   value_used="3e-4", range_used="logU 1e-4..1e-3", unit="1/s",
                   model_role="first-order removal (degradation) of misfolded monomer P",
                   leverage="MEDIUM")
add(row_id="L01a", **common_kdeg, cited_as="Belle 2006 PNAS (yeast)", source_record="PMID:16916930",
    record_citation="Proc Natl Acad Sci U S A|2006|103",
    source_location="abstract (PMC copy has no body text)", organism="Saccharomyces cerevisiae",
    strain="TAP-tagged collection", phase="exponential (per record: after translation inhibition)",
    medium="UNSPECIFIED_IN_RECORD", temperature="UNSPECIFIED_IN_RECORD",
    measured_in_source="half-lives of >3,750 tagged native proteins after translation inhibition",
    measurement_type="MEASURED", citation_check="MATCHED", citmatch_key="Belle2006_PNAS_103",
    value_located="NO", org_cond_match="ORGANISM_MISMATCH+QUANTITY_MISMATCH",
    verification_status="UNVERIFIED", miscitation_evidence="-",
    what_source_supports="(abstract only; the 1-10 h figure is not in it) yeast bulk native-protein half-lives; no misfolded-protein degradation rate and no E. coli data",
    note=kdr)
add(row_id="L01b", **common_kdeg, cited_as="Christiano 2014 'for bacteria'", source_record="PMID:25466257",
    record_citation="Cell Rep|2014|9", source_location="abstract; PMC4526151 Results and Fig. 1D",
    organism="Saccharomyces cerevisiae; Schizosaccharomyces pombe", strain="UNSPECIFIED_IN_RECORD",
    phase="steady-state growth", medium="UNSPECIFIED_IN_RECORD", temperature="UNSPECIFIED_IN_RECORD",
    measured_in_source="proteome-wide steady-state turnover (SILAC) in two yeasts",
    measurement_type="MEASURED", citation_check="MATCHED", citmatch_key="-", value_located="YES",
    org_cond_match="ORGANISM_MISMATCH+QUANTITY_MISMATCH", verification_status="MISCITED",
    miscitation_evidence="legacy attributes it to bacteria; the record title and abstract are the yeasts S. cerevisiae and S. pombe",
    what_source_supports=(f"median half-life 8.8 h (S. cerevisiae), i.e. k = {K['k_from_christiano_median_per_s']:.2e} /s, "
                          "14x slower than the legacy baseline; most proteins' abundance set by dilution, ~15% turned over rapidly"),
    note="full text gives S. pombe median as both 11.1 h (figure legend) and 12.0 h (text)")
add(row_id="L01c", **{**common_kdeg, "model_scope": "WO-06 candidate (not cited by legacy)"},
    cited_as="(WO-06 candidate) Goldberg 1972 PNAS", source_record="PMID:4551144",
    record_citation="Proc Natl Acad Sci U S A|1972|69", source_location="abstract (PMC copy scanned, no text)",
    organism="Escherichia coli", strain="UNSPECIFIED_IN_RECORD (incl. ram and missense-suppressor strains)",
    phase="growing and non-growing", medium="UNSPECIFIED_IN_RECORD", temperature="UNSPECIFIED_IN_RECORD",
    measured_in_source="degradation of puromycyl fragments, error-rich and analog-containing proteins vs normal proteins",
    measurement_type="MEASURED", citation_check="MATCHED", citmatch_key="-", value_located="NO",
    org_cond_match="MATCHED", verification_status="UNVERIFIED", miscitation_evidence="-",
    what_source_supports="E. coli degrades abnormal/mistranslated proteins faster than normal ones, energy-dependently, at similar rates in growing and non-growing cells; no rate constant located",
    note="right organism and right quantity class; a numeric k_d needs the full text")

# ---------------------------------------------------------------- k_obs_max
add(row_id="L02", param="k_obs_max", model_scope="legacy two-pool; paired_mc",
    used_at=f"{TP}:99, {TP}:401, {PM}:96", value_used="1e-2", range_used="logU 3e-3..8.4e-2", unit="1/s",
    model_role="maximal chaperone-assisted folding (rescue) rate of misfolded monomer",
    cited_as="Pierpaoli et al. 1997 EMBO J; 'DnaK R-state release rate for substrate turnover'",
    source_record="PMID:9843444", record_citation="Biochemistry|1998|37",
    source_location="abstract: 'kobs1 = 0.003-0.084 s-1 at pH 7.0 and 25 degreesC'",
    organism="Escherichia coli (purified DnaK)", strain="NA (in vitro)", phase="NA (in vitro)",
    medium="buffer pH 7.0", temperature="25 C",
    measured_in_source="apparent rate constant of the FIRST PHASE OF COMPLEX FORMATION between nucleotide-free R-state DnaK (1 uM) and 22-50 nM fluorescent peptides (7-22 residues)",
    measurement_type="MEASURED", citation_check="MISCITED", citmatch_key="Pierpaoli1997_EMBOJ",
    value_located="YES", org_cond_match="CONDITION_MISMATCH+QUANTITY_MISMATCH", verification_status="MISCITED",
    miscitation_evidence="no Pierpaoli EMBO J 1997 record (ecitmatch NOT_FOUND with author; control J Mol Biol 1997 fires); the range is verbatim in Pierpaoli, Gisler, Christen 1998 Biochemistry 37:16741. a 'release rate' would instead be in Gisler et al. 1998 J Mol Biol 279:833 (PMID 9642064, not fetched) or the T->R step of Pierpaoli 1997 J Mol Biol (L02b)",
    what_source_supports="an in-vitro peptide BINDING rate, not a release, turnover or folding rate; the record also quotes (from Gisler 1998, not measured there) DnaK-ATP binding kobs 0.001-7.9 s-1",
    leverage="HIGH", note="baseline 1e-2 is a choice inside the range; see L02b for a cycle rate")
add(row_id="L02b", param="k_obs_max;k_cat", model_scope="WO-06 candidate (not cited by legacy)",
    used_at="-", value_used="-", range_used="-", unit="1/s",
    model_role="rate-limiting step of the DnaK/DnaJ/GrpE cycle (candidate for k_cat)",
    cited_as="(WO-06 candidate) Pierpaoli 1997 J Mol Biol", source_record="PMID:9223639",
    record_citation="J Mol Biol|1997|269", source_location="abstract",
    organism="Escherichia coli (purified DnaK, DnaJ, GrpE)", strain="NA (in vitro)", phase="NA (in vitro)",
    medium="'conditions approximating those in the cell' (per abstract)", temperature="UNSPECIFIED_IN_RECORD",
    measured_in_source="apparent rate constants of DnaJ-triggered T->R (0.04 s-1) and GrpE-induced R->T (1.0 s-1) conversions, peptide substrates",
    measurement_type="MEASURED", citation_check="MATCHED", citmatch_key="-", value_located="YES",
    org_cond_match="CONDITION_MISMATCH", verification_status="MISMATCHED_CONDITION", miscitation_evidence="-",
    what_source_supports="the 1997 Pierpaoli paper the legacy probably meant: rate-limiting T->R 0.04 s-1 in vitro; a cycle step, not a client folding rate",
    leverage="-", note="the legacy baseline 1e-2 and the rebuild k_cat 1e-2 are 4x below this step rate")

# ---------------------------------------------------------------- C_tot
ctot = dict(param="C_tot_uM", model_scope="legacy two-pool; paired_mc; envelope scripts 09/11/12",
            used_at=f"{TP}:100, {TP}:402, {PM}:97", value_used="50", range_used="U 30..80", unit="uM",
            model_role="total effective chaperone pool (one species) in C_free and v_fold", leverage="HIGH")
add(row_id="L03a", **ctot, cited_as="Lorimer 1996: 'GroEL alone ~30 uM in exponentially growing E. coli'",
    source_record="PMID:8566548", record_citation="FASEB J|1996|10", source_location="abstract (no PMC full text)",
    organism="Escherichia coli", strain="UNSPECIFIED_IN_RECORD", phase="UNSPECIFIED_IN_RECORD",
    medium="UNSPECIFIED_IN_RECORD", temperature="UNSPECIFIED_IN_RECORD",
    measured_in_source="none: a quantitative assessment (review/calculation) from in-vitro rates and 'known quantities' of GroEL/GroES",
    measurement_type="ESTIMATED", citation_check="MATCHED", citmatch_key="Lorimer1996_FASEBJ_10",
    value_located="NO", org_cond_match="QUANTITY_MISMATCH", verification_status="UNVERIFIED", miscitation_evidence="-",
    what_source_supports="GroEL/GroES suffice for folding of 'no more than 5%' of E. coli proteins, i.e. GroEL is not a general pool",
    note=f"measured GroEL (Schmidt 2016): {sch('groL_uM_protomer','glucose')} uM protomer / {sch('groL_uM_oligomer','glucose')} uM 14-mer in glucose; among exponential conditions only 42C glucose (36.5 uM) exceeds 30 uM of protomers; LB 25.5, stationary 26.8-27.6")
add(row_id="L03b", **ctot, cited_as="'DnaK peaks at 30-50 uM' (no reference)", source_record="NONE",
    record_citation="-", source_location="-", organism=NS, strain=NS, phase=NS, medium=NS, temperature=NS,
    measured_in_source="nothing (no source)", measurement_type="ASSUMED", citation_check="NO_CITATION",
    citmatch_key="-", value_located="NA", org_cond_match="NA", verification_status="ASSUMED",
    miscitation_evidence="-",
    what_source_supports=(f"measured DnaK (Schmidt 2016, BW25113): exponential {sch('dnaK_uM_protomer','exp_min')}-"
                          f"{sch('dnaK_uM_protomer','exp_max')} uM, stationary {sch('dnaK_uM_protomer','stationary_1d')} uM; "
                          "30-50 uM is reached only in stationary phase"),
    note="an uncited number stated as fact is ASSUMED, not verified")
add(row_id="L03c", **ctot, cited_as="'combined active chaperone pool ~50 uM, up to 80 under stress' (no reference)",
    **{k: v for k, v in nosrc().items() if k not in ("cited_as",)} | {
        "what_source_supports": (f"sum of DnaK + GroEL protomers: glucose {pg['dnaK_plus_groL_protomers_uM']:.1f}, LB {pl['dnaK_plus_groL_protomers_uM']:.1f} uM; "
                                 f"as functional units (DnaK + GroEL14): glucose {pg['dnaK_plus_groEL14_uM']:.1f}, LB {pl['dnaK_plus_groEL14_uM']:.1f} uM"),
        "note": f"42C glucose {D['pools']['42°C glucose']['dnaK_plus_groL_protomers_uM']:.1f} and stationary 1 d {D['pools']['Stationary phase 1 day']['dnaK_plus_groL_protomers_uM']:.1f} uM of protomers exceed 50. adds machines with different clients and mixes protomers with oligomers; 50 uM is at the top of the protomer sum and 2-4x the functional-unit sum in exponential growth"})

# ---------------------------------------------------------------- K_d
add(row_id="L04", param="K_d_uM", model_scope="legacy two-pool; paired_mc; envelope scripts 11/12",
    used_at=f"{TP}:101, {TP}:403, {PM}:98", value_used="1", range_used="logU 0.06..2", unit="uM",
    model_role="equilibrium half-saturation of chaperone-client binding in C_free and v_fold",
    cited_as="Pierpaoli 1997 dissociation constants for DnaK-substrate peptide complexes",
    source_record="PMID:9843444", record_citation="Biochemistry|1998|37",
    source_location="abstract: 'R-state DnaK (1 microM) formed high-affinity complexes (Kd = 0.06-2 microM)'",
    organism="Escherichia coli (purified DnaK)", strain="NA (in vitro)", phase="NA (in vitro)",
    medium="buffer pH 7.0", temperature="25 C",
    measured_in_source="Kd of nucleotide-free (R-state) DnaK with 9 short fluorescent peptides; DnaK-ATP (T state) Kd 2.2-107 uM",
    measurement_type="MEASURED", citation_check="MISCITED", citmatch_key="Pierpaoli1997_EMBOJ",
    value_located="YES", org_cond_match="CONDITION_MISMATCH+QUANTITY_MISMATCH", verification_status="MISCITED",
    miscitation_evidence="attributed to Pierpaoli 1997 (EMBO J per LITERATURE_ANCHORS:17); the range is verbatim in the 1998 Biochemistry paper",
    what_source_supports="equilibrium Kd of a nucleotide-free state with peptides in vitro; the ATP state (dominant with ATP present) binds 1-2 orders weaker. WO-03: an ATP-driven cycle is set by K_M, not K_d",
    leverage="HIGH", note="the published range is used as the MC range unchanged")

# ---------------------------------------------------------------- k_agg
kagg = dict(param="k_agg_M_s", model_scope="legacy two-pool; paired_mc",
            used_at=f"{TP}:102, {TP}:404, {PM}:99", value_used="1e3", range_used="logU 3e2..3e3",
            unit="1/(M s)", model_role="second-order aggregation of misfolded monomer (drain = k_agg P^2 Prot_tot); k_nuc = k_agg by assumption",
            leverage="HIGH")
add(row_id="L05a", **kagg, cited_as="Upadhyay 2012 PLoS ONE: effective k_agg 1e2-1e3 M-1 s-1",
    source_record="PMID:22479486", record_citation="PLoS One|2012|7",
    source_location="abstract; full text PMC3315509 searched for rate constants / 'M-1' / reaction order: none",
    organism="Escherichia coli (host)", strain="UNSPECIFIED_IN_RECORD", phase="induced overexpression",
    medium="UNSPECIFIED_IN_RECORD", temperature="UNSPECIFIED_IN_RECORD",
    measured_in_source="inclusion-body size, denaturant/proteinase resistance and dye binding of overexpressed human growth hormone and asparaginase over induction time",
    measurement_type="MEASURED", citation_check="MATCHED", citmatch_key="Upadhyay2012_PLoSOne_7",
    value_located="NO", org_cond_match="CONDITION_MISMATCH+QUANTITY_MISMATCH", verification_status="NOT_SUPPORTED",
    miscitation_evidence="-", what_source_supports="qualitative IB growth kinetics of two recombinant proteins; no second-order constant",
    note="the attributed number is not in the paper")
add(row_id="L05b", **kagg, cited_as="Hoffmann & Rinas 2001: bulk misfold -> IB rates ~1e3 M-1 s-1",
    source_record="PMID:11135201", record_citation="Biotechnol Bioeng|2001|72",
    source_location="abstract (no PMC full text)", organism="Escherichia coli (host)",
    strain="UNSPECIFIED_IN_RECORD", phase="recombinant production, high-cell-density culture",
    medium="UNSPECIFIED_IN_RECORD (glucose-fed)", temperature="UNSPECIFIED_IN_RECORD",
    measured_in_source="partitioning of one recombinant human protein into soluble/insoluble fractions, fitted by a lumped kinetic model",
    measurement_type="MODEL_FIT", citation_check="AMBIGUOUS", citmatch_key="-", value_located="NO",
    org_cond_match="CONDITION_MISMATCH+QUANTITY_MISMATCH", verification_status="UNVERIFIED",
    miscitation_evidence="-",
    what_source_supports="abstract only: 'the irreversible aggregation step was found to follow first order kinetics' (exponentially decreasing production) and only transient aggregation (constant production with degradation); no second-order constant in the abstract; fitted constants may be in the full text, not retrieved",
    note="two 2001 Hoffmann/Rinas records exist: this one (Hoffmann, Posten, Rinas; IB kinetics) and 11745161 (metabolic burden). best match used")

add(row_id="L05c", **{**kagg, "value_used": "not used (1e5-1e6 excluded)", "range_used": "-",
    "model_role": "justification for NOT using amyloid elongation k+ as k_agg", "leverage": "LOW"},
    cited_as="'Cohen / Meisl / Knowles' amyloid k+ 1e5-1e6 M-1 s-1 (author names only)", source_record="NONE",
    record_citation="-", source_location="-", organism="UNKNOWN (no record)", strain="UNKNOWN (no record)",
    phase="UNKNOWN (no record)", medium="UNKNOWN (no record)", temperature="UNKNOWN (no record)",
    measured_in_source="unknown: no paper, year or journal given", measurement_type="UNKNOWN",
    citation_check="UNMATCHED", citmatch_key="-", value_located="NA", org_cond_match="NA",
    verification_status="UNVERIFIED", miscitation_evidence="-",
    what_source_supports="cannot be matched to a record without a paper identifier",
    note="sources no value used by the model; audited only because G6.2 covers every citation")

# ---------------------------------------------------------------- Prot_tot
ptn = (f"independent values: BNID 104726 = 4 mM; Milo 2013 2-4e6 proteins/um3 = {P['milo2013_mM_from_2to4e6_per_um3'][0]:.2f}-"
       f"{P['milo2013_mM_from_2to4e6_per_um3'][1]:.2f} mM; Schmidt 2016 glucose {P['schmidt_glucose_mM_wholecell']:.2f} mM. "
       f"legacy is {P['fold_low_vs_schmidt_glucose']:.1f}x (Schmidt) / {P['fold_low_vs_bnid']:.1f}x (BNID) / "
       f"{P['milo2013_mM_from_2to4e6_per_um3'][0]*1e3/300:.0f}-{P['milo2013_mM_from_2to4e6_per_um3'][1]*1e3/300:.0f}x (Milo) low")
ptot = dict(param="Prot_tot_uM", model_scope="legacy two-pool; paired_mc",
            used_at=f"{TP}:103, {TP}:405, {PM}:100", value_used="300", range_used="U 250..350", unit="uM",
            model_role="converts fraction P to concentration M (competition for chaperone) and scales the aggregation drain",
            leverage="HIGH")
add(row_id="L06a", **ptot, cited_as="Milo BioNumbers / Cell Biology by the Numbers (no BNID, no page)",
    source_record="NONE", record_citation="-", source_location="-", organism="Escherichia coli (as stated)",
    strain="UNSPECIFIED", phase="UNSPECIFIED", medium="UNSPECIFIED", temperature="UNSPECIFIED",
    measured_in_source="unknown: the citation names a database and a book without an entry",
    measurement_type="UNKNOWN", citation_check="UNMATCHED", citmatch_key="-", value_located="NA",
    org_cond_match="NA", verification_status="UNVERIFIED", miscitation_evidence="-",
    what_source_supports="the named database's own total-protein entry (BNID 104726, row L06b) gives 4 mM, 13x the legacy value. possible origins of '300' (NOT established): the BNID 104726 page points to BNID 104678 'total observed intracellular metabolite pool of 300 mM'; Zimmerman & Trach 1991 (PMID 1748995) give macromolecules 0.3-0.4 g/ml (a g/ml -> uM slip)",
    note=ptn)
add(row_id="L06b", **{**ptot, "model_scope": "WO-06 check of L06a"}, cited_as="(WO-06 check) BioNumbers BNID 104726",
    source_record="BNID:104726", record_citation="records/bionumbers_104726.html",
    source_location="BNID 104726 page, saved 2026-09-29",
    organism="Escherichia coli", strain="B/r", phase="balanced exponential growth", medium="aerobic glucose minimal medium",
    temperature="37 C", measured_in_source="calculated total protein concentration for an average cell (from Neidhardt 1996 table)",
    measurement_type="ESTIMATED", citation_check="MATCHED", citmatch_key="-", value_located="YES",
    org_cond_match="MATCHED", verification_status="VERIFIED", miscitation_evidence="-",
    what_source_supports="4 mM total protein", note="a calculation from composition data, labelled ESTIMATED, not a direct measurement")
add(row_id="L06c", **{**ptot, "model_scope": "WO-06 check of L06a"}, cited_as="(WO-06 check) Milo 2013 BioEssays",
    source_record="PMID:24114984", record_citation="Bioessays|2013|35", source_location="abstract",
    organism="bacteria, yeast, mammalian cells (general)", strain="NA (estimate)", phase="NA (estimate)",
    medium="NA (estimate)", temperature="NA (estimate)",
    measured_in_source="estimate from protein mass fraction and mean protein length: 2-4 million proteins per um3",
    measurement_type="ESTIMATED", citation_check="MATCHED", citmatch_key="-", value_located="YES",
    org_cond_match="CONDITION_MISMATCH", verification_status="MISMATCHED_CONDITION", miscitation_evidence="-",
    what_source_supports=f"generic cross-organism benchmark: {P['milo2013_mM_from_2to4e6_per_um3'][0]:.2f}-{P['milo2013_mM_from_2to4e6_per_um3'][1]:.2f} mM (unit conversion in audit_checks.py)",
    note="order-of-magnitude benchmark, not condition-resolved")
add(row_id="L06d", param="Prot_tot_uM;P_T", model_scope="WO-06 check of L06a and rebuild P_T",
    used_at="-", value_used="-", range_used="-", unit="uM", model_role="total protein concentration",
    leverage="HIGH",
    **{k: v for k, v in schmidt(
        phase="exponential (glucose is the independent value)", medium="M9 glucose", temperature="37 C",
        measured_in_source="summed MS copies/cell (total mass/cell measured for glucose only) / calculated whole-cell volume",
        org_cond_match="MATCHED", verification_status="VERIFIED",
        what_source_supports=(f"glucose {P['schmidt_glucose_mM_wholecell']:.2f} mM; 22-condition exponential "
                              f"{P['schmidt_exponential_range_mM'][0]:.2f}-{P['schmidt_exponential_range_mM'][1]:.2f} mM (whole-cell volume: lower bound on cytoplasmic)"),
        note="other conditions are not independent: Schmidt adjusted mass/cell assuming constant volumetric protein concentration; ~55% of genes quantified").items()
       if k not in ("model_scope", "used_at", "leverage")})

# ---------------------------------------------------------------- T_gen
add(row_id="L07", param="T_gen_s", model_scope="legacy two-pool; paired_mc; WO-05 mapping",
    used_at=f"{TP}:104, {TP}:406, {PM}:101", value_used="3600", range_used="U 1800..10800", unit="s",
    model_role="per-protein synthesis time in J_bare = f N (1-S) p / T_gen (WO-01: should be ln2/T_gen); legacy has no dilution",
    **nosrc(cited_as="'E. coli doubling time envelope, rich -> minimal medium' (no reference)",
            what_source_supports="measured doubling times (Schmidt 2016 Table S23, row S06): 0.3-2.7 h batch; the 30-180 min envelope is consistent, the 60 min baseline is a choice"),
    leverage="MEDIUM", note="phase-specific: meaningless for stationary phase (WO-02 self-review)")

# ---------------------------------------------------------------- N
add(row_id="L08a", param="N_prot;N", model_scope="legacy two-pool baseline; arithmetic baseline",
    used_at=f"{TP}:105, {AR}:56", value_used="300", range_used="fixed in deterministic runs", unit="codons (aa)",
    model_role="codons per protein in the error->flux mapping and in 1-(1-f)^N",
    cited_as="'median protein length (BioNumbers, Milo & Phillips)' (no BNID)", source_record="NONE",
    record_citation="-", source_location="-", organism="Escherichia coli (as stated)", strain="UNSPECIFIED",
    phase="NA (genome property)", medium="NA (genome property)", temperature="NA (genome property)",
    measured_in_source="unknown: no entry given", measurement_type="UNKNOWN", citation_check="UNMATCHED",
    citmatch_key="-", value_located="NA", org_cond_match="NA", verification_status="UNVERIFIED",
    miscitation_evidence="-",
    what_source_supports=f"the legacy's own UniProt file (row L08b): median {L['median']:.0f}, mean {L['mean']:.1f}; 300 is not the median",
    leverage="MEDIUM",
    note="D&W 2009 NRG (PMC2764353) states the average E. coli CDS is 335 codons; the burden-relevant length is synthesis-weighted, which none of these is")
add(row_id="L08b", param="N_prot;N", model_scope="legacy MC (two_pool, paired_mc, arithmetic); 06_translation_burden",
    used_at=f"{TP}:397-398, {PM}:88, {AR}:268", value_used="bootstrap", range_used=f"empirical, n={L['n']}, median {L['median']:.0f}, mean {L['mean']:.1f}",
    unit="codons (aa)", model_role="protein length distribution",
    cited_as="UniProt UP000000625 (E. coli K-12 MG1655) length table", source_record="LOCAL:proteostasis-P1/ecoli_proteome_lengths.tsv",
    record_citation="records/uniprot_UP000000625.json", source_location="Length column; UniProt proteome record geneCount 4403",
    organism="Escherichia coli", strain="K-12 MG1655", phase="NA (genome-encoded)", medium="NA (genome-encoded)",
    temperature="NA (genome-encoded)", measured_in_source="annotated protein lengths of the reference proteome",
    measurement_type="MEASURED", citation_check="LOCAL_DATA", citmatch_key="-", value_located="YES",
    org_cond_match="MATCHED", verification_status="VERIFIED", miscitation_evidence="-",
    what_source_supports="genome-encoded length distribution (unweighted by expression)", leverage="MEDIUM",
    note="file sha256 recorded in WO-00; n matches the UniProt proteome geneCount")

# ---------------------------------------------------------------- p_misfold
add(row_id="L09", param="p_baseline;p_misfold", model_scope="legacy two-pool; arithmetic; paired_mc",
    used_at=f"{TP}:106, {TP}:408 (U 0.1..0.5), {AR}:58, {AR}:273 (logU), {PM}:90", value_used="0.3",
    range_used="two-pool U 0.1..0.5; arithmetic/paired logU 0.1..0.5", unit="1 (probability)",
    model_role="probability that a nonsynonymous error yields a misfolded product",
    cited_as="Drummond & Wilke 2009 Cell", source_record="PMID:18662548;PMID:19763154",
    record_citation="Cell|2008|134;Nat Rev Genet|2009|10",
    source_location="both abstracts; full texts PMC2696314 and PMC2764353 searched for a misfolding fraction",
    organism="E. coli, yeast, worm, fly, mouse, human (comparative); review", strain="NA", phase="NA",
    medium="NA", temperature="NA",
    measured_in_source="Cell 2008: sequence-evolution/expression covariation + lattice-protein simulation; NRG 2009: review",
    measurement_type="THEORETICAL", citation_check="MISCITED", citmatch_key="DrummondWilke2009_Cell",
    value_located="YES", org_cond_match="QUANTITY_MISMATCH", verification_status="MISCITED",
    miscitation_evidence="no Drummond & Wilke Cell 2009 (ecitmatch NOT_FOUND); their Cell paper is 2008;134:341 and the 2009 paper is Nat Rev Genet 10:715",
    what_source_supports="the RANGE 0.1-0.5 matches Cell 2008 'Roughly ~10-50% of random substitutions disrupt protein function (Guo et al. 2004; Markiewicz et al. 1994)': a secondary citation about loss of FUNCTION from random MUTATIONS, not misfolding per mistranslation. NRG 2009 'up to 30% of newly synthesized proteins were rapidly degraded' is a different quantity (and a later study found 'at most a few percent'); NRG 2009 also says the failure rate 'remains essentially unknown'. Cell 2008 Fig. 6B plots the simulated foldable fraction of mistranslated chains (lattice model, figure only)",
    leverage="HIGH", note="enters J linearly; with S it is the (1-S)p factor that cancels in paired ratios (C09)")

# ---------------------------------------------------------------- S
add(row_id="L10", param="S_avg;S_syn", model_scope="legacy two-pool; arithmetic; paired_mc; WO-05 mapping (RAW only)",
    used_at=f"{TP}:107, {TP}:409 (U 0.25..0.35), {AR}:59, {AR}:272 (U 0.20..0.40), {PM}:91",
    value_used="0.30", range_used="two-pool U 0.25..0.35; arithmetic/paired U 0.20..0.40", unit="1 (fraction)",
    model_role="share of RAW decoding errors that are synonymous; applies only to raw error (WO-05)",
    **nosrc(cited_as="'wobble shielding' (no reference)",
            what_source_supports=(f"derived (audit_checks.py), standard code, uniform single-nt misreading, stops excluded: "
                                  f"all positions {T['all_positions_unweighted']:.3f} unweighted / {T['all_positions_usage_weighted']:.3f} usage-weighted; "
                                  f"third position only {T['third_position_unweighted']:.3f} / {T['third_position_usage_weighted']:.3f}")),
    leverage="HIGH",
    note="S is a property of the misreading SPECTRUM, which is not uniform (Stikeleather 2026: variant-specific third-position misreading). 0.30 matches neither limit")

# ---------------------------------------------------------------- k_clear
KC = D["k_clear"]
add(row_id="L11", param="k_clear", model_scope="legacy two-pool; paired_mc; rebuild k_dis (SCENARIO)",
    used_at=f"{TP}:108, {TP}:410, {TP}:373 (sweep), {PM}:105", value_used="4e-4", range_used="logU 1e-5..1e-3",
    unit="1/s", model_role="first-order clearance of aggregates (disaggregation + degradation lumped)",
    cited_as="Mogk et al. 1999 EMBO J: ~95% recovery of aggregated substrate in 2 h -> -ln(0.05)/7200",
    source_record="PMID:10601016", record_citation="EMBO J|1999|18",
    source_location="abstract; PMC1171757 is a scanned copy with no text; publisher PDF not retrievable",
    organism="Escherichia coli", strain="UNSPECIFIED_IN_RECORD", phase="heat stress and recovery (in vivo) and cell extracts",
    medium="UNSPECIFIED_IN_RECORD", temperature="heat shock (value not in record)",
    measured_in_source="aggregation of thermolabile proteins during heat stress; DnaK+ClpB solubilization, incl. after induced chaperone synthesis",
    measurement_type="MEASURED", citation_check="MATCHED", citmatch_key="Mogk1999_EMBOJ_18", value_located="NO",
    org_cond_match="CONDITION_MISMATCH", verification_status="UNVERIFIED", miscitation_evidence="-",
    what_source_supports="DnaK+ClpB 'quantitatively' solubilize heat-induced aggregates; 15-25% of detected protein species are thermolabile under heat stress",
    leverage="HIGH",
    note=(f"arithmetic -ln(0.05)/7200 = {KC['recomputed_minus_ln_0p05_over_7200']:.3g} /s is correct, but the 95%/2 h input is not located. "
          f"heat-shock recovery with induced ClpB; exponential ClpB6 is {sch('clpB_uM_oligomer','glucose')} uM (glucose)"))

# ---------------------------------------------------------------- A_half, k_nuc
add(row_id="L12", param="A_half", model_scope="legacy two-pool", used_at=f"{TP}:109", value_used="0.2",
    range_used="fixed", unit="1 (fraction of proteome)", model_role="saturation of the P->A drain, (1 - A/(A+A_half))",
    **nosrc(), leverage="MEDIUM", note="dropped in the conservative model (WO-02)")
add(row_id="L14", param="k_nuc", model_scope="legacy two-pool", used_at=f"{TP}:147", value_used="= k_agg",
    range_used="tied to k_agg", unit="1/(M s)", model_role="nucleation/drain constant set equal to k_agg",
    **nosrc(cited_as="'k_nuc = k_agg by assumption' (code comment)"), leverage="MEDIUM", note="-")

# ---------------------------------------------------------------- A_max
amax = dict(param="A_max", model_scope="legacy two-pool; paired_mc; envelope scripts 09/11/12",
            used_at=f"{TP}:92, {TP}:112, {TP}:411, {PM}:106", value_used="0.25", range_used="U 0.15..0.35",
            unit="1 (fraction of proteome aggregated)",
            model_role="death gate: sets the operational threshold in 99.94% of MC draws and at the headline point (WO-04)",
            leverage="HIGH")
add(row_id="L13a", **amax, cited_as="Ciryam et al. 2013 PNAS (110:E3453): ~0.2 fraction aggregated drives dysfunction",
    source_record="PMID:24183671", record_citation="Cell Rep|2013|5",
    source_location="Cell Rep abstract and full text PMC3883113 searched for a fraction/threshold: none",
    organism="Homo sapiens (proteome-wide analysis; C. elegans aggregation data)", strain="NA", phase="NA",
    medium="NA", temperature="NA",
    measured_in_source="supersaturation scores (abundance/solubility) of human proteins vs aggregation and disease pathways",
    measurement_type="THEORETICAL", citation_check="MISCITED", citmatch_key="Ciryam2013_PNAS_110_E3453",
    value_located="NO", org_cond_match="ORGANISM_MISMATCH+QUANTITY_MISMATCH", verification_status="MISCITED",
    miscitation_evidence="PNAS 2013;110:E3453 does not exist (ecitmatch NOT_FOUND); Ciryam's 2013 PNAS paper is 110:E132 (cotranslational folding, PMID 23256155); the supersaturation paper is Cell Rep 5:781",
    what_source_supports="supersaturated proteins are aggregation-prone; no aggregated-fraction viability threshold")
add(row_id="L13b", **amax, cited_as="Bednarska et al. 2013 Mol Cell (52:617): >20% aggregated proteome -> growth arrest",
    source_record="PMID:23894132;PMID:24239291", record_citation="Microbiology (Reading)|2013|159;Mol Cell|2013|52",
    source_location="both abstracts (no PMC full text of the review)", organism="bacteria (review)", strain="NA",
    phase="NA", medium="NA", temperature="NA",
    measured_in_source="Bednarska 2013 is a review; Mol Cell 52:617 is Aakre et al., a Caulobacter toxin-antitoxin paper",
    measurement_type="THEORETICAL", citation_check="MISCITED", citmatch_key="Bednarska2013_MolCell_52_617",
    value_located="NO", org_cond_match="QUANTITY_MISMATCH", verification_status="MISCITED",
    miscitation_evidence="Mol Cell 2013;52:617 resolves to PMID 24239291 (SocB toxin blocks replication), an unrelated paper; Bednarska's 2013 paper is Microbiology 159:1795",
    what_source_supports="aggregation reduces bacterial fitness (qualitative review); no 20% threshold in the abstract (review full text not checked)")
add(row_id="L13c", **amax, cited_as="Yamanaka et al. 2017 Curr Biol: chronic aggregation load >25%",
    source_record="NONE", record_citation="-", source_location="PubMed author+journal+year search: 0 records; ecitmatch NOT_FOUND",
    organism="UNKNOWN (no record)", strain="UNKNOWN (no record)", phase="UNKNOWN (no record)",
    medium="UNKNOWN (no record)", temperature="UNKNOWN (no record)", measured_in_source="unknown (no record)",
    measurement_type="UNKNOWN", citation_check="UNMATCHED", citmatch_key="Yamanaka2017_CurrBiol",
    value_located="NA", org_cond_match="NA", verification_status="UNVERIFIED", miscitation_evidence="-",
    what_source_supports="no record found", note="cannot be evaluated; not proof the paper does not exist")
add(row_id="L13d", **amax, cited_as="Stirling et al. 2018 Cell Rep (25:2242): yeast cytotoxic load ~20-30%",
    source_record="NONE", record_citation="-", source_location="PubMed journal+volume+page search: 0 records; ecitmatch NOT_FOUND; Crossref: no match",
    organism="UNKNOWN (no record)", strain="UNKNOWN (no record)", phase="UNKNOWN (no record)",
    medium="UNKNOWN (no record)", temperature="UNKNOWN (no record)", measured_in_source="unknown (no record)",
    measurement_type="UNKNOWN", citation_check="UNMATCHED", citmatch_key="Stirling2018_CellRep_25_2242",
    value_located="NA", org_cond_match="NA", verification_status="UNVERIFIED", miscitation_evidence="-",
    what_source_supports="no record found; as stated it would be yeast in any case")
add(row_id="L13e", **amax, **nosrc(cited_as="(summary) A_max as used",
    what_source_supports="none of the four citations supports a numeric threshold: two MISCITED, two UNMATCHED"),
    note="treated as ASSUMED in all later WOs")

# ---------------------------------------------------------------- arithmetic-only
add(row_id="A01", param="P_correct", model_scope="legacy arithmetic; paired_mc", used_at=f"{AR}:57, {AR}:270, {PM}:92",
    value_used="0.70", range_used="U 0.50..0.90", unit="1 (fraction)",
    model_role="minimum tolerable fraction of error-free-folded protein", **nosrc(), leverage="HIGH",
    note="sets the arithmetic threshold directly: -ln(P)/N")

# ---------------------------------------------------------------- error inputs
add(row_id="E01", param="mu_usage_weighted (f_codon input)", model_scope="envelope 06/09/11/12; WO-05 G5.5",
    used_at=f"{EPS}/06_translation_burden.py:37; {RB}/WO-05/run_wo05.py", value_used="6.334e-4",
    range_used="per-codon 3.3e-5..2.0e-2", unit="1/codon (MS-detected substitutions)",
    model_role="observed error rate at which headroom is evaluated",
    cited_as="Landerer et al. 2024 MBE Data_S2 'error detection rate' (usage-weighted)", source_record="PMID:38421032",
    record_citation="Mol Biol Evol|2024|41",
    source_location="Data_S2 xlsx sha256 08ef9455...; reproduced from Data_S4 in WO-05",
    organism="Escherichia coli", strain="MIXED (80 PRIDE datasets)", phase="MIXED (80 PRIDE datasets)",
    medium="MIXED (80 PRIDE datasets)", temperature="MIXED (80 PRIDE datasets)",
    measured_in_source="per-codon mean over datasets WITH >=1 detected substitution of substitution-PSMs/covering-PSMs",
    measurement_type="DERIVED_FROM_MEASUREMENT", citation_check="MATCHED", citmatch_key="-", value_located="YES",
    org_cond_match="CONDITION_MISMATCH+QUANTITY_MISMATCH", verification_status="MISMATCHED_CONDITION",
    miscitation_evidence="-",
    what_source_supports="an MS detection rate at substitution level, conditional on detection: upward-selected, lower bound on f_sub, pooled across conditions (WO-05)",
    leverage="HIGH", note="Landerer's model: 'on average 20% to 23% of proteins' carry >=1 misincorporation (abstract; stated for E. coli and S. cerevisiae together), vs legacy 16-18%")
add(row_id="E02", param="f_eTEL alternatives", model_scope="WO-05 G5.5",
    used_at=f"{RB}/WO-05/run_wo05.py:184-207", value_used="2.526e-4 (zero-inclusive); 1.414e-4 (PSM-pooled)",
    range_used="per-dataset 6.7e-6..3.29e-3", unit="1/codon (MS-detected substitutions)",
    model_role="selection-bias-free inputs for the recomputed headroom",
    cited_as="Landerer et al. 2024 MBE Data_S4 per-dataset counts", source_record="PMID:38421032",
    record_citation="Mol Biol Evol|2024|41", source_location="Data_S4 codon_counts.csv / substitution_errors.csv (WO-05)",
    organism="Escherichia coli", strain="MIXED (80 PRIDE datasets)", phase="MIXED (80 PRIDE datasets)",
    medium="MIXED (80 PRIDE datasets)", temperature="MIXED (80 PRIDE datasets)",
    measured_in_source="per-dataset base_count and error_count per codon",
    measurement_type="DERIVED_FROM_MEASUREMENT", citation_check="MATCHED", citmatch_key="-", value_located="YES",
    org_cond_match="CONDITION_MISMATCH", verification_status="MISMATCHED_CONDITION", miscitation_evidence="-",
    what_source_supports="aggregate eTEL statistics; no per-dataset growth-phase annotation used",
    leverage="HIGH", note="WO-07 needs per-dataset condition tags from PRIDE before any phase-matched use")
add(row_id="E03", param="f_stationary (Stikeleather)", model_scope="WO-05 (kept separate, no headroom)",
    used_at=f"{RB}/WO-05/run_wo05.py:144", value_used="1.82e-3 (SE 5.92e-5)", range_used="-",
    unit="1/codon (MS-detected substitutions)", model_role="condition-specific error rate; not combined with exponential parameters",
    cited_as="Stikeleather, Ali, Ho, Licknack, Lynch 2026 NAR 54(13) gkag674", source_record="PMID:42406629",
    record_citation="Nucleic Acids Res|2026|54",
    source_location="PMC13335486 Results: 'mean translation-error rate for the wild type (1.82 x 10-3 per codon, SE = 5.92 x 10-5)'; Methods: culture conditions",
    organism="Escherichia coli", strain="Xac [ara, delta(lac-proAB), gyrA, rpoB, argE(am)] (gift of H. Zaher)",
    phase="stationary (overnight)", medium="LB (Miller)", temperature="37 C",
    measured_in_source="total detected substitutions / total sites sampled, MS, I/L merged, 3 biological replicates",
    measurement_type="MEASURED", citation_check="MATCHED", citmatch_key="-", value_located="YES",
    org_cond_match="MATCHED", verification_status="VERIFIED", miscitation_evidence="-",
    what_source_supports="a stationary-phase WT rate; the SE exponent (1e-5) is now read in the saved PMC XML",
    leverage="HIGH", note="VERIFIED for its stated role (stationary). combined with exponential parameters it is MISMATCHED (WO-07)")
add(row_id="E04", param="f_window [1e-4, 1e-3]; f_codon=1e-4; f_obs=1e-4", model_scope="legacy two-pool crosscheck; paired_mc; envelope 09/11",
    used_at=f"{TP}:513, {PM}:145, {EPS}/11_headroom_sensitivity.py:47-54, {EPS}/09_supraadditivity.py:83",
    value_used="1e-4 and 1e-3", range_used="window", unit="1/codon (kind unstated: raw vs substitution, WO-05)",
    model_role="'quoted observed window' of E. coli error rates", cited_as="none in the scripts",
    source_record="NONE", record_citation="-", source_location="-", organism="UNSPECIFIED", strain="UNSPECIFIED",
    phase="UNSPECIFIED", medium="UNSPECIFIED", temperature="UNSPECIFIED",
    measured_in_source="unknown (no source)", measurement_type="UNKNOWN", citation_check="NO_CITATION",
    citmatch_key="-", value_located="NA", org_cond_match="NA", verification_status="UNVERIFIED",
    miscitation_evidence="-",
    what_source_supports="D&W 2008 Cell (PMC2696314) states errors 'occur at rates of one per 10^3-10^4 codons', citing reporter studies; single-codon, mixed organisms",
    leverage="MEDIUM", note="a secondary summary exists; the legacy never cites one")
add(row_id="E05", param="codon usage weights", model_scope="envelope 06; WO-05 usage weighting",
    used_at=f"{EPS}/06_translation_burden.py:31-35", value_used="61 sense-codon genomic counts",
    range_used="-", unit="codon counts", model_role="weights mu to a proteome mean",
    cited_as="global_codon_usage_ecoli.tsv (generator not in repo)",
    source_record="LOCAL:proteostasis_law/envelope-paper/data/raw/global_codon_usage_ecoli.tsv",
    record_citation="LOCAL:ecoli_k12_cds.fna (RefSeq NC_000913.3)",
    source_location="recount from ecoli_k12_cds.fna in audit_checks.py", organism="Escherichia coli",
    strain="K-12 MG1655 (NC_000913.3)", phase="NA (genome-encoded)", medium="NA (genome-encoded)",
    temperature="NA (genome-encoded)", measured_in_source="codon counts over annotated CDS",
    measurement_type="MEASURED", citation_check="LOCAL_DATA", citmatch_key="-", value_located="YES",
    org_cond_match="QUANTITY_MISMATCH", verification_status="MISMATCHED_CONDITION", miscitation_evidence="-",
    what_source_supports=(f"recount reproduces the table to max |diff| {D['codon_usage']['max_abs_count_diff']} counts "
                          f"in {D['codon_usage']['n_codons_differing']}/61 codons (tsv total {D['codon_usage']['tsv_total']}, recount {D['codon_usage']['fna_total_sense_excl_stop']})"),
    leverage="MEDIUM", note="genomic usage, not translation-weighted: a per-codon proteome error rate needs expression weighting")

# ---------------------------------------------------------------- envelope-script grids
add(row_id="E06", param="ANCHORINGS (C_tot, K_d)", model_scope="envelope 11_headroom_sensitivity",
    used_at=f"{EPS}/11_headroom_sensitivity.py:60-73", value_used="(50,1) (50,10) (5,1) (2,1) (1,1) (50,50)",
    range_used="grid", unit="uM", model_role="sensitivity grid for the chaperone arm", **nosrc(), leverage="MEDIUM",
    note="labelled as alternatives, not measurements")
add(row_id="E07", param="theta", model_scope="envelope 12_chaperone_availability", used_at=f"{EPS}/12_chaperone_availability.py:60",
    value_used="0, 0.5, 0.8, 0.9, 0.95, 0.98, 0.99", range_used="swept", unit="1 (fraction of pool committed)",
    model_role="fraction of chaperone committed elsewhere", **nosrc(), leverage="HIGH",
    note="WO-03: theta is an output of nascent-chain competition; not measured")
add(row_id="E08", param="error_factor; capacity_factor", model_scope="envelope 09_supraadditivity",
    used_at=f"{EPS}/09_supraadditivity.py:240", value_used="3, 3 (sweep 1)", range_used="1.5..20 grid",
    unit="1 (fold)", model_role="perturbation design", **nosrc(), leverage="LOW", note="design choice")

# ---------------------------------------------------------------- rebuild WO-02 SCENARIO
MD = f"{RB}/WO-02/model.py"
ill = lambda **kw: nosrc(measurement_type="ILLUSTRATIVE", verification_status="ILLUSTRATIVE", **kw)
add(row_id="R01", param="mu", model_scope="rebuild WO-02..05 SCENARIO", used_at=f"{MD}:118",
    value_used="ln2/3600 = 1.925e-4", range_used="-", unit="1/s", model_role="growth = dilution rate",
    **ill(cited_as="legacy T_gen 60 min", what_source_supports="inherits L07 (ASSUMED); measured batch doubling 0.3-2.7 h (row S06)"),
    leverage="HIGH", note="balanced growth only")
add(row_id="R02", param="P_T", model_scope="rebuild WO-02..05 SCENARIO", used_at=f"{MD}:119", value_used="3000",
    range_used="-", unit="uM", model_role="total protein",
    **ill(cited_as="'order of magnitude only' placeholder",
          what_source_supports=f"consistent with Schmidt glucose {P['schmidt_glucose_mM_wholecell']:.2f} mM and BNID 104726 4 mM (rows L06b-d)"),
    leverage="HIGH", note="placeholder happens to agree with the checks; stays ILLUSTRATIVE (not calibrated)")
add(row_id="R03", param="C_T", model_scope="rebuild WO-02..05 SCENARIO", used_at=f"{MD}:120", value_used="50",
    range_used="-", unit="uM", model_role="total effective chaperone", **ill(cited_as="legacy C_tot",
    what_source_supports="inherits L03a-c"), leverage="HIGH", note="-")
add(row_id="R04", param="k_on;k_off", model_scope="rebuild WO-02..05 SCENARIO", used_at=f"{MD}:121-122",
    value_used="1; 1", range_used="-", unit="1/(uM s); 1/s", model_role="chaperone-client binding",
    **ill(cited_as="chosen so k_off/k_on = legacy K_d = 1 uM", what_source_supports="inherits L04"),
    leverage="MEDIUM", note="-")
add(row_id="R05", param="k_cat", model_scope="rebuild WO-02..05 SCENARIO", used_at=f"{MD}:123", value_used="1e-2",
    range_used="-", unit="1/s", model_role="completion of one chaperone cycle",
    **ill(cited_as="legacy k_obs_max", what_source_supports="inherits L02 (a binding rate); rate-limiting DnaK cycle step in vitro 0.04 s-1 (L02b)"),
    leverage="HIGH", note="-")
add(row_id="R06", param="phi", model_scope="rebuild WO-02..05 SCENARIO", used_at=f"{MD}:124", value_used="1.0",
    range_used="-", unit="1 (probability)", model_role="productive fraction of completed cycles",
    **ill(cited_as="legacy has no partition"), leverage="MEDIUM", note="-")
add(row_id="R07", param="k_d", model_scope="rebuild WO-02..05 SCENARIO", used_at=f"{MD}:125", value_used="3e-4",
    range_used="-", unit="1/s", model_role="degradation of U", **ill(cited_as="legacy k_deg",
    what_source_supports="inherits L01a-c"), leverage="MEDIUM", note="-")
add(row_id="R08", param="k_a", model_scope="rebuild WO-02..05 SCENARIO", used_at=f"{MD}:126", value_used="1e-3",
    range_used="-", unit="1/(uM s)", model_role="aggregation, flux k_a U^2", **ill(cited_as="legacy k_agg 1e3 /(M s)",
    what_source_supports="inherits L05a-b (no verified second-order constant)"), leverage="HIGH", note="-")
add(row_id="R09", param="k_dis", model_scope="rebuild WO-02..05 SCENARIO", used_at=f"{MD}:127", value_used="4e-4",
    range_used="-", unit="1/s", model_role="disaggregation A -> U", **ill(cited_as="legacy k_clear",
    what_source_supports="inherits L11"), leverage="HIGH", note="-")
add(row_id="R10", param="k_dA;k_mis", model_scope="rebuild WO-02..05 SCENARIO", used_at=f"{MD}:128-129",
    value_used="0; 0", range_used="-", unit="1/s", model_role="aggregate degradation; spontaneous unfolding",
    **ill(cited_as="switched off"), leverage="MEDIUM", note="-")
add(row_id="R11", param="eps", model_scope="rebuild WO-02..04 SCENARIO", used_at=f"{MD}:130", value_used="0.04",
    range_used="-", unit="1 (probability)", model_role="share of new chains entering U",
    **ill(cited_as="placeholder; WO-05 derives it"), leverage="HIGH", note="-")
add(row_id="R12", param="s_P;s_C", model_scope="rebuild WO-02..05", used_at=f"{MD}:136-137",
    value_used="mu*P_T; mu*C_T", range_used="-", unit="uM/s", model_role="synthesis fluxes at balanced size",
    **ill(cited_as="derived from mu, P_T, C_T"), leverage="HIGH", note="assumes negligible turnover of native protein")

# ---------------------------------------------------------------- rebuild WO-03 / WO-04
add(row_id="R13", param="k_onT;k_offT;k_onD;k_offD;k_h;k_h0;k_ex", model_scope="rebuild WO-03 four-state cycle",
    used_at=f"{RB}/WO-03/run_wo03.py:52-53", value_used="10;100;0.01;0.01;10;0.01;0.1",
    range_used="-", unit="1/(uM s) or 1/s", model_role="DnaK-like nucleotide cycle (ultra-affinity demo)",
    **ill(cited_as="'illustrative rates only'",
          what_source_supports="K_dT 10 uM lies in the DnaK-ATP 2.2-107 uM range and K_dD 1 uM in the R-state 0.06-2 uM range (Pierpaoli 1998, L04); k_h = 10 s-1 is 250x the measured rate-limiting T->R 0.04 s-1 (L02b)"),
    leverage="LOW", note="only the existence of ultra-affinity was claimed, not its size")
add(row_id="R14", param="nu_c;k_onX;k_offX;k_catX;k_fX;k_xu", model_scope="rebuild WO-03 nascent competition",
    used_at=f"{RB}/WO-03/run_wo03.py:77-78", value_used="0..0.4;1;1;0.05;0.01;1e-3", range_used="-",
    unit="uM/s; 1/(uM s); 1/s", model_role="nascent-chain engagement of the same pool",
    **ill(cited_as="'illustrative nascent rates'"), leverage="LOW", note="theta values reported as not about E. coli")
add(row_id="R15", param="eps;mu;P_T;C_T;k_on;k_off;k_cat;phi;k_d;k_a;k_dis;k_dA;k_mis;k_onA;k_offA;k_dcat (WO-04 domain)",
    model_scope="rebuild WO-04 scans", used_at=f"{RB}/WO-04/bifurcation.py:170-180",
    value_used="log-uniform ranges", range_used="declared in bifurcation.py:sample", unit="as WO-01",
    model_role="parameter domain for steady-state counting", **nosrc(cited_as="declared domain"),
    leverage="LOW", note="a mathematical domain, not a biological claim")
add(row_id="R16", param="ultra-affinity (qualitative)", model_scope="rebuild WO-03 REPORT",
    used_at=f"{RB}/WO-03/REPORT.md (G3.2)", value_used="NA (qualitative)", range_used="-", unit="NA",
    model_role="attribution of the driven-cycle affinity argument",
    cited_as="De Los Rios & Barducci 2014 eLife (WO-03: 'UNVERIFIED until WO-06')", source_record="PMID:24867638",
    record_citation="Elife|2014|3", source_location="abstract",
    organism="Hsp70 (theory)", strain="NA", phase="NA", medium="NA", temperature="NA",
    measured_in_source="theoretical kinetic model of Hsp70 with ATP hydrolysis",
    measurement_type="THEORETICAL", citation_check="MATCHED", citmatch_key="-", value_located="YES",
    org_cond_match="MATCHED", verification_status="VERIFIED", miscitation_evidence="-",
    what_source_supports="'energy consumption can indeed decrease the dissociation constant ... by several orders of magnitude'",
    leverage="LOW", note="VERIFIED as an attribution of a theoretical claim, not as a measured parameter")

# ---------------------------------------------------------------- Schmidt reference values
for rid, par, col, role in [
        ("S01", "DnaK pool", "dnaK_uM_protomer", "DnaK concentration"),
        ("S02", "GroEL pool", "groL_uM_protomer", "GroEL concentration (protomer; 14-mer = /14)"),
        ("S03", "ClpB pool", "clpB_uM_protomer", "ClpB concentration (protomer; hexamer = /6)")]:
    s = SP[col]
    add(row_id=rid, param=par, model_scope="WO-06 reference (for WO-07)", used_at="-",
        value_used=f"glucose {g(s['glucose'])}; LB {g(s['LB'])}", range_used=(
            f"exponential {g(s['exp_min'])}..{g(s['exp_max'])}; 42C {g(s['42C_glucose'])}; "
            f"stationary 1d {g(s['stationary_1d'])}, 3d {g(s['stationary_3d'])}"),
        unit="uM (whole-cell volume)", model_role=role,
        **{k: v for k, v in schmidt(
            phase="exponential (20 conditions) and stationary (1 d, 3 d)", medium="M9 + carbon sources; LB; chemostat",
            temperature="37 C (one 42 C condition)",
            measured_in_source="MS copies/cell (SRM-calibrated) / calculated volume",
            org_cond_match="MATCHED", verification_status="VERIFIED",
            what_source_supports="condition-resolved pools; see schmidt_pools.json",
            note="phase must be matched in WO-07; volumes calculated, periplasm included").items()
           if k not in ("model_scope", "used_at", "leverage")}, leverage="HIGH")
add(row_id="S04", param="co-chaperones DnaJ, GrpE", model_scope="WO-06 reference (for WO-07)", used_at="-",
    value_used=f"glucose DnaJ {sch('dnaJ_uM_protomer','glucose')}, GrpE {sch('grpE_uM_protomer','glucose')}",
    range_used=f"DnaK:DnaJ glucose {pg['dnaK_over_dnaJ']:.0f}, LB {pl['dnaK_over_dnaJ']:.0f}", unit="uM (whole-cell volume)",
    model_role="co-chaperone limitation of the DnaK cycle",
    **{k: v for k, v in schmidt(phase="exponential and stationary", medium="M9 + carbon sources; LB", temperature="37 C",
                                 measured_in_source="MS copies/cell / calculated volume", org_cond_match="MATCHED",
                                 verification_status="VERIFIED",
                                 what_source_supports="DnaJ is ~25-30x substoichiometric to DnaK in exponential growth (Pierpaoli 1998 JBC assumed 10:1:3)",
                                 note="relevant to whether k_cat is set by DnaK or by DnaJ").items()
       if k not in ("model_scope", "used_at", "leverage")}, leverage="MEDIUM")
add(row_id="S06", param="doubling time / growth rate", model_scope="WO-06 reference (checks L07, R01)", used_at="-",
    value_used="glucose 0.58 /h; LB 1.9 /h", range_used="batch 0.26..1.9 /h (doubling 0.36..2.7 h); chemostat 0.12..0.5 /h",
    unit="1/h", model_role="dilution rate mu",
    **{k: v for k, v in schmidt(phase="exponential (batch, >=10 doublings) and chemostat", medium="M9 + carbon sources; LB",
                                 temperature="37 C", measurement_type="MEASURED",
                                 source_location="Supplementary Table S23 (records/schmidt2016/*.xlsx)",
                                 measured_in_source="growth rate from OD, triplicates",
                                 org_cond_match="MATCHED", verification_status="VERIFIED",
                                 what_source_supports="the legacy 30-180 min envelope lies within measured doubling times",
                                 note="-").items() if k not in ("model_scope", "used_at", "leverage")}, leverage="HIGH")


def main():
    for r in rows:
        missing = [c for c in COLS if c not in r]
        extra = [c for c in r if c not in COLS]
        assert not missing and not extra, (r.get("row_id"), missing, extra)
    with open(HERE / "parameter_audit.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: str(r[c]).replace("\t", " ").replace("\n", " ") for c in COLS})
    print(f"wrote {len(rows)} rows")


if __name__ == "__main__":
    main()
