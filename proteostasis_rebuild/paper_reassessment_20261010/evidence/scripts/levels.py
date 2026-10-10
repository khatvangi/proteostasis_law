"""Assign one evidentiary level per evidence row, plus the model's assumed inputs.

L1 DIRECT DATA   recomputed here from a local original data file we can hash
L2 LITERATURE    a published value read from the article, not recomputed here
L3 ASSUMED       a modelling convention or illustrative choice; no measurement
Levels are derived mechanically from provenance_status/analysis_status, never
hand-set per row. Assumed rows are enumerated explicitly because they appear in
no source table.
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]

# analysis_status -> level; recomputation from a local original file is the only L1 route
BY_ANALYSIS = {'recomputed_from_local_original_workbook': 'L1_direct_data',
               'recomputed_from_published_rounded_table': 'L2_literature',
               'reported_not_recomputed': 'L2_literature'}

ASSUMED = [
    ('ASSUME_obligate_share', 'obligate client flux share log-uniform [0.005,0.35]', 'model/model.py draw()'),
    ('ASSUME_essential_share', 'essential flux fraction within obligate uniform [0.02,0.4]', 'model/model.py draw()'),
    ('ASSUME_aggregate_share', 'aggregate flux share log-uniform [0.001,0.1]', 'model/model.py draw()'),
    ('ASSUME_overlap_share', 'upstream-overlap flux share uniform [0.05,0.45]', 'model/model.py draw()'),
    ('ASSUME_catalytic_multiplier', 'pool-to-capacity multipliers log-uniform [0.1,10] (sensitivity [0.01,100])', 'model/model.py draw()'),
    ('ASSUME_service_cost', 'per machine/client service cost log-uniform [0.5,2]', 'model/model.py draw()'),
    ('ASSUME_sink_capacity', 'generic protease and storage capacity 0.05 each', 'model/model.py draw()'),
    ('ASSUME_stage_atp', 'stage ATP coefficients TF 0, DnaK 1, GroEL 7, ClpB 2, protease 4', 'model/model.py architecture()'),
    ('ASSUME_serial_upstream', 'obligate clients require TF or KJE service upstream of GroEL', 'model/model.py architecture(); contradicted, see model/CONFRONTATION.md'),
    ('ASSUME_clpb_only_recovery', 'aggregate clients recovered only by ClpB with KJE', 'model/model.py architecture()'),
    ('ASSUME_radial_normalization', 'demand set to 0.9x the all-client productive radial margin', 'model/run.py; see model/CONFRONTATION.md'),
    ('ASSUME_overlap_no_groel', 'overlap clients have no GroEL alternative', 'model/model.py architecture()'),
]

# Direct reanalysis outputs produced in this run, outside evidence.tsv. Hashes of
# every input are in codon_bridge/source_hashes.json.
RECOMPUTED = [
    ('LAND_zero_structure', 'L2016', 'detection_convention',
     '4861 of 4880 covered cells; 3015 covered zeros dropped to reproduce published Data_S2 (max error 2.6e-9)'),
    ('LAND_usage_weighted', 'L2016', 'per_codon_error_rate',
     'usage-weighted E. coli rate 6.33e-4 conditional, 2.53e-4 zero-inclusive, 1.41e-4 pooled'),
    ('LAND_agc_contrast', 'L2016', 'synonymous_error_contrast',
     'Ser AGC above TCC/TCT by 12-17x pooled; depth-balanced Holm p 1e-5 and 1.1e-4; 96% of AGC errors are Ser->Asn'),
    ('LAND_sign_calibration', 'L2016', 'test_calibration',
     'all-data sign test family-wise false-hit rate 1.00 (500 Poisson simulations); balanced subset 0.01'),
    ('STIK_rate', 'S2026', 'per_codon_error_rate',
     'stationary phase 1.82e-3 +/- 0.06e-3 over 3 biological replicates; 0 Holm hits'),
    ('MODEL_serial_block', 'NONE', 'model_structural_fact',
     'TF+DnaK capacity zero gives zero essential flux in 200/200 draws at up to 100x GroEL'),
    ('MODEL_radial_inflation', 'NONE', 'model_structural_fact',
     'dropping the aggregate class raises the radial margin by a median 5.74x; >2x in 138/200 draws'),
]

def main():
    rows = list(csv.DictReader((HERE / 'evidence.tsv').open(), delimiter='\t'))
    out = []
    for r in rows:
        level = BY_ANALYSIS[r['analysis_status']]
        note = ''
        if level == 'L2_literature' and 'not_retrieved' in r['provenance_status'] or 'blocked' in r['provenance_status']:
            note = 'original table or supplement not retrieved'
        if r['analysis_status'] == 'recomputed_from_published_rounded_table':
            note = 'arithmetic on published rounded values; no replicates or uncertainty'
        out.append(dict(evidence_id=r['evidence_id'], source_id=r['source_id'],
                        evidence_class=r['evidence_class'], subject=r['subject'],
                        level=level, basis=r['analysis_status'], limit=note))
    for eid, sid, cls, subject in RECOMPUTED:
        out.append(dict(evidence_id=eid, source_id=sid, evidence_class=cls, subject=subject,
                        level='L1_direct_data',
                        basis='recomputed_in_this_run (codon_bridge/reanalyse.py, detection_controls.py, model/confront.py)',
                        limit='no folding, aggregation or growth readout in these inputs'))
    for eid, subject, where in ASSUMED:
        out.append(dict(evidence_id=eid, source_id='NONE', evidence_class='model_assumption',
                        subject=subject, level='L3_assumed', basis=where,
                        limit='not measured; no source supports this value'))
    p = HERE / 'evidence_levels.tsv'
    with p.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter='\t')
        w.writeheader(); w.writerows(out)
    counts = {}
    for r in out:
        counts[r['level']] = counts.get(r['level'], 0) + 1
    print(p, counts)

if __name__ == '__main__':
    main()
