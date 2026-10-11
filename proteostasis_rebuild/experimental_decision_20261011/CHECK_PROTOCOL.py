#!/usr/bin/env python3
"""Offline structural/provenance checks only. No fitting, simulation or assay evaluation."""
import csv
import hashlib
import itertools
from pathlib import Path
import re
import sys
ROOT = Path(__file__).resolve().parent
PROTECTED = {'EXPERIMENTAL_DECISION_REPORT.md': 'ab135311b3ecfd2055c6332aec3332204b0af0298afd4f38519748d83d209652', 'METHODS_SOURCE_AUDIT.md': 'ef75ff4fb8c17c6559b7d0dff217152d0e90429d57978e700f71945b1a53d839', 'STAGE1_STATUS.md': '23d69bae27a7b430ee8d35015ecfe299d77eab722b2a745a464ff7d39eb95869', 'SYSTEM_SELECTION.md': '9ffe2cbac162c86d65f783cda6464a734c8d89cc335c43eee9fc3a94038392e1', 'PREEXISTING_HASHES.tsv': '74dbdffc2924471d004f8545002fad631632701336064c939f33afdf761fb8b3'}
PHYSICAL_PARAGRAPHS = ['For culture c, let Bc(t) be background-subtracted, recovery-standardized heavy reporter recovered with DnaK, divided by its labelled input amount at chase zero. Define **Yc = trapezoidal integral from 0 to 40 min of Bc(t)**, in captured-cohort-fraction·min. Use equal cell equivalents relative to the original culture and correct growth dilution. Report raw recovery and B values; negative background-corrected estimates remain in the uncertainty calculation, not clipped. DnaK normalization corrects analytical recovery only: do not normalize away biological increases in DnaK abundance.', 'Let D(m,k) = [Y(ED,NVA,m,k)−Y(ED,NONE,m,k)]−[Y(REP,NVA,m,k)−Y(REP,NONE,m,k)]. The primary factorial estimand is **Theta = [D(HIGH,ADD)−D(LOW,ADD)]−[D(HIGH,EMPTY)−D(LOW,EMPTY)]**. Conditional hypothesis: functional extra ClpB lowers this error-attributable association integral, **Theta < −delta**. Transient transfer from aggregates can increase DnaK association; this is a test of the selected sign, not a universal prediction about ClpB. Val uses the same contrast as a secondary specificity check.', '**True occupied site-time primary outcome is not identifiable here.** IP detects co-associated clients, including possible indirect complexes; it does not count occupied substrate-binding sites or distinguish rebinding from long dwell. The selected mechanistic observable is only admissible if G3 passes. Even perfect cohort decay gives an apparent residence/clearance time, not microscopic dwell. Without calibrated binding multiplicity/capture, the only robust site-time bound is 0 to the integral of measured total available DnaK sites over the window, conditional on site-count calibration. No positive lower bound or infinite-window upper bound is claimed.']
COLUMNS = ['phase', 'block_id', 'lineage_id', 'culture_id', 'specimen_id', 'allele', 'exposure', 'cassette', 'inducer_level', 'reporter', 'time_id', 'planned_exposure_min', 'planned_chase_min', 'assay_bundle', 'medium_id', 'SOP_version', 'temperature_C', 'amino_acid_background', 'pulse_start_min', 'pulse_end_min', 'sample_role', 'allocation_seed', 'blinded_id', 'run_order', 'lineage_slot', 'exposure_dose_mM', 'IPTG_dose_uM', 'induction_duration_min', 'culture_date', 'culture_position', 'harvest_order', 'culture_vessel', 'fill_volume_mL', 'mixing_setting', 'culture_sampling_mode', 'cluster_id', 'actual_exposure_min', 'actual_chase_min', 'biomass', 'growth_rate', 'viable_count', 'morphology', 'total_withdrawal_volume_mL', 'missingness_reason', 'observation_status', 'INPUT__child_id', 'INPUT__parent_specimen_id', 'INPUT__assay', 'INPUT__aliquot_volume', 'INPUT__input_cell_equivalents', 'INPUT__collection_timestamp', 'INPUT__quench_timestamp', 'INPUT__storage_interval', 'INPUT__processing_interval', 'INPUT__extraction_batch', 'INPUT__instrument_batch', 'INPUT__raw_file', 'INPUT__sha256', 'INPUT__SOP_version', 'INPUT__standards_lot', 'INPUT__operator', 'INPUT__blinded_qc_decision', 'INPUT__reason', 'DNAK_IP__child_id', 'DNAK_IP__parent_specimen_id', 'DNAK_IP__assay', 'DNAK_IP__aliquot_volume', 'DNAK_IP__input_cell_equivalents', 'DNAK_IP__collection_timestamp', 'DNAK_IP__quench_timestamp', 'DNAK_IP__storage_interval', 'DNAK_IP__processing_interval', 'DNAK_IP__extraction_batch', 'DNAK_IP__instrument_batch', 'DNAK_IP__raw_file', 'DNAK_IP__sha256', 'DNAK_IP__SOP_version', 'DNAK_IP__standards_lot', 'DNAK_IP__operator', 'DNAK_IP__blinded_qc_decision', 'DNAK_IP__reason', 'ISOTYPE_IP__child_id', 'ISOTYPE_IP__parent_specimen_id', 'ISOTYPE_IP__assay', 'ISOTYPE_IP__aliquot_volume', 'ISOTYPE_IP__input_cell_equivalents', 'ISOTYPE_IP__collection_timestamp', 'ISOTYPE_IP__quench_timestamp', 'ISOTYPE_IP__storage_interval', 'ISOTYPE_IP__processing_interval', 'ISOTYPE_IP__extraction_batch', 'ISOTYPE_IP__instrument_batch', 'ISOTYPE_IP__raw_file', 'ISOTYPE_IP__sha256', 'ISOTYPE_IP__SOP_version', 'ISOTYPE_IP__standards_lot', 'ISOTYPE_IP__operator', 'ISOTYPE_IP__blinded_qc_decision', 'ISOTYPE_IP__reason', 'SOLUBLE__child_id', 'SOLUBLE__parent_specimen_id', 'SOLUBLE__assay', 'SOLUBLE__aliquot_volume', 'SOLUBLE__input_cell_equivalents', 'SOLUBLE__collection_timestamp', 'SOLUBLE__quench_timestamp', 'SOLUBLE__storage_interval', 'SOLUBLE__processing_interval', 'SOLUBLE__extraction_batch', 'SOLUBLE__instrument_batch', 'SOLUBLE__raw_file', 'SOLUBLE__sha256', 'SOLUBLE__SOP_version', 'SOLUBLE__standards_lot', 'SOLUBLE__operator', 'SOLUBLE__blinded_qc_decision', 'SOLUBLE__reason', 'PELLET__child_id', 'PELLET__parent_specimen_id', 'PELLET__assay', 'PELLET__aliquot_volume', 'PELLET__input_cell_equivalents', 'PELLET__collection_timestamp', 'PELLET__quench_timestamp', 'PELLET__storage_interval', 'PELLET__processing_interval', 'PELLET__extraction_batch', 'PELLET__instrument_batch', 'PELLET__raw_file', 'PELLET__sha256', 'PELLET__SOP_version', 'PELLET__standards_lot', 'PELLET__operator', 'PELLET__blinded_qc_decision', 'PELLET__reason', 'ACTIVITY__child_id', 'ACTIVITY__parent_specimen_id', 'ACTIVITY__assay', 'ACTIVITY__aliquot_volume', 'ACTIVITY__input_cell_equivalents', 'ACTIVITY__collection_timestamp', 'ACTIVITY__quench_timestamp', 'ACTIVITY__storage_interval', 'ACTIVITY__processing_interval', 'ACTIVITY__extraction_batch', 'ACTIVITY__instrument_batch', 'ACTIVITY__raw_file', 'ACTIVITY__sha256', 'ACTIVITY__SOP_version', 'ACTIVITY__standards_lot', 'ACTIVITY__operator', 'ACTIVITY__blinded_qc_decision', 'ACTIVITY__reason', 'STRESS__child_id', 'STRESS__parent_specimen_id', 'STRESS__assay', 'STRESS__aliquot_volume', 'STRESS__input_cell_equivalents', 'STRESS__collection_timestamp', 'STRESS__quench_timestamp', 'STRESS__storage_interval', 'STRESS__processing_interval', 'STRESS__extraction_batch', 'STRESS__instrument_batch', 'STRESS__raw_file', 'STRESS__sha256', 'STRESS__SOP_version', 'STRESS__standards_lot', 'STRESS__operator', 'STRESS__blinded_qc_decision', 'STRESS__reason', 'INPUT__peptide_sequence', 'INPUT__protein_accession', 'INPUT__substitution', 'INPUT__site', 'INPUT__isotope_state', 'INPUT__modification', 'INPUT__localization_score', 'INPUT__FDR', 'INPUT__peak_area', 'INPUT__response_factor', 'INPUT__precursor_enrichment', 'INPUT__amount', 'INPUT__uncertainty', 'INPUT__LOD', 'INPUT__LOQ', 'INPUT__censoring_flag', 'INPUT__error_count', 'INPUT__cognate_denominator', 'DNAK_IP__peptide_sequence', 'DNAK_IP__protein_accession', 'DNAK_IP__substitution', 'DNAK_IP__site', 'DNAK_IP__isotope_state', 'DNAK_IP__modification', 'DNAK_IP__localization_score', 'DNAK_IP__FDR', 'DNAK_IP__peak_area', 'DNAK_IP__response_factor', 'DNAK_IP__precursor_enrichment', 'DNAK_IP__amount', 'DNAK_IP__uncertainty', 'DNAK_IP__LOD', 'DNAK_IP__LOQ', 'DNAK_IP__censoring_flag', 'DNAK_IP__error_count', 'DNAK_IP__cognate_denominator', 'ISOTYPE_IP__peptide_sequence', 'ISOTYPE_IP__protein_accession', 'ISOTYPE_IP__substitution', 'ISOTYPE_IP__site', 'ISOTYPE_IP__isotope_state', 'ISOTYPE_IP__modification', 'ISOTYPE_IP__localization_score', 'ISOTYPE_IP__FDR', 'ISOTYPE_IP__peak_area', 'ISOTYPE_IP__response_factor', 'ISOTYPE_IP__precursor_enrichment', 'ISOTYPE_IP__amount', 'ISOTYPE_IP__uncertainty', 'ISOTYPE_IP__LOD', 'ISOTYPE_IP__LOQ', 'ISOTYPE_IP__censoring_flag', 'ISOTYPE_IP__error_count', 'ISOTYPE_IP__cognate_denominator', 'SOLUBLE__peptide_sequence', 'SOLUBLE__protein_accession', 'SOLUBLE__substitution', 'SOLUBLE__site', 'SOLUBLE__isotope_state', 'SOLUBLE__modification', 'SOLUBLE__localization_score', 'SOLUBLE__FDR', 'SOLUBLE__peak_area', 'SOLUBLE__response_factor', 'SOLUBLE__precursor_enrichment', 'SOLUBLE__amount', 'SOLUBLE__uncertainty', 'SOLUBLE__LOD', 'SOLUBLE__LOQ', 'SOLUBLE__censoring_flag', 'SOLUBLE__error_count', 'SOLUBLE__cognate_denominator', 'PELLET__peptide_sequence', 'PELLET__protein_accession', 'PELLET__substitution', 'PELLET__site', 'PELLET__isotope_state', 'PELLET__modification', 'PELLET__localization_score', 'PELLET__FDR', 'PELLET__peak_area', 'PELLET__response_factor', 'PELLET__precursor_enrichment', 'PELLET__amount', 'PELLET__uncertainty', 'PELLET__LOD', 'PELLET__LOQ', 'PELLET__censoring_flag', 'PELLET__error_count', 'PELLET__cognate_denominator', 'ACTIVITY__raw_signal', 'ACTIVITY__blank', 'ACTIVITY__dilution', 'ACTIVITY__linear_range', 'ACTIVITY__calibration', 'ACTIVITY__matrix_spike_recovery', 'ACTIVITY__activity_equivalents', 'MASS_BALANCE__total_recovery', 'MASS_BALANCE__soluble_recovery', 'MASS_BALANCE__pellet_recovery', 'MASS_BALANCE__cell_counts']
ASSAYS = ['INPUT', 'DNAK_IP', 'ISOTYPE_IP', 'SOLUBLE', 'PELLET', 'ACTIVITY', 'STRESS']
FREEZE_PARAMETERS = ['G0_decision', 'G1_decision', 'G2_decision', 'G3_decision', 'G4_decision', 'G5_decision', 'G6_decision', 'new_culture_provenance', 'pilot_access_log', 'confirmation_independence', 'sequence_identity', 'nascent_error_increase', 'ED_specific_error_increment', 'primary_spectrum_resolution', 'peptide_standard_identity', 'site_localization_limit', 'peptide_FDR_limit', 'response_factors', 'precursor_enrichment_limit', 'precursor_clearance_limit', 'pulse_exchange_SOP', 'late_pulse_requirement', 'error_LOD_LOQ', 'functional_capacity_increment', 'functional_ClpB_SOP', 'ClpB_abundance_plateau', 'machine_abundance_panel', 'growth_rate_margin', 'viable_fraction_loss_margin', 'synthesis_ratio_margin', 'NVA_error_frequency_margin', 'VAL_error_frequency_margin', 'spectrum_component_margins', 'REP_toxicity_margin', 'reporter_perturbation_margin', 'G2_comparison_family', 'G2_multiplicity', 'endogenous_synthesis_panel', 'capture_specificity_limit', 'capture_preservation_limit', 'capture_dead_time_limit', 'capture_linear_range', 'capture_recovery_limit', 'background_limit', 'mass_balance_limit', 'association_LOD_LOQ', 'DnaK_recovery_calibration', 'cell_equivalent_calibration', 'calibration_covariance', 'withdrawal_volume_limit', 'aliquot_volume_minima', 'sampling_mode', 'batch_QC_limits', 'NVA_dose', 'VAL_dose', 'HIGH_IPTG_dose', 'induction_duration', 'medium_formulation', 'temperature', 'aeration_vessel_fill_mixing', 'exposure_growth_window', 'pulse_timing', 'collection_grid', 'sigma_U', 'variance_bound_method', 'delta', 'planning_distance_d', 'power_and_precision', 'confirmation_n_blocks', 'feasible_n_limit', 'confirmation_lineage_balance', 'estimator_and_clustering', 'uncertainty_method', 'secondary_Holm_family', 'missing_and_censored_rule', 'technical_failure_rule', 'replacement_rule', 'randomization_seed', 'randomization_key_sha256', 'protocol_sha256', 'analysis_plan_sha256', 'source_list_sha256', 'sample_sheet_sha256', 'analysis_code_sha256', 'freeze_timestamp', 'analysts_and_signatures', 'decision_CI_rule', 'prior_gate_evidence', 'pulse_cohort_function_identifiability', 'specific_activity_calibration', 'activity_linearity_substrate', 'activity_matrix_spike_recovery', 'old_and_postpulse_reporter_contributions', 'new_native_yield_increase', 'competing_synthesis_loss_bound', 'competing_degradation_loss_bound', 'competing_death_loss_bound', 'competing_recovery_loss_bound', 'synthesis_stop_physiology_validation']

RESULTS = []
def check(name, condition, detail=''):
    passed = bool(condition)
    RESULTS.append(passed)
    print(('PASS' if passed else 'FAIL') + '\t' + name + ('\t' + str(detail) if detail else ''))

def table(name):
    with (ROOT / name).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = list(reader)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(name + ': missing/duplicate headers')
        if any(None in row or any(v is None or v == '' for v in row.values()) for row in rows):
            raise ValueError(name + ': ragged row or blank value')
        return reader.fieldnames, rows

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    required = ['IDENTIFIABILITY_MATRIX.tsv', 'SAMPLE_SHEET_TEMPLATE.tsv',
                'FREEZE_RECORD_TEMPLATE.tsv', 'GATE_SPEC.tsv', 'SOURCE_ANCHORS.tsv',
                'SELF_AUDIT.md', 'CHECK_PROTOCOL.py']
    check('required_files', all((ROOT / f).is_file() for f in required), len(required))
    fields, rows = table('SAMPLE_SHEET_TEMPLATE.tsv')
    check('sample_metadata_schema', fields == COLUMNS, str(len(fields)) + ' columns')
    check('sample_count', len(rows) == 672, str(len(rows)) + ' data rows; header excluded')
    blocks = ['P01', 'P02', 'P03', 'C_BLOCK_REPEAT']
    check('block_ids', {r['block_id'] for r in rows} == set(blocks), '4 blocks')
    combinations = set(itertools.product(['ED', 'REP'], ['NONE', 'NVA', 'VAL'], ['ADD', 'EMPTY'], ['LOW', 'HIGH']))
    factors = ['allele', 'exposure', 'cassette', 'inducer_level']
    times = {'B0': ('-1', 'NA'), 'C00': ('25', '0'), 'C02': ('27', '2'),
             'C05': ('30', '5'), 'C10': ('35', '10'), 'C20': ('45', '20'), 'C40': ('65', '40')}
    for block in blocks:
        subset = [r for r in rows if r['block_id'] == block]
        observed = [(tuple(r[f] for f in factors), r['time_id']) for r in subset]
        expected = set(itertools.product(combinations, times))
        check('balance_' + block, len(subset) == 168 and len(set(observed)) == 168 and set(observed) == expected,
              str(len(subset)) + ' rows; 24 factor combinations x 7 times')
    check('time_mapping', all(times.get(r['time_id']) == (r['planned_exposure_min'], r['planned_chase_min']) for r in rows), 'B0 plus six chase times')
    check('specimen_unique', len({r['specimen_id'] for r in rows}) == 672, '672 IDs')
    culture_groups = {}
    for r in rows:
        culture_groups.setdefault(r['culture_id'], []).append(r)
    check('cultures_and_repeated_sampling', len(culture_groups) == 96 and all(len(v) == 7 and len({r['time_id'] for r in v}) == 7 for v in culture_groups.values()), '96 cultures x 7 repeated specimens')
    child_ids = [r[a + '__child_id'] for r in rows for a in ASSAYS]
    all_ids = list(culture_groups) + [r['specimen_id'] for r in rows] + child_ids
    check('global_id_uniqueness', len(all_ids) == len(set(all_ids)) == 5472, '96 culture + 672 specimen + 4704 child IDs')
    ids_ok = True
    links_ok = True
    for r in rows:
        cid = 'EXP20261011_' + '_'.join([r['block_id']] + [r[f] for f in factors])
        sid = cid + '_' + r['time_id']
        ids_ok &= (r['culture_id'] == cid and r['specimen_id'] == sid)
        for a in ASSAYS:
            links_ok &= (r[a + '__child_id'] == sid + ':' + a and r[a + '__parent_specimen_id'] == sid and r[a + '__assay'] == a)
    check('stable_factor_ids', ids_ok)
    check('child_pairing', links_ok, '7 matched assay children per parent')
    check('phase_and_lineage_slots', all(r['phase'] == ('CONFIRMATION_PLACEHOLDER' if r['block_id'] == 'C_BLOCK_REPEAT' else 'PILOT') and r['lineage_id'] == 'PENDING' and r['lineage_slot'] == ('PENDING' if r['block_id'] == 'C_BLOCK_REPEAT' else r['block_id']) for r in rows), 'actual lineage identities unknown')
    known = set('phase block_id culture_id specimen_id allele exposure cassette inducer_level reporter time_id planned_exposure_min planned_chase_min assay_bundle temperature_C pulse_start_min pulse_end_min sample_role lineage_slot exposure_dose_mM IPTG_dose_uM'.split())
    known |= {a + '__' + f for a in ASSAYS for f in ['child_id', 'parent_specimen_id', 'assay']}
    unknown = set(fields) - known
    check('no_fake_observations', all(r[f] == 'PENDING' for r in rows for f in unknown), str(len(unknown)) + ' unknown fields per row explicitly PENDING')
    check('planned_doses_only', all(r['exposure_dose_mM'] == ('0' if r['exposure'] == 'NONE' else 'PENDING') and r['IPTG_dose_uM'] == ('0' if r['inducer_level'] == 'LOW' else 'PENDING') for r in rows), 'zeros are known design doses only')
    check('shared_design_settings', all(r['temperature_C'] == '30' and r['pulse_start_min'] == '20' and r['pulse_end_min'] == '25' and r['assay_bundle'] == ';'.join(ASSAYS) and r['reporter'] == 'single_copy_constitutive_firefly_luciferase' and r['sample_role'] == ('BASELINE_QC' if r['time_id'] == 'B0' else 'CHASE') for r in rows), 'planned settings; not measured')
    check('NA_only_baseline_chase', all(v != 'NA' or (k == 'planned_chase_min' and r['time_id'] == 'B0') for r in rows for k, v in r.items()))
    _, matrix = table('IDENTIFIABILITY_MATRIX.tsv')
    check('identifiability_matrix', len(matrix) == 12 and len({r['construct'] for r in matrix}) == 12, '12 distinct constructs')
    check('unidentified_constructs_explicit', all('UNIDENTIFIABLE' in next(r for r in matrix if r['construct'] == c)['prospective_identifiability'] for c in ['true_occupancy_site_time', 'microscopic_dwell_flux', 'pulse_cohort_functional_yield']))
    protocol = (ROOT / 'ONE_EXPERIMENT_PROTOCOL.md').read_text()
    plan = (ROOT / 'ANALYSIS_PLAN.md').read_text()
    _, gates = table('GATE_SPEC.tsv')
    check('gate_ids', [r['gate_id'] for r in gates] == ['G' + str(n) for n in range(7)], 'G0-G6 once each')
    extracted = re.findall(r'^\| (G[0-6]) (.*?) \| (.*?) \|$', protocol, re.M)
    check('gate_rules_verbatim', [(r['gate_id'], r['gate_name'], r['exact_rule']) for r in gates] == extracted, 'exact match to corrected protocol')
    check('no_gate_pass_fabricated', all(r['status'] == 'PENDING_LOCK' for r in gates))
    ffields, freeze = table('FREEZE_RECORD_TEMPLATE.tsv')
    check('freeze_parameter_coverage', [r['parameter'] for r in freeze] == FREEZE_PARAMETERS and len(FREEZE_PARAMETERS) == len(set(FREEZE_PARAMETERS)) and {r['gate_id'] for r in freeze} == {'G' + str(n) for n in range(7)}, str(len(freeze)) + ' gate/calibration/lock records')
    lock_fields = ['value', 'evidence_file', 'evidence_sha256', 'locked_at_utc', 'locked_by', 'approval_signature']
    check('freeze_unlocked', all(r[f] == 'PENDING_LOCK' for r in freeze for f in lock_fields), 'all values/evidence/signatures PENDING_LOCK')
    check('freeze_gate_rules', all(next(r for r in freeze if r['parameter'] == g['gate_id'] + '_decision')['requirement'] == g['exact_rule'] for g in gates))
    check('physical_correction_both_documents', all('UNIDENTIFIABLE' in text and 'unlabeled' in text.lower() and 'specific activity' in text and 'conditional synthesis-stop preserving physiology' in text for text in [protocol, plan]) and 'max(0,A−O)' not in protocol and 'min(A,N)' not in protocol)
    g6 = next(r['exact_rule'] for r in gates if r['gate_id'] == 'G6')
    check('G6_identifiability_required', 'no new-native yield claim' in g6 and 'UNIDENTIFIABLE' in g6 and 'specific activity' in g6)
    check('other_physical_equations_preserved', all(p in protocol for p in PHYSICAL_PARAGRAPHS) and 'n >= ((z_0.975 + z_0.90) sigma_U / d)^2' in plan, 'Y, D, Theta, site-time bound, power formula')
    _, sources = table('SOURCE_ANCHORS.tsv')
    check('source_ids', [r['source_id'] for r in sources] == ['S' + str(n) for n in range(1, 12)], 'S1-S11')
    for r in sources:
        p = Path(r['local_path'])
        start, end = int(r['line_start']), int(r['line_end'])
        lines = p.read_text().splitlines() if p.is_file() else []
        valid = 1 <= start <= end <= len(lines)
        anchor_hash = hashlib.sha256('\n'.join(lines[start-1:end]).encode()).hexdigest() if valid else ''
        check('source_link_' + r['source_id'], valid and anchor_hash == r['anchor_sha256'], str(p) + ':' + str(start) + '-' + str(end))
    doi_map = {r['source_id']: r['doi'] for r in sources}
    check('literature_identifiers', doi_map['S4'] == '10.1016/j.ijbiomac.2024.130068' and doi_map['S9'] == '10.1002/1873-3468.13578' and doi_map['S10'] == '10.1002/pro.3639' and doi_map['S11'] == '10.1093/nar/gkaf428' and all(r['url'] == ('NA' if r['doi'] == 'NA' else 'https://doi.org/' + r['doi']) for r in sources), 'DOI links recorded; no live URL access asserted')
    refs = set()
    for name in ['EXPERIMENTAL_DECISION_REPORT.md', 'ONE_EXPERIMENT_PROTOCOL.md', 'STAGE1_STATUS.md', 'SYSTEM_SELECTION.md']:
        content = (ROOT / name).read_text()
        refs |= set(re.findall(r'\bS\d+\b', content))
        for lo, hi in re.findall(r'S(\d+)[–-]S(\d+)', content):
            refs |= {'S' + str(n) for n in range(int(lo), int(hi) + 1)}
    check('document_source_refs_resolve', refs <= set(doi_map), ','.join(sorted(refs, key=lambda s: int(s[1:]))))
    failures = []
    with (ROOT / 'PREEXISTING_HASHES.tsv').open(newline='') as handle:
        preserved = list(csv.DictReader(handle, delimiter='\t'))
    for r in preserved:
        p = Path(r['path'])
        if not p.is_file() or p.stat().st_size != int(r['bytes']) or digest(p) != r['sha256']:
            failures.append(str(p))
    check('preexisting_hash_preservation', not failures, str(len(preserved) - len(failures)) + '/' + str(len(preserved)) + ' byte-identical files' + ('; failures=' + repr(failures) if failures else ''))
    check('protected_local_documents', all(digest(ROOT / name) == value for name, value in PROTECTED.items()), str(len(PROTECTED)) + ' unchanged documents/manifests')
    report = (ROOT / 'EXPERIMENTAL_DECISION_REPORT.md').read_text()
    check('prior_verdict_retained', 'VOID' in report and 'no confirmatory new test from currently identified public processed data' in report and 'no independent holdout is established' in report.lower())
    check('readiness_explicit', all('no wet lab, confirmation not ready' in (ROOT / name).read_text() for name in ['README.md', 'SELF_AUDIT.md']))
    broken = []
    for name in ['README.md', 'EXPERIMENTAL_DECISION_REPORT.md']:
        for target in re.findall(r'\]\(([^)]+)\)', (ROOT / name).read_text()):
            if '://' not in target and target != 'CHECK_PROTOCOL_RESULTS.txt' and not (ROOT / target.split('#')[0]).exists():
                broken.append(target)
    check('local_delivery_links', not broken, repr(broken) if broken else 'all targets present; results file generated by this run')

if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        check('checker_execution', False, type(exc).__name__ + ': ' + str(exc))
    passed = sum(RESULTS)
    failed = len(RESULTS) - passed
    print('TOTAL\tPASS=' + str(passed) + '\tFAIL=' + str(failed))
    print('no wet lab, confirmation not ready')
    sys.exit(1 if failed else 0)
