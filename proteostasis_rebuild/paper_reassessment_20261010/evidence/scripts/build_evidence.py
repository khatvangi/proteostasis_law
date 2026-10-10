"""Build explicit evidence records and limited comparisons; no simulated observations."""
import csv, hashlib, json, sys
from pathlib import Path
import xml.etree.ElementTree as ET
sys.dont_write_bytecode = True
OUT = Path(__file__).resolve().parents[1]
PROJECT = OUT.parents[1]

def tsv(name, rows):
    assert rows
    with (OUT/name).open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

sources={
 'K2005': ['10.1016/j.cell.2005.05.028','PMID:16051146;PII:S0092867405005441','https://pubmed.ncbi.nlm.nih.gov/16051146/','abstract_verified;original_table_not_retrieved','MC4100 derivatives (supplement search-index); assay-specific medium/temperature not verified'],
 'C2012': ['10.1016/j.celrep.2011.12.007','PMID:22832197;PII:S2211124711000179','https://pure.mpg.de/rest/items/item_1602419_3/component/file_1836221/content','main_PDF_web_text_verified;original_supplement_tables_not_retrieved','MC4100 derivatives; conditions specified per row'],
 'B2014': ['10.1371/journal.pgen.1004120','PMCID:PMC3930520;PMID:24586182','https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004120','main_article_verified;XLS_not_retrieved','Hfr-C donor; Keio deletion and hypomorphic recipients; LB+Kan+Cm;32 C;36 h'],
 'F2010': ['10.1038/emboj.2010.52','PMCID:PMC3212837;PMID:20360681','https://pmc.ncbi.nlm.nih.gov/articles/PMC3212837/','main_article_web_verified','MGM100 (MG1655 derivative);LB+DAP;2 h GroE depletion;glucose versus arabinose'],
 'G1989': ['10.1128/jb.171.3.1379-1385.1989','PMCID:PMC209756;PMID:2563997','https://pmc.ncbi.nlm.nih.gov/articles/PMC209756/?page=2','primary_abstract_verified','E. coli groE heterodiploids;17-30 C; exact strain/medium unresolved'],
 'E1989': ['10.1128/jb.171.5.2748-2755.1989','PMCID:PMC209960;PMID:2651417','https://pubmed.ncbi.nlm.nih.gov/2651417/','primary_abstract_verified','E. coli;grpE deletion/complementation; exact strain/medium/temp panel unresolved'],
 'H2021': ['10.3389/fmolb.2021.653073','PXD010014','https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2021.653073/full','main_article_web_verified','W3110;LB;30 C permissive for dnaKJ deletion;37 C for other mutant precultures'],
 'I2015': ['10.3389/fmicb.2015.00051','PMCID:PMC4316790;PMID:25699035','https://pmc.ncbi.nlm.nih.gov/articles/PMC4316790/','primary_article_search_index_verified;direct_open_challenge','MG1655;LB;37 C;nlpI plasmid with/without 0.5 mM IPTG'],
 'P1999': ['10.1128/jb.181.12.3674-3680.1999','PMCID:PMC93843;PMID:10368140','https://pmc.ncbi.nlm.nih.gov/articles/PMC93843/','Table_1_and_2_primary_article_search_index_verified;direct_open_challenge','MG1655 derivatives;L broth preculture 42 C;L agar plating;20-45 C'],
 'L1991': ['10.1128/jb.173.14.4254-4262.1991','PMID:2066329','https://pubmed.ncbi.nlm.nih.gov/2066329/','primary_abstract_verified','E. coli clpB null;44 C growth;50 C survival;exact background unresolved'],
 'D2004': ['10.1074/jbc.M313635200','PMID:14729669','https://d-nb.info/1106330919/34','primary_methods_search_index_verified;PDF_download_blocked','MC4100 GK2;complete dnaK ORF deletion;DnaJ retained;LB;precise permissive temperature unresolved'],
 'S2016': ['10.1038/nbt.3418','PMCID:PMC4888949;PMID:26641532','https://pmc.ncbi.nlm.nih.gov/articles/PMC4888949/bin/NIHMS65833-supplement-Supplementary_tables.xlsx','local_original_workbook_checksum_verified','BW25113;glucose exponential;Schmidt S6 and calculated volume S23'],
}

manifest=[]
for key,(doi,acc,url,status,conditions) in sources.items():
    manifest.append(dict(source_id=key,doi=doi,accession_ids=acc,exact_url=url,artifact='article_or_workbook',labels='see evidence.tsv',strain_conditions=conditions,retrieval_status=status,local_artifact='WO-06/records/schmidt2016/NIHMS65833-supplement-Supplementary_tables.xlsx' if key=='S2016' else 'sources/*_web.json and related tool snapshots',sha256=sha(PROJECT/'WO-06/records/schmidt2016/NIHMS65833-supplement-Supplementary_tables.xlsx') if key=='S2016' else '',notes='No dataset accession invented; article IDs are not raw-data accessions'))
for label,suffix in [('Table S1','s022'),('Table S2','s023'),('Table S6','s027'),('Table S7','s028'),('Table S8','s029'),('Table S9','s030'),('Table S10','s031'),('Table S16','s037'),('Protocol S2','s007'),('Protocol S16','s021')]:
    manifest.append(dict(source_id='B2014_'+suffix,doi='10.1371/journal.pgen.1004120.'+suffix,accession_ids='PMCID:PMC3930520',exact_url='https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1004120.'+suffix+'&type=supplementary',artifact='publisher_supplement',labels=label,strain_conditions=sources['B2014'][4],retrieval_status='link_identity_verified;bytes_not_retrieved',local_artifact='',sha256='',notes='S2=high-confidence pairs;S9=profile correlation;S10=GI overlap, not physical clients;unattempted individual URLs remain discovery only'))
manifest.extend([
 dict(source_id='K2005_SUPP',doi=sources['K2005'][0],accession_ids='CiteSeerX document hash:4087a93abb0c94d87af0e09c2c904f50c46c3d48',exact_url='https://citeseerx.ist.psu.edu/document?doi=4087a93abb0c94d87af0e09c2c904f50c46c3d48&repid=rep1&type=pdf',artifact='original_supplement_mirror_candidate',labels='Supplemental Experimental Procedures;Table S1;client classification Table S3 (secondary pointer only)',strain_conditions=sources['K2005'][4],retrieval_status='indexed_methods_only;download_failed',local_artifact='',sha256='',notes='Mirror title/authors match; cannot authenticate bytes or reconstruct clients'),
 dict(source_id='C2012_SUPP',doi=sources['C2012'][0],accession_ids='PII:S2211124711000179',exact_url='https://doi.org/10.1016/j.celrep.2011.12.007',artifact='supplement_landing',labels='S1 interactors;S2 abundance;S3 aggregation;S6 TF loss;S7 GroE depletion;S8 combined loss',strain_conditions=sources['C2012'][4],retrieval_status='labels_verified_from_article;direct_table_URL_unresolved',local_artifact='',sha256='',notes='Do not label guessed mmc filenames as verified supplement URLs'),
 dict(source_id='C2012_TRANCHE',doi=sources['C2012'][0],accession_ids='RKuFJcuu8iBnZHlBN4Nth0pz+HQgPNvV0zJvdnMN8wyAN8l7ifufEnmW6j2Cwn0msEakFC6eulEqlYv+7B+dFALonrgAAAAAAAAKDA==',exact_url='https://proteomecommons.org/',artifact='legacy_MaxQuant_repository',labels='Experimental Procedures p262',strain_conditions=sources['C2012'][4],retrieval_status='article_identifier_transcribed;repository_not_accessible',local_artifact='',sha256='',notes='Published public passphrase: vqXGUF93rGLCraqm1yIl;line-wrap reconstruction not validated against repository;raw spectra have additional codes in unavailable extended procedures'),
])
tsv('source_manifest.tsv',manifest)

ev=[]
def add(eid,src,cls,subject,label,metric,value,unit,source_fact,interpretation,conditions=None,confidence='medium',operator='=',status='reported_not_recomputed'):
    doi,acc,url,access,cond=sources[src]
    ev.append(dict(evidence_id=eid,source_id=src,evidence_class=cls,subject=subject,strain_conditions=conditions or cond,figure_table=label,metric=metric,value=value,operator=operator,unit=unit,source_fact=source_fact,interpretation=interpretation,doi=doi,accession_ids=acc,exact_url=url,provenance_status=access,analysis_status=status,confidence=confidence))

add('K01','K2005','physical_interactome','GroEL','Abstract','client_identities',250,'proteins','GroEL-associated proteins','Binding is not obligate dependence or demand flux',operator='approx')
add('K02','K2005','inferred_obligate_dependence','GroEL','Abstract','classified_obligate_candidates',85,'proteins','Original dependence estimate','Not 85 individually validated in-vivo dependencies',operator='approx')
add('K03','K2005','client_essentiality_annotation','GroEL clients','Abstract','essential_clients_in_original_class',13,'proteins','Original essential-client annotation','Not 13 essential chaperone machines')
add('K04','K2005','occupancy_capacity_share_estimate','GroEL','Abstract','reported_capacity_share',75,'percent','Share attributed to obligate class','No time denominator; not folding throughput',operator='>')
add('C01','C2012','physical_interactome','DnaK','S1;Fig1','interactors',674,'proteins','674 with cofactors;672 excluding DnaJ/GrpE','Co-isolation;not universal direct binary binding',conditions='MC4100 KHis;M63 SILAC;37 C;exponential OD600~1')
add('C02','C2012','perturbed_physical_interactome','TF deletion','S6A/S6B;Fig4','increased_association',71,'proteins','71 increased;74 decreased','Rerouting interpretation;not flux',conditions='MC4100 delta-tig/KHis;37 C;M63 SILAC')
add('C03','C2012','perturbed_physical_interactome','GroE depletion','S7A/S7B;Fig5','increased_association',92,'proteins','92 increased;54 decreased','Retention does not establish productive folding',conditions='MC4100 LS-/KHis;37 C;3.5 h glucose shift;~97% depletion')
add('C04','C2012','perturbation_aggregation','dnaKJ versus dnaKJ tig','S8;Fig6A','aggregating_proteins_dnaKJ',474,'proteins','Combined deletion:1087;delta-tig:15','Count contrasts;not biomass or flux',conditions='MC4100 derivatives;30 C;SILAC')
add('C05','C2012','kinetic_association_not_absolute_flux','DnaK','Fig3;Experimental Procedures p262','chase_times','2;4;8','min','Pulse-chase association measurements','No matched molecules/cell/s client demand recovered',conditions='MC4100 KHis;M63 SILAC;37 C')
add('B01','B2014','deletion_growth_genetic_interaction','22 stability factors','S2;Fig6','high_confidence_GIs',3816,'gene_pairs','At least one of22 factors;|S|>=3;P<=0.05','Not 3816 machine-machine edges;not binding',confidence='high')
add('F01','F2010','tested_in_vivo_dependence','Kerner Class III','Fig2C;S3/S4','dependent',49,'proteins','49 dependent;34 independent;Ypt1 untested','ClassIII is not uniformly obligatory',confidence='high')
add('H01','H2021','deletion_growth_suppression','dnaKJ htpG','Fig1A','test_temperatures','37;39','C','htpG removal improves dnaKJ-mutant growth','Competes with monotone capacity-only deletion model',conditions='W3110;LB agar;overnight spots',confidence='high')

machines=[('TF','tig','P0A850'),('DnaK','dnaK','P0A6Y8'),('DnaJ','dnaJ','P08622'),('GrpE','grpE','P09372'),('GroEL','groL','P0A6F5'),('GroES','groS','P0A6F9'),('HtpG','htpG','P0A6Z3'),('ClpB','clpB','P63284'),('IbpA','ibpA','P0C054'),('IbpB','ibpB','P0C058'),('Lon','lon','P0A9M0'),('ClpX','clpX','P0A6H1'),('ClpP','clpP','P0A6G7'),('HslU','hslU','P0A6H5'),('HslV','hslV','P0A7B8')]
# Each row states exactly which perturbation was observed. UNKNOWN is not viable.
calls={
 'tig':('C2012','viable','single_deletion','Fig1B;Fig4E','delta-tig/KHis','MC4100;M63;37 C','Functional tagged DnaK background'),
 'dnaK':('D2004','viable','single_deletion','Experimental Procedures:Strains and Plasmids','GK2 delta-dnaK','MC4100;LB;temperature unresolved','DnaJ expression retained;growth effect not quantified'),
 'dnaJ':('H2021','UNKNOWN_single_deletion','combined_deletion','Fig1A','delta-dnaK dnaJ','W3110;LB;30 C','Combined mutant grows;cannot assign clean single-deletion call'),
 'grpE':('E1989','essential_without_suppressor','deletion_complementation','Abstract;genetic replacement experiment','grpE deletion/omega marker','E. coli;temperature panel and strain unresolved','Viable with complement or extragenic suppression;not universal lethality'),
 'groL':('G1989','essential_in_tested_background','inactivation_transduction','Abstract;heterodiploid experiments','groEL inactivation','E. coli;17-30 C;strain unresolved','Complement dependence;not clean Keio deletion'),
 'groS':('G1989','essential_in_tested_background','inactivation_transduction','Abstract;heterodiploid experiments','groES inactivation','E. coli;17-30 C;strain unresolved','Complement dependence;not clean Keio deletion'),
 'htpG':('H2021','viable','single_deletion','Fig1A','delta-htpG','W3110;LB;30/37 C','No major defect at30 C'),
 'clpB':('L1991','viable_growth_defect_heat','null_allele','Abstract;growth/survival experiment','clpB null','E. coli;44 C growth;50 C survival;strain unresolved','Null allele existence and growth;exact deletion boundaries not verified'),
 'ibpA':('I2015','viable','single_deletion','Table1;Fig3D','delta-ibpA','MG1655;LB;37 C','nlpI plasmid assay;not matched glucose'),
 'ibpB':('I2015','viable','single_deletion','Table1;FigS2B as described in Results','delta-ibpB','MG1655;LB;37 C','Supplement figure not independently viewed'),
 'lon':('H2021','viable','single_deletion','Methods:Bacterial Strains;Viability Assays','delta-lon','W3110;LB;37 C preculture','Constructed single mutant;not quantitative fitness'),
 'clpX':('P1999','UNKNOWN_single_deletion','combined_deletion','Table1;Table2','KY2347 delta(clpPX-lon)','MG1655;L agar;42 C reference','Combined deletion viable;single-gene inference withheld'),
 'clpP':('H2021','viable','single_deletion','Methods:Bacterial Strains;Viability Assays','delta-clpP','W3110;LB;37 C preculture','Constructed single mutant;not quantitative fitness'),
 'hslU':('P1999','UNKNOWN_single_deletion','combined_deletion','Table1;Table2','KY2966 delta-hslVU','MG1655;L agar;42 C reference','Operon deletion viable;single-gene inference withheld'),
 'hslV':('H2021','viable','single_deletion','Methods:Bacterial Strains;Viability Assays','delta-hslV','W3110;LB;37 C preculture','Constructed single mutant;not quantitative fitness'),
}
deletions=[]
for machine,gene,acc in machines:
    src,call,pert,label,allele,cond,note=calls[gene]
    doi,ids,url,access,_=sources[src]
    deletions.append(dict(machine=machine,gene=gene,aliases={'groL':'groEL','groS':'groES'}.get(gene,''),uniprot=acc,single_gene_call=call,observed_perturbation=pert,allele=allele,strain_conditions=cond,experiment_location=label,source_id=src,doi=doi,exact_url=url,provenance_status=access,confidence='medium',matched_BW25113_glucose='no',limitations=note))
    add('DEL_'+gene,src,'deletion_essentiality' if not call.startswith('UNKNOWN') else 'combined_deletion_not_single_essentiality',gene,label,'condition_specific_call',call,'categorical',allele,note,cond)
tsv('machine_deletion_essentiality.tsv',deletions)

# A real, small experimental table, transcribed from the primary article index.
# Values are relative to each strain at42 C, NOT normalized to WT at each temperature.
temps=[20,30,37,42,44,45]
matrix={
 'MG1655':([1.5,1.5,1.4,1,0.88,1.4],'WT'),
 'KY2966':([0.8,1.3,0.88,1,None,0.97],'delta-hslVU'),
 'KY2347':([0.72,0.22,0.41,1,1.5,1.6],'delta(clpPX-lon)'),
 'KY2350':([1.4e-5,3.3e-5,1e-2,1,1.5,1.6e-4],'delta-hslVU delta(clpPX-lon)')}
plate=[]
for strain,(vals,geno) in matrix.items():
    for temp,val in zip(temps,vals):
        plate.append(dict(strain=strain,genotype=geno,temperature_C=temp,relative_plating_efficiency='' if val is None else val,missing_status='ND' if val is None else '',normalization='same strain at42 C',footnote='b' if strain in ['KY2347','KY2350'] and temp in [20,30] else '',source_id='P1999',table='Table2',exact_url=sources['P1999'][2],provenance='manual_transcription_of_primary_table_via_search_index;no_raw_colony_counts'))
tsv('published_plating_table.tsv',plate)
comparison=[]
for i,t in enumerate(temps):
    wt,a,b,ab=[matrix[k][0][i] for k in ['MG1655','KY2966','KY2347','KY2350']]
    ratio=ab*wt/(a*b) if a is not None else None
    comparison.append(dict(temperature_C=t,relative_interaction_ratio='' if ratio is None else ratio,formula='EOP_AB*EOP_WT/(EOP_A*EOP_B)',interpretation='temperature contrast relative to42 C interaction;not an absolute epistasis score',status='computed_from_published_rounded_values' if ratio is not None else 'not_computable_ND',p_value='',n_replicates='',source='P1999 Table2',confound='multi-gene deletions;strain-specific42 C denominators;no uncertainty'))
tsv('plating_comparison.tsv',comparison)
add('P01','P1999','deletion_growth_genetic_interaction','combined protease deletions','Table2','relative_interaction_ratio_at30C',comparison[1]['relative_interaction_ratio'],'dimensionless','Published EOP ratio contrast','Relative to42 C;not calibrated GI S-score',confidence='medium',status='recomputed_from_published_rounded_table')

sys.path.insert(0,str(PROJECT/'network_failure_test'))
import pools
poolvalues=pools.findPools()
sys.path.insert(0,str(PROJECT/'WO-06'))
import schmidt_pools as sp
copies,mass,vol,gr=sp.load()
poolrows=[]
for machine,gene,acc in machines:
    subunits=(pools.WO07_ACC|pools.SCHMIDT_EXTRA)[machine][1]
    n=float(copies.loc[acc,'Glucose'])
    poolrows.append(dict(machine=machine,gene=gene,uniprot=acc,copies_per_cell=n,volume_fL=vol['Glucose'],subunits_per_unit=subunits,unit_pool_uM=poolvalues[machine],source_id='S2016',label='S6 Glucose;S23 BW25113 Glucose',meaning='abundance converted using calculated cell volume;not capacity or throughput'))
    add('POOL_'+gene,'S2016','abundance_not_flux',gene,'S6;S23','unit_pool',poolvalues[machine],'uM functional units','Copies/cell converted with calculated volume','Per-unit throughput absent',confidence='high',status='recomputed_from_local_original_workbook')
tsv('machine_abundance.tsv',poolrows)
tsv('evidence.tsv',ev)

with (PROJECT/'network_failure_test/results/machine_knockouts.csv').open() as f:
    model={r['machine']:r for r in csv.DictReader(f) if r['setting']=='bio_measured_W2'}
tsv('model_vs_experiment.tsv',[dict(machine=r['machine'],gene=r['gene'],model_setting='bio_measured_W2',model_fraction_structural_zero=model[r['machine']]['frac_ko_structural_zero'],model_fraction_KO_lethal_rho05=model[r['machine']]['frac_ko_lethal_rho0.5'],experimental_call=r['single_gene_call'],source_id=r['source_id'],matched_condition='no',conclusion='qualitative_constraint_only;no_accuracy_or_p_value') for r in deletions])

inventory=[]
for r in json.loads((OUT/'original_files_baseline.json').read_text()):
    p=r['path']
    if p.endswith(('.tsv','.csv','.xlsx','.xml','.json')):
        kind='synthetic_model_output' if p.startswith('network_failure_test/results/') or p.startswith('network_failure_test/audit/results/') else 'original_experimental_supplement' if p.endswith('.xlsx') else 'bibliographic_or_fulltext_record' if p.startswith('WO-06/records/') else 'derived_audit_or_parameter_data'
        inventory.append(dict(**r,evidence_role=kind))
tsv('local_data_inventory.tsv',inventory)
print(json.dumps({'evidence_rows':len(ev),'machines':len(deletions),'unresolved_single_deletions':[r['gene'] for r in deletions if r['single_gene_call'].startswith('UNKNOWN')],'plating_rows':len(plate),'ratio30C':comparison[1]['relative_interaction_ratio'],'local_data_files':len(inventory)},indent=2))
