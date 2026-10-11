#!/usr/bin/env python3
"""Stage I only: bounded retrieval, hashes, schemas, first-three-row previews.
No contrasts, fits, imputation, ranks, enrichment or hypothesis searches.
Writes only beside this script. Python standard library; curl for --fetch.
"""
import argparse
import concurrent.futures
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT.parents[1] / 'proteostasis-law-empirical-phase3'
RAW = ARCHIVE / 'data/raw/outcomes/pranjic_2024'

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()

def tsv(name, rows, fields):
    with (ROOT/name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

URLS = [
 ('pranjic_px','Pranjic2024','metadata','https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD041341-1'),
 ('pranjic_proxi','Pranjic2024','metadata','https://proteomecentral.proteomexchange.org/api/proxi/v0.1/datasets/PXD041341'),
 ('pranjic_jpost','Pranjic2024','metadata','https://repository.jpostdb.org/entry/JPST002117'),
 ('pranjic_pubmed','Pranjic2024','metadata','https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:38340920&resultType=core&format=json'),
 ('pranjic_doi','Pranjic2024','article_landing','https://doi.org/10.1016/j.ijbiomac.2024.130068'),
 ('evans_pmc','Evans2019','article','https://pmc.ncbi.nlm.nih.gov/articles/PMC6878130/'),
 ('evans_xml','Evans2019','article_and_supplement_links','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6878130/fullTextXML'),
 ('evans_publisher','Evans2019','article_and_supplement_links','https://febs.onlinelibrary.wiley.com/doi/10.1002/1873-3468.13578'),
 ('mcdonald_pmc','McDonald2025','article','https://pmc.ncbi.nlm.nih.gov/articles/PMC12082455/'),
 ('mcdonald_xml','McDonald2025','article_and_supplement_links','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12082455/fullTextXML'),
 ('mcdonald_supp','McDonald2025','supplement_candidate','https://pmc.ncbi.nlm.nih.gov/articles/instance/12082455/bin/gkaf428_supplemental_file.pdf'),
 ('mcdonald_figshare','McDonald2025','code_metadata','https://api.figshare.com/v2/articles/28688966'),
 ('mcdonald_geo','McDonald2025','reused_dataset_metadata','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174145&form=text&view=full'),
 ('mcdonald_geo_counts','McDonald2025','processed_counts_NOT_MS_RAW','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174145/suppl/GSE174145_MistranslatingtRNA_RNAseq_RawCounts.txt.gz'),
]
for stem in ('aggregates','proteomes','pulldown','wtcontrol'):
    URLS.append(('pranjic_'+stem+'_https','Pranjic2024','processed_candidate_https_mirror',f'https://ftp.jpostdb.org/JPST002117/proteinGroups_{stem}.txt'))
    URLS.append(('pranjic_'+stem+'_ftp','Pranjic2024','processed_archived_source_url',f'ftp://ftp.jpostdb.org/JPST002117/proteinGroups_{stem}.txt'))

def fetch(item):
    ident,study,role,url=item
    folder=ROOT/'retrieval'; folder.mkdir(exist_ok=True)
    body=folder/(ident+'.body'); hdr=folder/(ident+'.headers'); err=folder/(ident+'.stderr')
    now=dt.datetime.now(dt.timezone.utc).isoformat()
    # Maximum body 30 MiB: no RAW MS / giant archives / peptide tables.
    cmd=['curl','--silent','--show-error','--location','--connect-timeout','5',
         '--max-time','25','--max-filesize',str(30*1024*1024),
         '--dump-header',str(hdr),'--output',str(body),
         '--write-out','%{http_code}\t%{size_download}\t%{url_effective}\t%{content_type}',url]
    p=subprocess.run(cmd,capture_output=True,text=True)
    err.write_text(p.stderr)
    parts=p.stdout.split('\t'); parts += ['']*(4-len(parts))
    code=parts[0]
    size=body.stat().st_size if body.exists() else 0
    return dict(source_id=ident,study=study,role=role,source_url=url,access_route='curl_live',
      observed_at_utc=now,http_status=code if url.startswith('http') and code!='000' else 'NA_NO_HTTP_RESPONSE',
      transport_code=code,curl_exit=p.returncode,actual_bytes=size,sha256=sha(body) if body.exists() else '',
      local_path=str(body.relative_to(ROOT)) if body.exists() else '',
      final_url=parts[2],content_type=parts[3],status='TRANSFER_COMPLETE_INSPECT_BODY' if p.returncode==0 else 'TRANSPORT_FAILED',
      note=p.stderr.strip().replace('\n',' | '),historical_sha256='',historical_match='')

def definition(c):
    if c in ('Reverse','Potential contaminant','Only identified by site'):
        return 'MaxQuant flag; + marks flagged group; retain in preview, exclude before inference'
    if c.startswith('iBAQ peptides'): return 'Theoretical peptide-count denominator; NOT a replicate intensity'
    if c.startswith('iBAQ'): return 'Intensity-based absolute quantification proxy in instrument units; not cell counts, mass fraction or flux'
    if c.startswith('LFQ intensity'): return 'MaxQuant label-free relative abundance; not synthesis rate or throughput'
    if c.startswith('Ratio H/L normalized'): return 'Normalized heavy/light DnaK pull-down/background enrichment; PD condition key unresolved; not residence time'
    if c.startswith('Ratio'): return 'MaxQuant isotope-channel ratio or ratio statistic; consult publisher processing methods'
    if c in ('Protein IDs','Majority protein IDs'): return 'Semicolon-delimited protein-group accessions; not necessarily one unique protein'
    if c=='Gene names': return 'Semicolon-delimited gene annotations; not a unique join key'
    if c=='Fasta headers': return 'Search-database descriptions; not independent wet-lab strain evidence'
    if c.startswith('Intensity'): return 'MS signal intensity in arbitrary units'
    return 'Verbatim MaxQuant annotation/QC field; unused in Stage I decision; no biological meaning inferred'

def inspect_tables():
    schema=[]; tables=[]; samples=[]
    for path in sorted(RAW.glob('proteinGroups_*.txt')):
        with path.open(newline='') as f:
            reader=csv.DictReader(f,delimiter='\t'); cols=reader.fieldnames
            preview=[]; n=0; malformed=0
            for row in reader:
                n+=1
                if None in row or any(v is None for v in row.values()): malformed+=1
                if n<=3: preview.append(row)
        tsv(path.stem+'_PREVIEW.tsv',preview,cols)
        tables.append(dict(file=path.name,rows=n,columns=len(cols),malformed_rows=malformed,preview_rows=len(preview),bytes=path.stat().st_size,sha256=sha(path)))
        for k,c in enumerate(cols,1):
            schema.append(dict(file=path.name,column_index=k,column=c,definition=definition(c),preview_nonempty=sum(r[c]!='' for r in preview)))
            m=re.fullmatch(r'iBAQ ([ABCD])_(ed_dclpB|wt_dclpB|ed)_(blank|Nva|Val)',c)
            if m:
                rep,genotype,treat=m.groups()
                samples.append(dict(file=path.name,column=c,replicate_label=rep,genotype=genotype,treatment=treat,
                   identity_basis='explicit_header',replication='four labeled samples per arm; biological independence not established by header alone',
                   pairing='same labels across fractions; culture/block pairing not independently authenticated'))
            m=re.fullmatch(r'iBAQ (control|Nva|Val) ([123])',c)
            if m:
                samples.append(dict(file=path.name,column=c,replicate_label=m[2],genotype='WT_according_to_archive',treatment=m[1],identity_basis='header_and_archive',replication='three labeled samples per arm; independent cultures not reverified',pairing='separate WT experiment; do not pair to A-D'))
            if re.fullmatch(r'Ratio H/L normalized PD[1-6]',c):
                samples.append(dict(file=path.name,column=c,replicate_label=c.split()[-1],genotype='ED_DnaK_His6_vs_untagged_ED',treatment='UNRESOLVED_PD_KEY',identity_basis='header; archived methods say three control plus three Val',replication='six mixtures; no Nva arm',pairing='H=tagged_pull_down L=untagged_background; archived methods'))
    tsv('COLUMN_DEFINITIONS.tsv',schema,['file','column_index','column','definition','preview_nonempty'])
    tsv('TABLE_INVENTORY.tsv',tables,['file','rows','columns','malformed_rows','preview_rows','bytes','sha256'])
    tsv('SAMPLE_COLUMNS.tsv',samples,['file','column','replicate_label','genotype','treatment','identity_basis','replication','pairing'])

def manifest():
    rows=[]
    old=json.loads((ARCHIVE/'data/manifest.json').read_text())
    for a in old['artifacts']:
        if a['study_id']=='pranjic_2024' and a.get('local_path'):
            p=ARCHIVE/a['local_path']; digest=sha(p)
            rows.append(dict(source_id=a['artifact_id'],study='Pranjic2024',role=a.get('role',''),source_url=a['source_url'],
              access_route='local_archive_rehashed',observed_at_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
              http_status='NOT_REQUESTED_LOCAL_READ',actual_bytes=p.stat().st_size,sha256=digest,local_path=str(p),
              status='AVAILABLE_VERIFIED_BYTES',historical_sha256=a['sha256'],historical_match=digest==a['sha256'],
              note='Historical acquisition '+a['retrieved_at_utc']+'; historical HTTP '+str(a.get('http_status','not recorded (FTP/cache)'))))
    if (ROOT/'retrieval/attempts.json').exists(): rows+=json.loads((ROOT/'retrieval/attempts.json').read_text())
    if (ROOT/'WEB_ACCESS.tsv').exists():
        with (ROOT/'WEB_ACCESS.tsv').open() as f: rows+=list(csv.DictReader(f,delimiter='\t'))
    fields=['source_id','study','role','source_url','access_route','observed_at_utc','http_status','transport_code','curl_exit','actual_bytes','sha256','local_path','final_url','content_type','status','historical_sha256','historical_match','note']
    tsv('SOURCE_DATA_MANIFEST.tsv',rows,fields)

def source_hashes(verify=False):
    output=ROOT/'ORIGINAL_SOURCE_HASHES.tsv'
    if verify:
        with output.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
        bad=[r['path'] for r in rows if sha(Path(r['path']))!=r['sha256']]
        assert not bad, bad
        return len(rows)
    paths=list(RAW.glob('*'))+[ARCHIVE/x for x in ('README.md','GATE2A_STATUS.md','analysis/GATE3_STATUS.md','data/manifest.json','data/DATA_ACCESS_AUDIT.md')]
    paths+=list((ROOT.parent/'scope_reset_20261011').glob('*'))
    paths+=list((ROOT.parent/'experiment_grounded_rebuild_20261010/novelty_audit/sources').rglob('*'))
    rows=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(set(paths)) if p.is_file()]
    tsv(output.name,rows,['path','bytes','sha256'])
    return len(rows)

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--fetch',action='store_true'); ap.add_argument('--verify',action='store_true'); args=ap.parse_args()
    if args.verify:
        print('Original source files unchanged:',source_hashes(True))
    else:
        if not (ROOT/'ORIGINAL_SOURCE_HASHES.tsv').exists(): source_hashes()
        if args.fetch:
            with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool: results=list(pool.map(fetch,URLS))
            (ROOT/'retrieval/attempts.json').write_text(json.dumps(results,indent=2)+'\n')
        inspect_tables(); manifest()
        print((ROOT/'TABLE_INVENTORY.tsv').read_text())
