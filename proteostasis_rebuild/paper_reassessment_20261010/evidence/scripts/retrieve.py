"""Download public sources; keep failures distinct from usable experimental files."""
import csv, hashlib, json, subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
TARGETS = [
 ('babu_article.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3930520/fullTextXML','article_XML','verified_article_identifier'),
 ('babu_S2.xls','https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1004120.s023&type=supplementary','genetic_interaction_table','publisher_supplement_DOI_verified'),
 ('babu_S1.xls','https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1004120.s022&type=supplementary','query_strains','publisher_supplement_DOI_verified'),
 ('babu_S10.xls','https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1004120.s031&type=supplementary','GI_overlap_not_client_overlap','publisher_supplement_DOI_verified'),
 ('babu_S16.xls','https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1004120.s037&type=supplementary','strains','publisher_supplement_DOI_verified'),
 ('babu_protocol_S2.pdf','https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1004120.s007&type=supplementary','scoring_protocol','publisher_supplement_DOI_verified'),
 ('babu_protocol_S16.pdf','https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1004120.s021&type=supplementary','screen_protocol','publisher_supplement_DOI_verified'),
 ('calloni_main.pdf','https://pure.mpg.de/rest/items/item_1602419_3/component/file_1836221/content','article_PDF','verified_web_readable'),
 ('kerner_supp_mirror.pdf','https://citeseerx.ist.psu.edu/document?doi=4087a93abb0c94d87af0e09c2c904f50c46c3d48&repid=rep1&type=pdf','supplement_PDF','search_index_identifies_original_supplement_mirror_not_authenticated'),
 ('kerner_mmc1_candidate.pdf','https://ars.els-cdn.com/content/image/1-s2.0-S0092867405005441-mmc1.pdf','supplement_candidate','filename_pattern_candidate_not_verified'),
 ('calloni_mmc1_candidate.pdf','https://ars.els-cdn.com/content/image/1-s2.0-S2211124711000179-mmc1.pdf','supplement_candidate','filename_pattern_candidate_not_verified'),
 ('keio_article.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC1681482/fullTextXML','deletion_collection','verified_article_identifier'),
 ('fujiwara_article.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3212837/fullTextXML','dependence_followup','verified_article_identifier'),
 ('calloni_tranche.html','https://proteomecommons.org/','legacy_repository','URL_in_original_article_no_modern_accession_verified'),
]

def fetch(t):
    name,url,kind,identity=t
    path=ROOT/'sources'/name
    if path.exists():
        raise RuntimeError(f'Refusing to replace source {path}')
    r=subprocess.run(['curl','--location','--fail','--silent','--show-error','--connect-timeout','8','--max-time','20','--output',str(path),'--write-out','%{http_code}\t%{url_effective}\t%{content_type}',url],capture_output=True,text=True)
    data=path.read_bytes() if path.exists() else b''
    signature_ok = (data.startswith(b'%PDF') if name.endswith('.pdf') else data.startswith(bytes.fromhex('d0cf11e0a1b11ae1')) or data.startswith(b'PK') if name.endswith('.xls') else b'<article' in data if name.endswith('.xml') else bool(data))
    valid=r.returncode==0 and signature_ok
    row=dict(filename=name,url=url,kind=kind,url_identity=identity,attempted_utc=datetime.now(timezone.utc).isoformat(),curl_exit=r.returncode,response=r.stdout,error=r.stderr.strip(),status='downloaded_signature_checked' if valid else 'failed_or_invalid',bytes=len(data),sha256=hashlib.sha256(data).hexdigest() if data else '',usable_experimental_table=int(valid and kind=='genetic_interaction_table'))
    if data and not valid:
        path.rename(path.with_suffix(path.suffix+'.invalid'))
    elif path.exists() and not data:
        path.unlink()
    return row

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=6) as pool:
        rows=list(pool.map(fetch,TARGETS))
    (ROOT/'retrieval_attempts.json').write_text(json.dumps(rows,indent=2)+'\n')
    with (ROOT/'retrieval_attempts.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    print(json.dumps({'attempts':len(rows),'successful_downloads':sum(r['status']=='downloaded_signature_checked' for r in rows),'errors':sorted(set(r['error'] for r in rows))},indent=2))
