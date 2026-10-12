import urllib.request, urllib.error, hashlib, csv, datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parent
URLS=[('martens_pmc.html','https://pmc.ncbi.nlm.nih.gov/articles/PMC12262517/'),('martens_europepmc.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12262517/fullTextXML'),('geo_record.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104303&targ=self&form=text&view=full'),('geo_supp_listing.html','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104303/suppl/'),('martens_supp.pdf','https://pmc.ncbi.nlm.nih.gov/articles/instance/12262517/bin/NIHPP2025.06.25.661590V1-supplement-1.pdf')]
URLS += [(n, 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104303/suppl/'+n) for n in ['GSE104303_mRNA_mapped_read_counts.tsv.gz','GSE104303_mapped_read_counts.tsv.gz']]
LIMIT=5*1024*1024
with (ROOT/'HTTP_RETRIEVAL.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['utc','url','http_status','bytes','sha256','saved_as','note'])
 for name,url in URLS:
  status='NO_HTTP_RESPONSE';data=b'';note=''
  try:
   with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'ProxyEvidenceAudit/1.0'}),timeout=18) as r:
    status=r.status
    if int(r.headers.get('Content-Length','0'))>LIMIT: raise ValueError('size cap exceeded')
    data=r.read(LIMIT+1)
    if len(data)>LIMIT: raise ValueError('size cap exceeded')
  except urllib.error.HTTPError as e: status=e.code;data=e.read(LIMIT);note=str(e)
  except Exception as e: note=repr(e)
  saved=''
  if data: (ROOT/name).write_bytes(data);saved=name
  w.writerow([datetime.datetime.now(datetime.timezone.utc).isoformat(),url,status,len(data),hashlib.sha256(data).hexdigest() if data else 'NA_NO_BYTES',saved,note]);f.flush()
  print(name,status,len(data),note)
