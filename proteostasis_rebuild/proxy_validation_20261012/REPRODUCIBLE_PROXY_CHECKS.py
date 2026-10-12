#!/usr/bin/env python3
"""Standard-library algebra/provenance checks. Arbitrary rational witnesses, NOT biology.
Run from any directory. Writes CHECK_RESULTS.txt only inside this directory.
No network, random draws, biological simulation, threshold estimation, or statistical fit.
"""
from pathlib import Path
from fractions import Fraction as Q
import csv, hashlib, re, sys
R=Path(__file__).resolve().parent
results=[]
def check(name, fn):
 try:
  assert fn(), 'condition false'
  results.append(('PASS',name))
 except Exception as e: results.append(('FAIL',name+': '+repr(e)))
def rows(name):
 with (R/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))

# Conservative U,C,A,N state transfers: capture, release, productive folding,
# aggregation and disaggregation are each a single equal-and-opposite transfer.
def mass_balance():
 X=[Q(2),Q(3),Q(5),Q(7)];mu=Q(1,10);d=[Q(1,7),Q(1,8),Q(1,9),Q(1,11)];s=Q(13)
 dx=[s-mu*X[0]-d[0]*X[0]]+[-mu*X[i]-d[i]*X[i] for i in range(1,4)]
 for a,b,j in [(0,1,Q(2,3)),(1,0,Q(1,4)),(1,3,Q(2,7)),(0,2,Q(3,5)),(2,0,Q(1,9)),(3,0,Q(1,13))]:dx[a]-=j;dx[b]+=j
 return sum(dx)==s-mu*sum(X)-sum(a*b for a,b in zip(d,X))
check('internal client transfers cancel exactly',mass_balance)
def steady_mass():
 X=[Q(2),Q(3),Q(5),Q(7)];mu=Q(1,10);loss=sum(Q(1,i+2)*x for i,x in enumerate(X));s=mu*sum(X)+loss
 return s-mu*sum(X)-loss==0 and sum(X)<=s/mu
check('steady-state mass balance implies s >= mu P',steady_mass)
check('zero dilution and zero losses forbid positive-source steady state',lambda: all(Q(2)-Q(0)*p-Q(0)!=0 for p in [Q(0),Q(1),Q(100)]))
# Exact dimensions as exponents (chain, site, visit, volume, time).
chain=(1,0,0,-1,0);site=(0,1,0,-1,0);flux=(1,0,0,-1,-1)
mu=(0,0,0,0,-1);b=(0,0,1,-1,-1);nu=(0,1,-1,0,0);tau=(0,0,0,0,1)
def mul(*xs):return tuple(map(sum,zip(*xs)))
check('source and dilution units agree',lambda:mul(mu,chain)==flux)
check('visit rate times sites/visit times residence equals site concentration',lambda:mul(b,nu,tau)==site)
check('reject adding aggregate concentration to processing flux',lambda:chain!=flux)
check('reject equating dimensionless induction with occupancy',lambda:(0,0,0,0,0)!=site)
check('binding rate units k F U equal chain flux',lambda:mul((0,-1,0,1,-1),site,chain)==flux)
check('competition activity a=gamma*k*U/q dimensionless',lambda:mul((-1,1,0,0,0),(0,-1,0,1,-1),chain,(0,0,0,0,1))==(0,0,0,0,0))
# Arbitrary exact stationary example: b in visits, q in inverse time, nu sites/visit.
bb=[Q(2),Q(3)];qq=[Q(2),Q(3)];nn=[Q(1),Q(2)];cc=[x/y for x,y in zip(bb,qq)];M=Q(5)
check('admitted visit residence identity',lambda:sum(n*c for n,c in zip(nn,cc))==sum(x*n/q for x,n,q in zip(bb,nn,qq)))
check('finite pool leaves nonnegative free sites',lambda:M-sum(n*c for n,c in zip(nn,cc))==2)
check('bound state equations stationary',lambda:all(x-q*c==0 for x,q,c in zip(bb,qq,cc)))
# Fixed pool competition exact difference and symbolic derivative witness.
a=[Q(2),Q(3)];MM=Q(12)
def occupancy(aa):return [MM*x/(1+sum(aa)) for x in aa]
check('shared-pool conservation at equilibrium',lambda:sum(occupancy(a))+MM/(1+sum(a))==MM)
check('competitor increase displaces other fixed client',lambda:occupancy([a[0],a[1]+1])[0]<occupancy(a)[0])
check('strict derivative sign requires occupied focal client',lambda:-MM*a[0]/(1+sum(a))**2<0 and -MM*0/(1+sum(a))**2==0)
def separate_pool(focal,competitor):return MM*focal/(1+focal)
check('separate pools give zero cross displacement',lambda:separate_pool(a[0],a[1])==separate_pool(a[0],a[1]+1))
check('feedback can reverse fixed-pool sign inference',lambda:2*MM*a[0]/(1+sum(a)+1)>occupancy(a)[0])

# Source hashes protect cited read-only originals, not a claim of project-wide snapshot.
ledger=rows('SOURCE_LEDGER.tsv')
local=[x for x in ledger if x['source_path_or_url'].startswith('/')]
check('cited local originals and registration retain SHA256',lambda:all(hashlib.sha256(Path(x['source_path_or_url']).read_bytes()).hexdigest()==x['sha256'] for x in local))
nums=rows('PUBLISHED_NUMBERS.tsv');rates=[x for x in nums if x['quantity']=='spontaneous_refolding_k1']
# Verify transcription against source strings; not an additional biological test.
source=(R/'reingewertz_readable.txt').read_text()
def transcription():
 for x in rates:
  v=x['condition'];v='Wild‐type' if v=='Wild-type' else v
  pattern=re.escape(v)+r'\s+([\d.]+)\s+±\s+([\d.]+)'
  m=re.search(pattern,source)
  if not m or Q(m[1])!=Q(x['value']) or Q(m[2])!=Q(x['uncertainty']):return False
 return len(rates)==12
check('all 12 Table 1 k1 and SE cells transcribed exactly',transcription)
check('matrix includes all required claim classes',lambda:{'proven/math','observed','extrapolated','unavailable'}<=set(x['claim_class'] for x in rows('PROXY_EVIDENCE_MATRIX.tsv')))
report=(R/'ONE_INTEGRATED_REPORT.md').read_text()
check('integrated report within requested 2000-3000 words',lambda:2000<=len(report.split())<=3000)
check('old VOID and missing arm stated in report',lambda:'remains VOID' in report and 'DnaK–norvaline is absent' in report)
check('preprint status explicit',lambda:'preprint, not peer reviewed' in report)
retr=rows('HTTP_RETRIEVAL.tsv')
check('small processed GEO files both attempted',lambda:all(any(n in x['url'] for x in retr) for n in ['GSE104303_mRNA_mapped_read_counts.tsv.gz','GSE104303_mapped_read_counts.tsv.gz']))
check('no checksum invented for absent HTTP bytes',lambda:all(x['sha256']=='NA_NO_BYTES' for x in retr if x['bytes']=='0'))
required=['README.md','ONE_INTEGRATED_REPORT.md','PROXY_EVIDENCE_MATRIX.tsv','SOURCE_LEDGER.tsv','ADVERSARIAL_PROXY_AUDIT.md','PREREGISTERED_PROXY_RELATION.md','PUBLISHED_NUMBERS.tsv']
check('required dossier files exist',lambda:all((R/f).is_file() for f in required))
failed=sum(s=='FAIL' for s,n in results)
text='PROXY-ONLY CHECK RESULTS\nStandard Python; exact algebraic witnesses are not biological simulations.\n\n'+'\n'.join(f'{s}\t{n}' for s,n in results)
text+=f'\n\nAutomated checks: {len(results)-failed} PASS / {failed} FAIL / 0 ERROR\nReport word count: {len(report.split())}\n'
text+='Biological Prediction A: NOT EXECUTABLE (matched published numerical protein cells not retrieved/authenticated).\nPrediction B: published descriptive Table 1 reproduction target; no new inferential reanalysis.\nOld Gate 3: VOID (unchanged).\nNo biological threshold/headroom test passed or was estimated.\n'
(R/'CHECK_RESULTS.txt').write_text(text);print(text);sys.exit(bool(failed))
