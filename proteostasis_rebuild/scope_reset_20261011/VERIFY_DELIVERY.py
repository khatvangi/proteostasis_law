#!/usr/bin/env python3
"""Read-only integrity and consistency audit of this delivery (not biological tests)."""
from pathlib import Path
import csv
import hashlib
import re
D=Path(__file__).resolve().parent
checks=[]
def check(name,condition):
    checks.append((name,bool(condition)))
    print(('PASS' if condition else 'FAIL')+' '+name)
required=['COMPLETE_INTEGRATED_REPORT.md','REVISED_MANUSCRIPT.md','MATHEMATICAL_CHECKS.py',
          'CHECK_RESULTS.txt','CLAIM_AND_EVIDENCE_MATRIX.tsv','ADVERSARIAL_AUDIT.md','README.md']
check('all requested deliverables exist',all((D/x).is_file() and (D/x).stat().st_size for x in required))
original=D.parent.parent/'envelope-paper/manuscript/MANUSCRIPT.md'
check('frozen original byte-identical',original.read_bytes()==(D/'ORIGINAL_FROZEN_MANUSCRIPT.md').read_bytes())
m=(D/'REVISED_MANUSCRIPT.md').read_text()
check('canonical authors preserved',original.read_text().splitlines()[2] in m)
words=len(m.split())
check('manuscript 2500-4500 words including references',2500<=words<=4500)
print('Manuscript word count:',words)
with (D/'SOURCE_MANIFEST.tsv').open() as f: sources=list(csv.DictReader(f,delimiter='\t'))
check('all source hashes unchanged',all(hashlib.sha256(Path(r['source_file']).read_bytes()).hexdigest()==r['sha256'] for r in sources))
with (D/'CLAIM_AND_EVIDENCE_MATRIX.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
valid=True
for r in rows:
    ls=Path(r['source_file']).read_text().splitlines(); a=int(r['line_start']); b=int(r['line_end'])
    valid &= 1<=a<=b<=len(ls) and ls[a-1][:240]==r['anchor_excerpt']
check('all matrix line ranges and anchor excerpts valid',valid)
check('matrix uses allowed statuses',all(r['status'] in {'observed','derived','conditional','unsupported','withdrawn'} for r in rows))
check('all manuscript claim IDs covered',set(re.findall(r'\[C(\d\d)\]',m))<=set(r['claim_id'][1:] for r in rows))
check('22 distinct claim IDs and 34 rows',len(rows)==34 and len(set(r['claim_id'] for r in rows))==22)
body,refs=m.split('## References',1)
refids=set(re.findall(r'\*\*([A-Z][A-Za-z0-9]+)\*\*',refs))
cited=set()
for group in re.findall(r'\[([^\]]+)\]',body):
    for token in group.split(', '):
        if re.fullmatch(r'[A-Z][A-Za-z0-9]+',token) and not re.fullmatch(r'C\d\d',token): cited.add(token)
check('literature IDs resolve and all listed references cited',cited==refids)
links_ok=True
for p in D.glob('*.md'):
    # Frozen original retains historical relative figure paths; do not rewrite it.
    if p.name=='ORIGINAL_FROZEN_MANUSCRIPT.md': continue
    for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
        if target.startswith(('http:','https:','#')): continue
        links_ok &= (p.parent/target.split('#')[0]).exists()
check('all active local Markdown links resolve',links_ok)
readme=(D/'README.md').read_text()
check('start-here map covers every other delivery file',all(p.name in readme for p in D.iterdir() if p.is_file() and p.name!='README.md'))
counts=re.search(r'COUNTS: (\d+) run; (\d+) passed; (\d+) failed; (\d+) errors; (\d+) skipped', (D/'CHECK_RESULTS.txt').read_text())
vals=tuple(map(int, counts.groups())) if counts else ()
check('mathematical run has no failures; optional symbolic availability explicit', vals in {(20,20,0,0,0),(20,18,0,0,2)})
report=(D/'COMPLETE_INTEGRATED_REPORT.md').read_text()
check('report closes scientific scope',report.rstrip().endswith('**STOP / FINISHED SCIENTIFIC SCOPE: bounded reconstruction complete; empirical validation outstanding, not silently claimed complete.**'))
check('current report and manuscript state not submission-ready',all('not submission-ready' in t for t in [report,m]))
print(f'PACKAGE COUNTS: {len(checks)} checks; {sum(v for _,v in checks)} passed; {sum(not v for _,v in checks)} failed')
print('SOURCE COUNTS:',len(sources),'files;',len(rows),'evidence rows')
print('Biological validation: not performed. Manual scientific audit: ADVERSARIAL_AUDIT.md.')
raise SystemExit(0 if all(v for _,v in checks) else 1)
