#!/usr/bin/env python3
"""Bounded descriptive analysis; see frozen PROTOCOL.md. No inferential protein n."""
from pathlib import Path
import hashlib, json, datetime, sys
import numpy as np
import pandas as pd
OUT = Path(__file__).resolve().parent
RAW = OUT.parent.parent / 'proteostasis-law-empirical-phase3/data/raw/outcomes/pranjic_2024'
META = OUT.parent / 'experimental_decision_20261011'
LABELS = 'ABCD'
GENOS = ('ed','ed_dclpB')
FRACTIONS = ('aggregates','proteomes')
GATES = []
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(df, name): df.to_csv(OUT/name, sep='\t',index=False, float_format='%.15g')
def gate(name, ok, detail, fatal=True):
    GATES.append(dict(gate=name,status='PASS' if ok else ('BLOCKED' if fatal else 'LIMITATION'),detail=detail))
    save(pd.DataFrame(GATES),'GATES.tsv')
    if not ok and fatal: raise RuntimeError(f'{name}: {detail}')
def cols(treatments):
    return [f'iBAQ {r}_{g}_{t}' for g in GENOS for t in treatments for r in LABELS]
def eligible(d):
    ids=d['Protein IDs'].fillna('').str.strip()
    reasons=pd.DataFrame(index=d.index)
    reasons['reverse']=d['Reverse'].fillna('').str.strip().ne('') | ids.str.contains('REV__',regex=False)
    reasons['contaminant']=d['Potential contaminant'].fillna('').str.strip().ne('') | ids.str.contains('CON__',regex=False)
    reasons['site_only']=d['Only identified by site'].fillna('').str.strip().ne('')
    reasons['ambiguous']=ids.eq('') | ids.str.contains(';',regex=False) | pd.to_numeric(d['Number of proteins'],errors='coerce').ne(1)
    reasons['duplicate_accession']=ids.duplicated(keep=False)
    return ~reasons.any(axis=1), reasons

def normalized(d, ids, treatments):
    v=d.loc[ids,cols(treatments)].astype(float)
    assert np.isfinite(v.to_numpy()).all() and (v.to_numpy()>0).all()
    s=v/v.sum(axis=0)
    assert np.allclose(s.sum(axis=0),1,rtol=0,atol=1e-12)
    return s

def ys(s, treatments, la=LABELS, lp=LABELS):
    y={}
    for g in GENOS:
        for t in treatments:
            a=s['aggregates'][[f'iBAQ {r}_{g}_{t}' for r in la]].mean(axis=1)
            p=s['proteomes'][[f'iBAQ {r}_{g}_{t}' for r in lp]].mean(axis=1)
            y[g+'_'+t]=np.log2(a/p)
    return pd.DataFrame(y)
def interaction(y,t):
    return (y['ed_dclpB_'+t]-y['ed_dclpB_blank'])-(y['ed_'+t]-y['ed_blank'])
def summary(name,v,reference=None):
    v=np.asarray(v,dtype=float)
    assert len(v) and np.isfinite(v).all()
    q=np.quantile(v,[.25,.5,.75],method='linear')
    z=dict(analysis=name,n=len(v),q25=q[0],median=q[1],q75=q[2],positive_fraction=float(np.mean(v>0)),minimum=v.min(),maximum=v.max())
    z['sign_agreement_with_primary']=float(np.mean(np.sign(v)==np.sign(reference))) if reference is not None else np.nan
    return z

def main():
    freeze=json.loads((OUT/'PROTOCOL_FREEZE.json').read_text())
    gate('protocol_unchanged',sha(OUT/'PROTOCOL.md')==freeze['sha256'],'Protocol hash matches freeze record')
    receipt=OUT/'OUTCOME_READ_RECEIPT.json'
    if not receipt.exists(): receipt.write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol_sha256':freeze['sha256'],'stage':'before full raw table reads; header inspection occurred after freeze'},indent=2)+'\n')
    manifest=pd.read_csv(META/'SOURCE_DATA_MANIFEST.tsv',sep='\t')
    md=pd.read_csv(META/'SAMPLE_COLUMNS.tsv',sep='\t')
    hashes=[]; data={}; qc=[]; exclusions=[]; sampleqc=[]
    for f in FRACTIONS:
        path=RAW/f'proteinGroups_{f}.txt'
        expected=manifest.loc[manifest.source_id.eq('pranjic_proteinGroups_'+f),'sha256'].iloc[0]
        h=sha(path); hashes.append(dict(path=str(path),bytes=path.stat().st_size,sha256=h,expected_sha256=expected))
        gate(f+'_archive_hash',h==expected,'Actual bytes match historical metadata SHA256')
        d=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False)
        required=['Protein IDs','Number of proteins','Reverse','Potential contaminant','Only identified by site']+cols(['blank','Nva','Val'])
        gate(f+'_schema',set(required)<=set(d),'Required identity, QC and 24 iBAQ columns present')
        sub=md[md.file.eq(path.name)&md.genotype.isin(GENOS)&md.treatment.isin(['blank','Nva','Val'])]
        mapping_ok=len(sub)==24 and set(sub.column)==set(cols(['blank','Nva','Val']))
        mapping_ok=mapping_ok and all(row.column==f'iBAQ {row.replicate_label}_{row.genotype}_{row.treatment}' for row in sub.itertuples())
        gate(f+'_metadata',mapping_ok,'A-D labels and condition identities agree with SAMPLE_COLUMNS.tsv')
        keep,reasons=eligible(d)
        qc.append(dict(fraction=f,stage='raw_rows',n=len(d)))
        for reason in reasons:
            qc.append(dict(fraction=f,stage='excluded_'+reason+'_overlapping',n=int(reasons[reason].sum())))
        qc.append(dict(fraction=f,stage='excluded_union',n=int((~keep).sum())))
        qc.append(dict(fraction=f,stage='eligible_unique',n=int(keep.sum())))
        for i in d.index[~keep]: exclusions.append(dict(fraction=f,protein_id=d.loc[i,'Protein IDs'],reasons=';'.join(reasons.columns[reasons.loc[i]])))
        e=d.loc[keep].copy();e.index=e['Protein IDs'].str.strip()
        for c in cols(['blank','Nva','Val']):
            e[c]=pd.to_numeric(e[c],errors='coerce')
            v=e[c].to_numpy(float)
            sampleqc.append(dict(fraction=f,column=c,eligible_n=len(e),missing_or_nonfinite=int((~np.isfinite(v)).sum()),nonpositive_finite=int((np.isfinite(v)&(v<=0)).sum()),positive_n=int((np.isfinite(v)&(v>0)).sum()),eligible_positive_total=float(v[np.isfinite(v)&(v>0)].sum())))
        data[f]=e
    for path in [META/'SAMPLE_COLUMNS.tsv',META/'SOURCE_DATA_MANIFEST.tsv']:
        hashes.append(dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path),expected_sha256='not independently pinned'))
    save(pd.DataFrame(hashes),'SOURCE_HASHES.tsv');save(pd.DataFrame(exclusions),'EXCLUSIONS.tsv')
    save(pd.DataFrame(qc),'FILTER_QC.tsv');save(pd.DataFrame(sampleqc),'SAMPLE_QC.tsv')
    shared=sorted(set(data['aggregates'].index)&set(data['proteomes'].index))
    coverage=[]; totals=[]; summaries=[]; universes={}; normalized_sets={}
    for name,ts in [('primary',['blank','Nva']),('secondary',['blank','Val']),('common',['blank','Nva','Val'])]:
        masks={}
        for f in FRACTIONS:
            v=data[f][cols(ts)].to_numpy(float)
            masks[f]=set(data[f].index[(np.isfinite(v)&(v>0)).all(axis=1)])
        ids=sorted(masks['aggregates']&masks['proteomes'])
        gate(name+'_nonempty',len(ids)>0,f'{len(ids)} complete-positive matched proteins')
        universes[name]=ids
        coverage.append(dict(universe=name,n=len(ids),eligible_shared=len(shared),fraction_shared=len(ids)/len(shared),aggregate_complete=len(masks['aggregates']),proteome_complete=len(masks['proteomes']),required_measurements_per_protein=len(cols(ts))*2))
        s={f:normalized(data[f],ids,ts) for f in FRACTIONS};normalized_sets[name]=s
        for f in FRACTIONS:
            for c in cols(ts):
                v=data[f][c].to_numpy(float);den=float(v[np.isfinite(v)&(v>0)].sum());ret=float(data[f].loc[ids,c].sum())
                totals.append(dict(universe=name,fraction=f,column=c,matched_total=ret,eligible_positive_total=den,retained_signal_fraction=ret/den,share_sum=float(s[f][c].sum())))
            x=s[f].copy();x.insert(0,'protein_id',ids);save(x,f'{name}_{f}_shares.tsv')
        y=ys(s,ts); result=y.copy()
        for t in ts[1:]:
            result['interaction_'+t]=interaction(y,t)
            summaries.append(summary(name+'_'+t,result['interaction_'+t]))
        if name=='common':
            result['Nva_minus_Val']=result.interaction_Nva-result.interaction_Val
            summaries.append(summary('common_Nva_minus_Val',result.Nva_minus_Val))
        result.insert(0,'protein_id',ids);save(result,f'{name}_PROTEIN_RESULTS.tsv')
    s=normalized_sets['primary'];base=interaction(ys(s,['blank','Nva']),'Nva')
    robustness=pd.DataFrame(index=universes['primary'])
    for scope in ['aggregates','proteomes','both']:
        for label in LABELS:
            rest=LABELS.replace(label,'');la=rest if scope in ['aggregates','both'] else LABELS;lp=rest if scope in ['proteomes','both'] else LABELS
            key=f'LOO_{scope}_{label}';v=interaction(ys(s,['blank','Nva'],la,lp),'Nva');robustness[key]=v
            summaries.append(summary(key,v,base))
    for labelset in ['AB','CD']:
        v=interaction(ys(s,['blank','Nva'],labelset,labelset),'Nva');robustness['primary_'+labelset]=v
        summaries.append(summary('primary_'+labelset,v,base))
    v=robustness.primary_AB-robustness.primary_CD;robustness['primary_AB_minus_CD']=v;summaries.append(summary('primary_AB_minus_CD',v))
    robustness.insert(0,'protein_id',robustness.index);save(robustness,'ROBUSTNESS_PROTEINS.tsv')
    controls=[]
    for t in ['blank','Nva','Val']:
        name='secondary' if t=='Val' else 'primary';s=normalized_sets[name]
        ab=ys(s,[t],'AB','AB');cd=ys(s,[t],'CD','CD')
        v=(ab['ed_dclpB_'+t]-cd['ed_dclpB_'+t])-(ab['ed_'+t]-cd['ed_'+t])
        summaries.append(summary('negative_control_'+t,v))
        controls.extend(dict(treatment=t,universe=name,protein_id=i,split_interaction=float(val)) for i,val in v.items())
    save(pd.DataFrame(controls),'NEGATIVE_CONTROL_PROTEINS.tsv');save(pd.DataFrame(coverage),'COVERAGE.tsv');save(pd.DataFrame(totals),'DENOMINATOR_QC.tsv');save(pd.DataFrame(summaries),'DISTRIBUTION_SUMMARY.tsv')
    gate('normalization',True,'All sample shares sum to one within 1e-12; finite strictly positive input and denominators')
    gate('pairing_and_biological_n',False,'Culture independence and pairing unauthenticated; descriptive fraction means only',False)
    gate('absolute_burden_capacity_headroom',False,'Not identified by compositional iBAQ; no aggregate mass or site-time calibration',False)
    gate('preregistration_confirmation',False,'Exploratory existing archives; historical gate remains VOID',False)
    (OUT/'RUN_ENVIRONMENT.json').write_text(json.dumps({'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2)+'\n')
    print(pd.DataFrame(summaries).to_string(index=False))
if __name__=='__main__':
    try: main()
    except Exception as exc:
        (OUT/'ANALYSIS_BLOCKER.txt').write_text(f'{type(exc).__name__}: {exc}\n')
        raise
