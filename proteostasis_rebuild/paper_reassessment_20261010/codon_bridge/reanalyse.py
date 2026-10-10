"""Read-only source reanalysis; writes only beside this script. No flux calibration.
Dependencies: numpy, pandas, scipy, openpyxl. Seed: 20261010.
Exact sign tests concern dataset-level detected-rate direction, not true fidelity.
"""
from pathlib import Path
from itertools import combinations, product
import hashlib, json, sys
import numpy as np
import pandas as pd
from scipy.stats import binomtest

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
LAND = Path('/storage/kiran-stuff/triplet-proof/reviewer_response/position_rates/raw/landerer2024')
STIK = ROOT.parent / 'proteostasis-law-empirical-phase3/data/raw/error_rates/stikeleather_2026'
BASES = 'TCAG'
CODE = dict(zip((''.join(c) for c in product(BASES, repeat=3)),
                'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
PAIRS = [(aa,a,b) for aa in sorted(set(CODE.values())-{'*'})
         for a,b in combinations(sorted(c for c in CODE if CODE[c]==aa),2)]

def holm(p):
    p=np.asarray(p,float); order=np.argsort(p); out=np.empty(len(p))
    out[order]=np.minimum(1,np.maximum.accumulate(p[order]*np.arange(len(p),0,-1)))
    return out

def sign_p(d):
    d=np.asarray(d); d=d[d!=0]
    return binomtest(int((d>0).sum()),len(d),0.5).pvalue if len(d) else 1.0

def contrasts(cc):
    rate=cc.pivot(index='ds',columns='codon',values='rate')
    base=cc.pivot(index='ds',columns='codon',values='base')
    rows=[]; rng=np.random.default_rng(20261010)
    for aa,a,b in PAIRS:
        good=rate[a].notna() & rate[b].notna()
        d=(rate.loc[good,a]-rate.loc[good,b]).to_numpy()
        ratio=base.loc[good,a]/base.loc[good,b]
        # Predefined sensitivity: both depths >=10,000 and within a factor of 2.
        balanced=(base.loc[good,a]>=10000)&(base.loc[good,b]>=10000)&ratio.between(.5,2)
        db=d[balanced.to_numpy()]
        boot=d[rng.integers(0,len(d),size=(4000,len(d)))].mean(axis=1)
        lo,hi=np.quantile(boot,[.025,.975])
        mean=d.mean(); leave=(d.sum()-d)/(len(d)-1) if len(d)>1 else np.array([np.nan])
        rows.append(dict(aa=aa,codon_a=a,codon_b=b,n_pairs=len(d),n_nonzero=int((d!=0).sum()),
            positive=int((d>0).sum()),negative=int((d<0).sum()),mean_delta=mean,
            min_delta=d.min(),max_delta=d.max(),bootstrap_lo=lo,bootstrap_hi=hi,
            sign_p=sign_p(d),lodo_sign_stable=bool(np.all(np.sign(leave)==np.sign(mean))),
            balanced_n=len(db),balanced_delta=db.mean() if len(db) else np.nan,
            balanced_sign_p=sign_p(db)))
    result=pd.DataFrame(rows)
    # One conservative family covers both all-data and balanced-depth tests.
    adj=holm(np.r_[result.sign_p,result.balanced_sign_p])
    result['holm_p']=adj[:len(result)]; result['balanced_holm_p']=adj[len(result):]
    return result

def load_landerer():
    counts=[]; subs=[]; files=[]
    for f in sorted((LAND/'S4/ecoli_filtered').glob('*_codon_counts.csv')):
        d=pd.read_csv(f); d['ds']=f.name.split('_')[0]; counts.append(d); files.append(f)
        s=f.with_name(f.name.replace('_codon_counts','_substitution_errors'))
        t=pd.read_csv(s); t['ds']=d.ds.iloc[0]; subs.append(t); files.append(s)
    c=pd.concat(counts,ignore_index=True); s=pd.concat(subs,ignore_index=True)
    c=c[c.codon.map(CODE).ne('*')].copy()
    assert not c.duplicated(['ds','codon']).any()
    assert np.allclose(c.detection_rate,c.error_count/c.base_count,atol=1e-15)
    observed=s.groupby(['ds','codon']).size()
    reconstructed=c.set_index(['ds','codon']).error_count
    assert (observed.reindex(reconstructed.index,fill_value=0)==reconstructed).all()
    c=c.rename(columns={'base_count':'base','error_count':'err','detection_rate':'rate'})
    c.loc[c.base<=0,'rate']=np.nan
    s2=pd.read_excel(LAND/'Data_S2_error_detection_rate.xlsx',sheet_name='E. coli').set_index('Codon')
    g=c[c.base>0].groupby('codon')
    per=g.rate.agg(['mean','count']).rename(columns={'mean':'zero_inclusive','count':'n_covered'})
    per['conditional']=c[c.err>0].groupby('codon').rate.mean()
    per['pooled']=g.err.sum()/g.base.sum()
    per['n_detected']=c[c.err>0].groupby('codon').size()
    per['err']=g.err.sum(); per['base']=g.base.sum()
    assert (per.n_detected.reindex(s2.index)==np.rint((s2.sd/s2.se)**2)).all()
    residual=float(abs(per.conditional.reindex(s2.index)-s2['mean']).max())
    assert residual<3e-9
    usage=pd.read_csv(ROOT.parent/'envelope-paper/data/raw/global_codon_usage_ecoli.tsv',sep='\t')
    usage=usage[usage.codon.map(CODE).ne('*')].set_index('codon')['count']; usage=usage/usage.sum()
    summary=dict(datasets=int(c.ds.nunique()),cells_possible=80*61,covered=int((c.base>0).sum()),
        zero=int(((c.base>0)&(c.err==0)).sum()),positive=int((c.err>0).sum()),
        absent=int(80*61-(c.base>0).sum()),error_psms=int(c.err.sum()),base=int(c.base.sum()),
        s2_max_error=residual,global_pooled=float(c.err.sum()/c.base.sum()),
        usage_weighted={k:float((usage*per[k]).sum()) for k in ['conditional','zero_inclusive','pooled']})
    return c,s,per,summary,files+[LAND/'Data_S2_error_detection_rate.xlsx']

def load_stikeleather():
    f=STIK/'218_stats.csv'; d=pd.read_csv(f); d['codon']=d.codon.str.replace('U','T')
    assert len(d)==61 and d.codon.nunique()==61
    assert (d.codon.map(CODE)==d.AA).all()
    rows=[]; residual=0
    for j in range(1,4):
        base=d[f'codon_total{j}']; rate=d[f'error_rate{j}']; e=base*rate
        residual=max(residual,float(abs(e-e.round()).max()))
        assert np.all(abs(e-e.round())<1e-7)
        rows.append(pd.DataFrame(dict(ds=str(j),codon=d.codon,base=base,err=e.round().astype(int),rate=rate)))
    cc=pd.concat(rows,ignore_index=True)
    assert np.allclose(d.mean_error_rate,d[[f'error_rate{j}' for j in range(1,4)]].mean(axis=1),atol=1e-14)
    g=cc.groupby('ds')[['err','base']].sum(); rates=g.err/g.base
    return cc,dict(datasets=3,integer_count_residual=residual,replicate_rates=rates.to_dict(),
        mean_rate=float(rates.mean()),se=float(rates.std(ddof=1)/np.sqrt(3)),
        errors_by_rep=g.err.to_dict(),bases_by_rep=g.base.to_dict()),[f]

def flux(s,h,r,a,tau,direct_degradation=0):
    """Partition closure: probabilities h, r, a and d0; no fitted constants.
    h=client entry fraction; r=rescue; a=aggregate/unresolved among clients;
    d0=direct terminal loss outside client entry, disjoint from h.
    """
    if s<0 or tau<0 or not (0<=h<=1 and 0<=direct_degradation<=1-h and 0<=r<=1 and 0<=a<=1-r):
        raise ValueError('invalid flux partition')
    y=1-direct_degradation-h*(1-r)
    return dict(yield_=y,folded=s*y,groel_demand=s*h*tau,
        protease=s*(direct_degradation+h*(1-r-a)),unresolved=s*h*a)

def main():
    c,s,per,lsummary,files=load_landerer(); st,ssummary,sfiles=load_stikeleather()
    c.to_csv(OUT/'landerer_cells.tsv',sep='\t',index=False)
    per.to_csv(OUT/'landerer_codon_rates.tsv',sep='\t')
    st.to_csv(OUT/'stikeleather_cells.tsv',sep='\t',index=False)
    lc=contrasts(c); sc=contrasts(st)
    lc.to_csv(OUT/'landerer_same_aa.tsv',sep='\t',index=False)
    sc.to_csv(OUT/'stikeleather_same_aa.tsv',sep='\t',index=False)
    spectra=pd.crosstab(s.codon,s.destination).reindex(sorted(c.codon.unique()),fill_value=0)
    spectra.to_csv(OUT/'landerer_destination_counts.tsv',sep='\t')
    tv=[]
    for aa,a,b in PAIRS:
        x=spectra.loc[a].to_numpy(); y=spectra.loc[b].to_numpy()
        tv.append(dict(aa=aa,codon_a=a,codon_b=b,err_a=int(x.sum()),err_b=int(y.sum()),
            detected_spectrum_TV=float(abs(x/x.sum()-y/y.sum()).sum()/2)))
    pd.DataFrame(tv).to_csv(OUT/'landerer_spectrum_contrasts.tsv',sep='\t',index=False)
    for df,summary in [(lc,lsummary),(sc,ssummary)]:
        summary['same_aa_pairs']=len(df)
        summary['holm_sign_hits']=int((df.holm_p<.05).sum())
        summary['balanced_holm_sign_hits']=int((df.balanced_holm_p<.05).sum())
        summary['proline_CCA_minus_CCG']=df[(df.codon_a=='CCA')&(df.codon_b=='CCG')].to_dict('records')[0]
    # Condition-dependent expression comparisons are observational, not folding data.
    levels=pd.read_csv(STIK/'codon_error_rates_by_expr_levels.csv')
    mismatch=~np.isclose(levels.error_rate,levels.total_mut_count/levels.total_count,atol=1e-12)
    ssummary['unused_expression_table_inconsistent_rows']=int(mismatch.sum())
    ssummary['unused_expression_table_zero_count_nonzero_rate']=int(((levels.total_mut_count==0)&(levels.error_rate>0)).sum())
    # Save provenance from the existing manifest; verify every input hash independently.
    manifest=json.loads((ROOT.parent/'proteostasis-law-empirical-phase3/data/manifest.json').read_text())
    sources=[]
    for f in files+sfiles+[STIK/'codon_error_rates_by_expr_levels.csv']:
        sources.append(dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    (OUT/'source_hashes.json').write_text(json.dumps(sources,indent=2)+'\n')
    # Extract all original manifest entries mentioning numerical source names, without rewriting them.
    def entries(obj):
        if isinstance(obj,dict):
            if any(isinstance(v,str) and ('218_stats.csv' in v or 'codon_error_rates_by_expr_levels.csv' in v) for v in obj.values()):
                yield obj
            else:
                for v in obj.values(): yield from entries(v)
        elif isinstance(obj,list):
            for v in obj: yield from entries(v)
    (OUT/'stikeleather_manifest_entries.json').write_text(json.dumps(list(entries(manifest)),indent=2)+'\n')
    # Pre-existing report/model arithmetic checks, not a biological simulation.
    from scipy.optimize import brentq
    x=brentq(lambda x:1+4/(1+x)**2-.3*x,0,20)
    mathcheck=dict(raw_threshold=(1-.7**(1/300))/(.7*.3),sub_threshold=(1-.7**(1/300))/.3,
        reduced_fold_x=x,reduced_fold_lambda=x+4*x/(1+x)-.15*x*x)
    (OUT/'numerical_results.json').write_text(json.dumps(dict(landerer=lsummary,stikeleather=ssummary,arithmetic=mathcheck),indent=2)+'\n')
    print(json.dumps(dict(landerer=lsummary,stikeleather=ssummary,arithmetic=mathcheck),indent=2))

if __name__=='__main__': main()
