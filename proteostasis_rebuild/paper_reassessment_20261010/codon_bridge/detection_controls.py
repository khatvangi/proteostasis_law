"""Diagnostic nulls; synthetic counts are NOT biological evidence or fitted fluxes.
Equal rate for every synonym within each dataset and amino acid, actual exposure.
Independent Poisson count sampling tests whether unequal exposure alone defeats
naive sign inference on observed rates. PSM independence is a diagnostic assumption.
"""
import json
import numpy as np
import pandas as pd
from scipy.stats import binom
import reanalyse as r

def p_sign(k,n):
    return np.minimum(1,2*binom.cdf(np.minimum(k,n-k),n,.5))

def calibration(cc,simulations=500):
    base=cc.pivot(index='ds',columns='codon',values='base').fillna(0)
    err=cc.pivot(index='ds',columns='codon',values='err').fillna(0)
    cols=list(base.columns); n=base.to_numpy(); e=err.to_numpy(); rates=np.zeros_like(n,dtype=float)
    for aa in sorted(set(r.CODE.values())-{'*'}):
        ix=[j for j,c in enumerate(cols) if r.CODE[c]==aa]
        den=n[:,ix].sum(axis=1); num=e[:,ix].sum(axis=1)
        rates[:,ix]=np.divide(num,den,out=np.zeros_like(num,dtype=float),where=den>0)[:,None]
    rng=np.random.default_rng(20261010)
    x=rng.poisson((rates*n)[None,:,:],size=(simulations,*n.shape))
    obs=np.divide(x,n[None,:,:],out=np.zeros_like(x,dtype=float),where=n[None,:,:]>0)
    ps=[]
    for aa,a,b in r.PAIRS:
        i,j=cols.index(a),cols.index(b); mask=(n[:,i]>0)&(n[:,j]>0)
        d=obs[:,mask,i]-obs[:,mask,j]
        ps.append(p_sign((d>0).sum(axis=1),(d!=0).sum(axis=1)))
        bal=mask&(n[:,i]>=10000)&(n[:,j]>=10000)&(n[:,i]>=.5*n[:,j])&(n[:,i]<=2*n[:,j])
        db=obs[:,bal,i]-obs[:,bal,j]
        ps.append(p_sign((db>0).sum(axis=1),(db!=0).sum(axis=1)))
    p=np.array(ps).T; hit=p<.05/174
    # split the family: even columns are all-data tests, odd columns balanced-depth tests
    per=[dict(aa=aa,codon_a=a,codon_b=b,null_hit_all=float(hit[:,2*k].mean()),
        null_hit_balanced=float(hit[:,2*k+1].mean())) for k,(aa,a,b) in enumerate(r.PAIRS)]
    return dict(simulations=simulations,any_false_sign_hit=float(hit.any(axis=1).mean()),
        any_false_sign_hit_all_data=float(hit[:,0::2].any(axis=1).mean()),
        any_false_sign_hit_balanced=float(hit[:,1::2].any(axis=1).mean()),
        mean_false_hits_all_data=float(hit[:,0::2].sum(axis=1).mean()),
        per_pair=per,
        null='Equal detected-event rates within dataset/amino-acid; observed unequal exposures; independent Poisson sampling',
        inference='Diagnostic only: sign-test equality of signs is not equality of underlying detection rates')

def main():
    c,s,_,_,_=r.load_landerer(); st,_,_=r.load_stikeleather()
    checks={name:calibration(d) for name,d in [('landerer',c),('stikeleather',st)]}
    # Filtering retained numerator rows cannot recover undetected peptides or a decoy FDR.
    quality=s.PeptideProphet_Probability
    checks['landerer_psm_quality']={'min':float(quality.min()),'below_0_99':int((quality<.99).sum()),
        'duplicate_ds_spectrum_codon':int(s.duplicated(['ds','Spectrum','codon']).sum()),
        'warning':'PeptideProphet score is not a substitution-localization/decoy error guarantee.'}
    (r.OUT/'detection_controls.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(checks,indent=2))

if __name__=='__main__': main()
