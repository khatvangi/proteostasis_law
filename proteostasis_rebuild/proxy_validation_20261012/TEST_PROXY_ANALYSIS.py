#!/usr/bin/env python3
"""Synthetic invariants plus independent raw-archive reconstruction."""
import csv, math, unittest, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
import PROXY_ANALYSIS as a

class ProxyTests(unittest.TestCase):
    def test_filtering(self):
        d=pd.DataFrame({'Protein IDs':['good','REV__x','CON__y','p;q','dup','dup','site','count',''], 'Number of proteins':['1']*7+['2','1'],'Reverse':['']*9,'Potential contaminant':['']*9,'Only identified by site':['']*6+['+','','']})
        keep,_=a.eligible(d)
        self.assertEqual(list(d.loc[keep,'Protein IDs']),['good'])
        d.loc[0,'Reverse']='+';self.assertFalse(a.eligible(d)[0].any())
    def test_contrast_sign(self):
        y=pd.DataFrame({'ed_blank':[1.],'ed_Nva':[3.],'ed_dclpB_blank':[2.],'ed_dclpB_Nva':[7.]})
        self.assertEqual(a.interaction(y,'Nva').iloc[0],3.)
    def test_scale_invariance_and_zero_rejection(self):
        rng=np.random.default_rng(17);cs=a.cols(['blank','Nva'])
        d=pd.DataFrame(rng.uniform(.1,30,(8,len(cs))),columns=cs)
        s=a.normalized(d,list(d.index),['blank','Nva'])
        scaled=d*np.geomspace(.01,100,len(cs))
        np.testing.assert_allclose(s,a.normalized(scaled,list(d.index),['blank','Nva']))
        d.iloc[0,0]=0
        with self.assertRaises(AssertionError): a.normalized(d,list(d.index),['blank','Nva'])
        d.iloc[0,0]=np.nan
        with self.assertRaises(AssertionError): a.normalized(d,list(d.index),['blank','Nva'])
    def test_independent_fraction_averaging(self):
        s={f:pd.DataFrame({f'iBAQ {r}_{g}_blank':[v] for g in a.GENOS for r,v in zip('ABCD',vs)}) for f,vs in [('aggregates',[.1,.2,.3,.8]),('proteomes',[.4,.1,.2,.1])]}
        y=a.ys(s,['blank'])
        self.assertAlmostEqual(y.ed_blank.iloc[0],math.log2(.35/.2))
        # Permuting labels in only one fraction cannot change the full mean ratio.
        p=s['proteomes'].copy();p.columns=[c.replace('A_','X_').replace('D_','A_').replace('X_','D_') for c in p]
        np.testing.assert_allclose(y,a.ys({'aggregates':s['aggregates'],'proteomes':p},['blank']))
        self.assertNotAlmostEqual(y.ed_blank.iloc[0],float(np.mean(np.log2(np.array([.1,.2,.3,.8])/np.array([.4,.1,.2,.1])))))
    def test_protocol_and_input_hashes(self):
        self.assertEqual(a.sha(a.OUT/'PROTOCOL.md'),json.loads((a.OUT/'PROTOCOL_FREEZE.json').read_text())['sha256'])
        for row in csv.DictReader((a.OUT/'SOURCE_HASHES.tsv').open(),delimiter='\t'):
            self.assertEqual(a.sha(Path(row['path'])),row['sha256'])
    def test_independent_raw_reconstruction(self):
        # csv and scalar math independently check selection, sums and ALL primary protein estimates.
        tables={}
        for f in a.FRACTIONS:
            with (a.RAW/f'proteinGroups_{f}.txt').open() as handle: rows=list(csv.DictReader(handle,delimiter='\t'))
            counts={}
            for r in rows: counts[r['Protein IDs'].strip()]=counts.get(r['Protein IDs'].strip(),0)+1
            clean={}
            for r in rows:
                i=r['Protein IDs'].strip()
                if not i or ';' in i or 'REV__' in i or 'CON__' in i or counts[i]!=1: continue
                if any(r[c].strip() for c in ['Reverse','Potential contaminant','Only identified by site']):continue
                if float(r['Number of proteins'])!=1:continue
                try: v={c:float(r[c]) for c in a.cols(['blank','Nva'])}
                except ValueError:continue
                if all(math.isfinite(x) and x>0 for x in v.values()):clean[i]=v
            tables[f]=clean
        ids=sorted(set(tables['aggregates'])&set(tables['proteomes']))
        result=pd.read_csv(a.OUT/'primary_PROTEIN_RESULTS.tsv',sep='\t').set_index('protein_id')
        self.assertEqual(ids,list(result.index))
        denom={f:{c:math.fsum(tables[f][i][c] for i in ids) for c in a.cols(['blank','Nva'])} for f in a.FRACTIONS}
        for i in ids:
            y={}
            for g in a.GENOS:
                for t in ['blank','Nva']:
                    means={f:math.fsum(tables[f][i][f'iBAQ {r}_{g}_{t}']/denom[f][f'iBAQ {r}_{g}_{t}'] for r in 'ABCD')/4 for f in a.FRACTIONS}
                    y[g+'_'+t]=math.log2(means['aggregates']/means['proteomes'])
                    self.assertAlmostEqual(y[g+'_'+t],result.loc[i,g+'_'+t],places=11)
            value=y['ed_dclpB_Nva']-y['ed_dclpB_blank']-y['ed_Nva']+y['ed_blank']
            self.assertAlmostEqual(value,result.loc[i,'interaction_Nva'],places=11)
    def test_summary_and_specificity(self):
        sums=pd.read_csv(a.OUT/'DISTRIBUTION_SUMMARY.tsv',sep='\t').set_index('analysis')
        for name,ts in [('primary',['Nva']),('secondary',['Val']),('common',['Nva','Val'])]:
            d=pd.read_csv(a.OUT/f'{name}_PROTEIN_RESULTS.tsv',sep='\t')
            for t in ts:
                v=d['interaction_'+t];row=sums.loc[name+'_'+t]
                self.assertEqual(len(v),row.n)
                self.assertAlmostEqual(v.median(),row['median'],places=12)
                self.assertAlmostEqual((v>0).mean(),row.positive_fraction,places=12)
                np.testing.assert_allclose(v.quantile([.25,.75]),row[['q25','q75']].astype(float))
            if name=='common':np.testing.assert_allclose(d.Nva_minus_Val,d.interaction_Nva-d.interaction_Val,atol=1e-12)
    def test_controls_and_robustness(self):
        s={f:pd.read_csv(a.OUT/f'primary_{f}_shares.tsv',sep='\t').set_index('protein_id') for f in a.FRACTIONS}
        rr=pd.read_csv(a.OUT/'ROBUSTNESS_PROTEINS.tsv',sep='\t').set_index('protein_id')
        for scope in ['aggregates','proteomes','both']:
            for label in 'ABCD':
                rest='ABCD'.replace(label,'')
                la=rest if scope in ['both','aggregates'] else 'ABCD';lp=rest if scope in ['both','proteomes'] else 'ABCD'
                np.testing.assert_allclose(rr[f'LOO_{scope}_{label}'],a.interaction(a.ys(s,['blank','Nva'],la,lp),'Nva'),atol=1e-10)
        nc=pd.read_csv(a.OUT/'NEGATIVE_CONTROL_PROTEINS.tsv',sep='\t')
        blank=nc[nc.treatment=='blank'].set_index('protein_id').split_interaction
        nva=nc[nc.treatment=='Nva'].set_index('protein_id').split_interaction
        np.testing.assert_allclose(nva-blank,rr.primary_AB_minus_CD,atol=1e-10)

if __name__=='__main__': unittest.main(verbosity=2)
