import unittest
import numpy as np
from model import (Network,architecture,constraints,solve,evaluate,radial_capacity,
                   draw,perturbations,pools)

class Tests(unittest.TestCase):
    def test_pool_anchor(self):
        np.testing.assert_allclose(pools(),[14.5782,10.8803,.201287,1.72666,.949851,1.23648,.0116582])

    def test_series_not_or(self):
        n=Network(np.array([[1.],[1.]]),np.ones((1,1)),np.ones(1),np.array(['rescue']),np.ones(1),['serial'])
        x,v,p,e=solve(n,np.array([1.]),np.array([1.,0.]),'architecture',np.ones(1))
        self.assertEqual(v,0)
        self.assertEqual(solve(n,np.array([1.]),np.array([1.,0.]),'scalar',np.ones(1))[1],.5)

    def test_dual_shadow_price(self):
        n=Network(np.ones((1,1)),np.ones((1,1)),np.ones(1),np.array(['rescue']),np.ones(1),['route'])
        _,v,p,e=solve(n,np.array([2.]),np.array([1.]),'architecture',np.ones(1))
        v2=solve(n,np.array([2.]),np.array([1.00001]),'architecture',np.ones(1))[1]
        self.assertAlmostEqual((v2-v)/.00001,p[0],places=8)
        self.assertLess(e,1e-12)

    def test_hall_cut(self):
        # k0 can use m0,m1; k1 only m0. Enumerated min-cuts: 3,3,2,3.
        n=Network(np.array([[1,0,1],[0,1,0.]]),np.array([[1,1,0],[0,0,1.]]),np.ones(3),
                  np.array(['rescue']*3),np.ones(3),['a','b','c'])
        _,v,_,_=solve(n,np.array([1.,2.]),np.array([1.,1.]),'architecture',np.ones(3))
        self.assertEqual(v,2.)

    def test_disposal_never_essential_rescue(self):
        n=architecture(); d=np.array([.1,.1,.1,.1,.1]); C=np.ones(9)*100; C[4]=0
        r=evaluate(n,d,C,'architecture')
        self.assertEqual(r['essential_fraction'],0.)
        self.assertAlmostEqual(r['unserved'],0.)
        self.assertGreater(r['degradation']+r['sequestration'],0.)

    def test_zero_capacity_no_miracles(self):
        d=np.ones(5); r=evaluate(architecture(),d,np.zeros(9),'architecture')
        self.assertEqual(r['folded'],0.); self.assertEqual(r['unserved'],5.)

    def test_finite_storage(self):
        C=np.zeros(9); C[8]=.25
        r=evaluate(architecture(),np.ones(5),C,'architecture')
        self.assertAlmostEqual(r['sequestration'],.25)
        self.assertAlmostEqual(r['unserved'],4.75)

    def test_deterministic_draws(self):
        for a,b in zip(draw(10,3),draw(10,3)): np.testing.assert_array_equal(a,b)
        self.assertFalse(np.array_equal(draw(10,3)[0],draw(10,4)[0]))

    def test_matched_and_nonnegative(self):
        for i in range(5):
            _,C,_,_=draw(10,i); pp=perturbations(C)
            for c in pp.values(): self.assertGreaterEqual(c.min(),0.)
            np.testing.assert_allclose([pp['matched_'+s].sum() for s in ('GroEL','DnaK','ClpB')],
                                       pp['matched_GroEL'].sum(),rtol=0,atol=1e-14)
            for pair in ('TF_DnaK','DnaK_GroEL','DnaK_ClpB','GroEL_ClpB'):
                np.testing.assert_allclose([pp['pair_'+pair+'_'+s].sum() for s in ('ab','aa','bb')],
                                           pp['pair_'+pair+'_ab'].sum(),rtol=0,atol=1e-14)

    def test_nesting_and_baseline(self):
        for i in range(5):
            s,C,costs,_=draw(11,i); n=architecture(costs); d=.9*radial_capacity(n,s,C)*s
            r=evaluate(n,d,C,'architecture')
            self.assertAlmostEqual(r['folded_fraction'],1.)
            self.assertAlmostEqual(r['essential_fraction'],1.)
            vals=[]
            for kind in ('architecture','independent','scalar'):
                x,v,_,_=solve(n,d,C,kind,n.productive); vals.append(v)
                if kind=='architecture':
                    for relaxed in ('independent','scalar'):
                        A,b=constraints(n,C,relaxed); self.assertLessEqual(np.max(A@x-b),1e-8)
            self.assertTrue(np.all(np.diff(vals)>-1e-8))

    def test_scalar_target_invariance(self):
        s,C,costs,_=draw(13,2); n=architecture(costs); d=s
        vals=[evaluate(n,d,c,'scalar')['folded'] for name,c in perturbations(C).items() if name.startswith('matched')]
        np.testing.assert_allclose(vals,vals[0],atol=1e-10)

    def test_label_permutation_negative_control(self):
        s,C,costs,_=draw(17,2); n=architecture(costs)
        perm=np.array([4,3,2,1,0,6,5,8,7]); other=Network(n.A[perm],n.H,n.productive,n.sink,n.atp,n.names)
        a=evaluate(n,s,C,'architecture'); b=evaluate(other,s,C[perm],'architecture')
        for k in a:
            self.assertAlmostEqual(a[k],b[k],places=7)

    def test_ablation_has_finite_bypass(self):
        d=np.ones(5)*.01; C=np.ones(9); C[4]=0
        r=evaluate(architecture(control='remove_obligate'),d,C,'architecture')
        self.assertAlmostEqual(r['essential_fraction'],1.)
        C[1]=0
        self.assertEqual(evaluate(architecture(control='remove_obligate'),d,C,'architecture')['essential_fraction'],0.)

    def test_nonidentification_witness(self):
        H=np.eye(2); C=np.array([1.,0.]); d=np.array([.5,.5]); reward=np.array([1.,0.])
        def net(A):return Network(np.array(A),H,np.ones(2),np.array(['rescue']*2),np.ones(2),['e','o'])
        a,b=net([[1.,1.],[0.,0.]]),net([[0.,1.],[1.,0.]])
        self.assertEqual(solve(a,d,C,'architecture',reward)[1],.5)
        self.assertEqual(solve(b,d,C,'architecture',reward)[1],0.)
        self.assertEqual(solve(a,d,C,'scalar',reward)[1],solve(b,d,C,'scalar',reward)[1])

if __name__=='__main__':unittest.main()
