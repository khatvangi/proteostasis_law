import unittest
import numpy as np
import pandas as pd
import reanalyse as r

class Checks(unittest.TestCase):
    def test_code_and_family(self):
        self.assertEqual(len(r.PAIRS),87)
        self.assertEqual(r.CODE['CCG'],'P')
        self.assertEqual(r.CODE['TGA'],'*')
        self.assertTrue(all(r.CODE[a]==r.CODE[b]==aa for aa,a,b in r.PAIRS))
    def test_holm_known_and_symmetry(self):
        np.testing.assert_allclose(r.holm([.01,.03,.2]),[.03,.06,.2])
        self.assertEqual(r.sign_p([1,1,1]),.25)
        self.assertEqual(r.sign_p([0,0]),1)
        self.assertEqual(r.sign_p([1,2,-1]),r.sign_p([-1,-2,1]))
    def test_data_reconstruction(self):
        c,s,p,summary,_=r.load_landerer()
        self.assertEqual(summary['zero'],3015)
        self.assertEqual(summary['absent'],19)
        self.assertEqual(summary['error_psms'],58990)
        self.assertTrue((p.conditional>=p.zero_inclusive).all())
        self.assertTrue((c.base>0).all())
    def test_stikeleather_counts_and_summary(self):
        c,s,_=r.load_stikeleather()
        self.assertLess(s['integer_count_residual'],1e-7)
        self.assertAlmostEqual(s['mean_rate'],.00181810928199195,places=14)
        self.assertTrue((c.err<=c.base).all())
    def test_mass_accounting_and_signs(self):
        # Algebraic property tests, synthetic numbers solely as software fixtures.
        rng=np.random.default_rng(12)
        for _ in range(1000):
            h,rp,ap,tau=rng.random(4); ap*=1-rp; d0=(1-h)*rng.random()
            a=r.flux(2,h,rp,ap,tau,d0)
            self.assertAlmostEqual(a['folded']+a['protease']+a['unresolved'],2)
            b=r.flux(2,h/2,rp,ap,tau,d0)
            self.assertGreaterEqual(b['yield_'],a['yield_'])
            self.assertLessEqual(b['groel_demand'],a['groel_demand'])
            self.assertLessEqual(b['protease'],a['protease'])
        with self.assertRaises(ValueError): r.flux(1,.8,.8,.4,1)
    def test_no_zero_coverage_conversion(self):
        c,_,_=r.load_stikeleather()
        p=r.contrasts(c)
        self.assertGreaterEqual(p.sign_p.min(),.25)
        self.assertEqual((p.holm_p<.05).sum(),0)

if __name__=='__main__': unittest.main()
