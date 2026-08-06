"""Sections 3.3, 7 and 9: every number recomputed, then found in the manuscript.

Two claims entered the paper in this pass and neither had a generator behind it
before:

  * Section 3.3's exclusion of Haldane homeostasis, now a sign condition on
    `G_af` rather than the false assertion that all four chain-rule products
    carry the same sign (D067);
  * Sections 7 and 9's account of what happens INSIDE an instability window,
    from the staged integration of C5a (D068).

The direction of assertion matters. Each quantity is recomputed here from
`data/computed/` and the FORMATTED value is then required to appear in the
manuscript text. A test that pinned a literal would keep passing after the data
moved; this one fails, which is the only arrangement that catches a stale
sentence. It is the same discipline that caught a test still pinning a claim
Task B3 had disproved.

The window numbers are asserted against `hopf_window_fate.tsv` and not against
the JSON summary beside it: a serialised summary is not a generator.
"""

from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np
import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parents[2]
_COMPUTED = _REPO_ROOT / "data" / "computed"
_MS = _REPO_ROOT / "manuscript" / "MANUSCRIPT_BMB_v5.md"

_EVAL = ("cycle", "divergent", "fixed_same", "fixed_other")


def _manuscript() -> str:
    return _MS.read_text()


def _grouped():
    """the kinetic-box fold states, labelled by their Section 7 group."""
    G = pd.read_csv(_COMPUTED / "gu_sign.tsv", sep="\t")
    G = G[G["population"] == "kinetic_box"]
    Z = pd.read_csv(_COMPUTED / "hopf_zero_counts.tsv", sep="\t")
    grp = dict(zip(Z["name"], Z["group"]))
    return G.assign(group=G["name"].map(grp).fillna("non_crossing"))


def _window():
    return pd.read_csv(_COMPUTED / "hopf_window_fate.tsv", sep="\t")


class TestSection33IsTheMeasuredSignCondition(unittest.TestCase):

    def setUp(self):
        self.G = _grouped()
        self.text = _manuscript()

    def testTheGroupSizesAreTheOnesQuoted(self):
        n = self.G["group"].value_counts()
        self.assertEqual(int(n["window"]), 47)
        self.assertEqual(int(n["terminal"]), 61)
        self.assertEqual(int(n["non_crossing"]), 2659)
        for c in ("47 window networks", "the 2659 that never cross"):
            self.assertIn(c, self.text)

    def testTheFractionSatisfyingTheSufficientConditionIsWhatIsPrinted(self):
        for g, want in (("window", "93.6%"), ("terminal", "86.9%"),
                        ("non_crossing", "25.3%")):
            s = self.G[self.G["group"] == g]
            frac = float((s["G_af"] >= 0.0).mean())
            self.assertEqual(f"{100 * frac:.1f}%", want,
                             f"{g} moved to {100 * frac:.1f}%")
            self.assertIn(want, self.text)

    def testSection7QuotesTheSameOverlapAsACount(self):
        cross = self.G[self.G["group"].isin(("window", "terminal"))]
        self.assertEqual(int((cross["G_af"] >= 0.0).sum()), 97)
        self.assertEqual(len(cross), 108)
        self.assertIn("97 of these 108 networks", self.text)

    def testTheTermMediansAndTheirSuppressionFactors(self):
        med = {}
        for g in ("window", "non_crossing"):
            s = self.G[self.G["group"] == g]
            base = s["t_growth"].abs().clip(lower=1e-300)
            med[g] = ((s["t_disagg"].abs() / base).median(),
                      (s["t_clear"].abs() / base).median())
        (wd, wc), (nd, nc) = med["window"], med["non_crossing"]
        self.assertEqual(f"{100 * wd:.1f}%", "0.5%")
        self.assertEqual(f"{100 * wc:.1f}%", "6.7%")
        self.assertEqual(f"{100 * nd:.0f}%", "76%")
        self.assertEqual(f"{100 * nc:.0f}%", "71%")
        # the suppression the manuscript names, disaggregation the larger of the
        # two -- which Section 7's own descriptor does not mention
        self.assertEqual(f"{nd / wd:.0f}", "139")
        self.assertEqual(f"{nc / wc:.1f}", "10.6")
        self.assertGreater(nd / wd, nc / wc)
        for c in ("medians of 0.5% and 6.7%", "against 76% and 71%",
                  "suppressed 139-fold and the clearance term 10.6-fold"):
            self.assertIn(c, self.text)

    def testTheResidualRunsTheRightWay(self):
        """G_u is smallest where the proof reaches, not where it does not."""
        G = self.G
        self.assertEqual(int((G["G_u"] <= 0.0).sum()), 0)
        L = pd.read_csv(_COMPUTED / "gu_sign.tsv", sep="\t")
        self.assertEqual(int((L["G_u"] <= 0.0).sum()), 0,
                         "a nonpositive G_u appeared in some population")

        neg = G["G_af"] < 0.0
        self.assertEqual(int(neg.sum()), 1998)
        lo_ok = float(G[~neg]["G_u"].min())
        lo_no = float(G[neg]["G_u"].min())
        self.assertLess(lo_ok, lo_no,
                        "the smallest G_u is no longer where the proof applies")
        self.assertEqual(f"{lo_ok:.1e}", "3.4e-11")
        self.assertEqual(f"{lo_no:.1e}", "7.6e-07")
        self.assertEqual(f"{lo_no / lo_ok:.1e}", "2.3e+04")
        for c in ("3.4×10⁻¹¹", "7.6×10⁻⁷", "2.3×10⁴", "the 1998 states"):
            self.assertIn(c, self.text)

    def testTheChainRuleAgreesWithAnIndependentGradient(self):
        """the four-product form is a check only because the routes differ."""
        G = pd.read_csv(_COMPUTED / "gu_sign.tsv", sep="\t")
        rel = ((G["chain_sum"] - G["G_u"]).abs()
               / G["G_u"].abs().clip(lower=1e-300)).max()
        self.assertLess(float(rel), 1e-8)


class TestTheWindowInteriorWasIntegratedNotInterpolated(unittest.TestCase):

    def setUp(self):
        self.D = _window()
        self.C = self.D[self.D["fate"] == "cycle"]
        self.text = _manuscript()

    def testTheEvaluableAndSettledCountsAreTheOnesPrinted(self):
        ev = self.D[self.D["fate"].isin(_EVAL)]
        self.assertEqual(len(ev), 237)
        self.assertEqual(len(self.C), 227)
        self.assertEqual(int((self.D["fate"] == "divergent").sum()), 0,
                         "a divergent point reappeared; check the horizon")
        for c in ("Of 237 evaluable points 227 carry a settled oscillation",
                  "none diverges"):
            self.assertIn(c, self.text)

    def testEveryInteriorPointCyclesIn34OfThe38Windows(self):
        def verdict(s):
            e = s[s.isin(_EVAL)]
            if e.empty:
                return "none_evaluable"
            if (e == "cycle").all():
                return "all_evaluable_cycle"
            return "no_cycle" if (e == "cycle").sum() == 0 else "mixed"

        v = self.D.groupby("name")["fate"].apply(verdict)
        self.assertEqual(int((v == "all_evaluable_cycle").sum()), 34)
        self.assertEqual(len(v), 38)
        self.assertEqual(int((v == "no_cycle").sum()), 0,
                         "a no-cycle network reappeared")
        self.assertIn("34 of the 38 windows", self.text)

    def testTheOrbitsAreSettledAndThePeriodExceedsTheLinearPrediction(self):
        cv = float(self.C["period_cv"].median())
        self.assertEqual(f"{cv:.1e}", "8.2e-03")
        self.assertIn("median 8.2×10⁻³", self.text)
        r = self.C["period"] / self.C["period_predicted"]
        self.assertEqual(f"{100 * (float(r.median()) - 1.0):.1f}%", "5.0%")
        self.assertIn("by a median of 5.0%", self.text)

    def testAmplitudeVariesSmoothlyRatherThanJumping(self):
        spread = []
        for _, g in self.C.sort_values("frac").groupby("name"):
            a = g["amp_late"].to_numpy(float)
            if a.size >= 2:
                spread.append(a.max() / max(a.min(), 1e-300))
        spread = np.asarray(spread, float)
        self.assertEqual(f"{np.median(spread):.2f}", "1.78")
        self.assertEqual(f"{spread.max():.1f}", "6.0")
        self.assertIn("a median factor of 1.78 and at most 6.0", self.text)

    def testTheUnsettledAndDepartingPointsAreReportedNotAbsorbed(self):
        n = self.D["fate"].value_counts()
        self.assertEqual(int(n.get("not_converged", 0)), 29)
        self.assertEqual(
            self.D[self.D["fate"] == "not_converged"]["name"].nunique(), 11)
        other = self.D[self.D["fate"] == "fixed_other"]
        self.assertEqual(len(other), 10)
        self.assertEqual(other["name"].nunique(), 3)
        # they sit in the opening half of their windows, never the closing half
        self.assertLessEqual(float(other["frac"].max()), 0.5)
        for c in ("29 at 11 networks", "10 at three networks",
                  "draw4627"):
            self.assertIn(c, self.text)

    def testNothingIsAnywhereNearTheDetectionFloor(self):
        """the amplitude-floor account of the non-cycle points cannot hold."""
        self.assertGreater(float(self.C["amp_late"].min()), 1e-3)
        other = self.D[self.D["fate"] == "fixed_other"]
        # not small oscillation -- a different state entirely
        self.assertLess(float(other["a_mean_over_eq"].max()), 1e-2)


class TestSection9MeanBurdenIsBoundedAndIsOneQuantity(unittest.TestCase):

    def setUp(self):
        D = _window()
        self.C = D[D["fate"] == "cycle"]
        self.text = _manuscript()

    def testTheMeanBurdenOverTheCycleIsWhatIsPrinted(self):
        m = self.C["a_mean_over_eq"]
        self.assertEqual(f"{float(m.median()):.2f}", "1.15")
        self.assertEqual(f"{float(m.quantile(0.95)):.2f}", "3.04")
        self.assertEqual(f"{float(m.max()):.2f}", "4.40")
        for c in ("median 1.15 times the equilibrium",
                  "p95 3.04 and at most 4.40"):
            self.assertIn(c, self.text)

    def testAmplitudeAndMeanAreOneQuantitySeenTwice(self):
        rho = float(self.C[["a_amp_rel", "a_mean_over_eq"]]
                    .corr(method="spearman").iloc[0, 1])
        self.assertEqual(f"{rho:.2f}", "0.69")
        self.assertIn("Spearman 0.69", self.text)
        k = max(1, int(round(0.05 * len(self.C))))
        top_mean = set(self.C.nlargest(k, "a_mean_over_eq")["name"])
        top_amp = set(self.C.nlargest(k, "a_amp_rel")["name"])
        self.assertEqual(len(top_mean & top_amp), 4)
        self.assertIn("the upper 5% of each is the same four networks",
                      self.text)

    def testThePeakToPeakSwingIsTheOneQuoted(self):
        self.assertEqual(f"{float(self.C['a_amp_rel'].median()):.1f}", "2.4")
        self.assertIn("median 2.4 times that equilibrium", self.text)


if __name__ == "__main__":
    unittest.main()
