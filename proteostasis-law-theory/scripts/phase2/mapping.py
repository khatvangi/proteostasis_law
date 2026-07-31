"""boron <-> nitrogen parameter mapping along the sequestration ladder epsilon.

this module is the single source of truth for the matched-implementation
protocol.  the derivation, and the distinction between the EXACT identity at
epsilon = 0 and the O(epsilon) APPROXIMATION at finite epsilon, is written up in
`theory/MATCHED_IMPLEMENTATION_PROTOCOL.md`.  the short version:

boron (`scripts/proteostasis/model.py`) carries TOTAL burdens (u, a) and solves
an implicit rapid-equilibrium closure

    uf = u / su,   su = 1 + cf/kappa_cu + df/kappa_du
    af = a / sa,   sa = 1 + cf/kappa_ca + df/kappa_da
    cf = c_tot / (1 + nu + uf/kappa_cu + af/kappa_ca)
    df = d_tot / (1 +      uf/kappa_du + af/kappa_da)

nitrogen (`.../nitrogen-check/src/proteostasis_model.py`) carries FREE burdens,
i.e. it asserts su = sa = 1 and keeps only the resource-side occupancy.

the ladder ties every binding constant to one scalar

    kappa_cu = kappa_ca = c_tot / epsilon        c_u = c_a = epsilon / c_tot
    kappa_du = kappa_da = d_tot / epsilon        d_u = d_a = epsilon / d_tot

so that epsilon = max(c_tot/kappa_cu, d_tot/kappa_du, c_tot/kappa_ca,
d_tot/kappa_da) exactly, and, crucially, BOTH sides use the same resource-side
occupancy at every epsilon.  the two models then differ ONLY by the substrate
sequestration factors su, sa.  under this ladder kappa_ca = kappa_cu and
kappa_da = kappa_du, hence

    su = sa = 1 + epsilon * (cf/c_tot + df/d_tot)

is a single scalar, which is what makes the discrepancy cleanly first order.

gauge.  boron's phase-1 nondimensionalization fixes the concentration scale as
phi = c_tot + d_tot, so c_tot + d_tot == 1 there.  the matched benchmark instead
adopts the NITROGEN gauge: time factor mu = 1 and d_tot ≡ `D_TOT_GAUGE` = 1, so
that ref_capacity maps onto c_tot directly and rho_U, rho_A map onto the two
degradation capacities.  this is a rescaling of the concentration unit, not a
change of dynamics, but it does mean the boron parameter vectors produced here
are NOT in boron's phase-1 phi-normalized gauge and must not be compared to
phase-1 parameter tables coordinate by coordinate.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, Tuple

# --- the epsilon ladder ----------------------------------------------------
# 1e-6 is the identity anchor; 1.0 and 2.0 are the strongly sequestered end.
# boron's own phase-1 baseline (kappa_cu = 0.3, c_tot = 0.6) sits at
# epsilon = 2.0 on the chaperone axis and epsilon = 0.8 on the degradation axis,
# so the ladder brackets it rather than reproducing it exactly.
EPSILON_LADDER: Tuple[float, ...] = (1e-6, 1e-3, 1e-2, 1e-1, 0.3, 1.0, 2.0)

#: epsilon values used by the T0 right-hand-side identity test.
T0_EPSILONS: Tuple[float, ...] = (1e-6, 1e-3, 1e-2, 1e-1)

#: gauge constants (see module docstring).
MU_GAUGE = 1.0          # time-scale factor; 1.0 means boron and nitrogen clocks agree
D_TOT_GAUGE = 1.0       # boron d_tot; not identifiable from the nitrogen groups

#: nitrogen `Parameters` defaults that `lhs_sweep.py` never touches.  these are
#: the pinned coordinates of the matched benchmark.
PINNED_NITROGEN: Dict[str, float] = {
    "j_fold": 0.0,
    "j_other": 0.0,
    "o_load": 0.0,
    "c_tot": 1.0,               # nitrogen normalises free resources to fractions
    "d_tot": 1.0,
    "ref_k": 0.2,               # K_ref / U0
    "deg_u_k": 0.4,
    "deg_a_k": 0.5,
    "disaggregation": 0.35,
    "disaggregation_k": 0.3,
    "m": 2.0,
    "threshold_u": 1.0,
    "threshold_a": 1.0,
}

#: the seven coordinates `lhs_sweep.parameters_from_row` samples, in order,
#: with (transform, lo, hi).  transcribed from
#: nitrogen:.../src/lhs_sweep.py:41-48.
LHS_COORDINATES = (
    ("j_u",            "log",    0.001, 1.0),
    ("n_load",         "linear", 0.0,   4.0),
    ("nucleation",     "log",    0.005, 0.5),
    ("growth",         "log",    0.01,  1.0),
    ("ref_capacity",   "log",    0.25,  4.0),
    ("deg_u_capacity", "log",    0.1,   3.0),
    ("deg_a_capacity", "log",    0.1,   3.0),
)


@dataclass(frozen=True)
class NitrogenParams:
    """dimensionless parameters in nitrogen's convention (free substrate).

    field names and defaults mirror
    nitrogen:.../src/proteostasis_model.py::Parameters exactly, so a value of
    this dataclass can be handed to the real nitrogen module unchanged.
    """

    j_u: float = 0.02
    j_fold: float = 0.0
    j_other: float = 0.0
    n_load: float = 0.5
    o_load: float = 0.0
    c_u: float = 1.0
    c_a: float = 1.0
    d_u: float = 1.0
    d_a: float = 1.0
    c_tot: float = 1.0
    d_tot: float = 1.0
    ref_capacity: float = 0.6
    ref_k: float = 0.2
    deg_u_capacity: float = 1.0
    deg_a_capacity: float = 0.7
    deg_u_k: float = 0.4
    deg_a_k: float = 0.5
    nucleation: float = 0.08
    m: float = 2.0
    growth: float = 0.25
    disaggregation: float = 0.35
    disaggregation_k: float = 0.3
    threshold_u: float = 1.0
    threshold_a: float = 1.0

    @property
    def j(self) -> float:
        return self.j_u + self.j_fold + self.j_other

    def asDict(self) -> Dict[str, float]:
        return asdict(self)

    def removalCeiling(self) -> float:
        """supremum of total removal flux; the free-limit analogue of
        `model.removalCeiling`.  cf <= c_tot, df <= d_tot and every saturating
        factor is < 1."""
        return (self.ref_capacity * self.c_tot
                + (self.deg_u_capacity + self.deg_a_capacity) * self.d_tot)

    def residualScale(self) -> float:
        """matches `equilibria.residualScale` so the two protocols apply the
        same relative residual tolerance to both model forms."""
        return max(self.j, 1e-30, 1e-6 * self.removalCeiling())


# ---------------------------------------------------------------------------
# the mapping itself
# ---------------------------------------------------------------------------


def bindingConstants(epsilon: float, c_tot: float, d_tot: float) -> Dict[str, float]:
    """the four boron kappas and the four nitrogen reciprocals at one epsilon.

    epsilon = 0 is rejected: boron's `Params.validate` requires strictly
    positive, finite kappas, so the exact identity is reached as a limit in the
    code even though it is an algebraic identity on paper.  use
    `freeLimitParams` for the epsilon-free model form instead.
    """
    if not epsilon > 0.0:
        raise ValueError("epsilon must be strictly positive; the epsilon = 0 "
                         "model form is `nitrogen_limit`, not a boron Params")
    return {
        "kappa_cu": c_tot / epsilon, "kappa_ca": c_tot / epsilon,
        "kappa_du": d_tot / epsilon, "kappa_da": d_tot / epsilon,
        "c_u": epsilon / c_tot, "c_a": epsilon / c_tot,
        "d_u": epsilon / d_tot, "d_a": epsilon / d_tot,
    }


def nitrogenToBoron(q: NitrogenParams, epsilon: float) -> Dict[str, float]:
    """nitrogen parameter vector -> boron `Params` kwargs at sequestration epsilon.

    mu = 1 gauge, d_tot pinned to `D_TOT_GAUGE`.  every line below is an
    algebraic identity between the two flux expressions, verified term by term
    in `theory/MATCHED_IMPLEMENTATION_PROTOCOL.md` section 3.
    """
    c_tot = q.ref_capacity * q.c_tot / MU_GAUGE     # ref_capacity = mu * c_tot
    d_tot = D_TOT_GAUGE
    kb = bindingConstants(epsilon, c_tot, d_tot)
    return {
        "j": q.j / MU_GAUGE,
        "nu": q.n_load + q.o_load,
        "c_tot": c_tot,
        "d_tot": d_tot,
        "rho_U": q.deg_u_capacity * q.d_tot / (MU_GAUGE * d_tot),
        "rho_A": q.deg_a_capacity * q.d_tot / (MU_GAUGE * d_tot),
        "alpha_n": q.nucleation / MU_GAUGE,
        "alpha_g": q.growth / MU_GAUGE,
        "alpha_d": q.disaggregation * q.c_tot / (MU_GAUGE * c_tot),
        "m": q.m,
        "kappa_ref": q.ref_k,
        "kappa_u": q.deg_u_k,
        "kappa_a": q.deg_a_k,
        "kappa_dis": q.disaggregation_k,
        "kappa_cu": kb["kappa_cu"], "kappa_ca": kb["kappa_ca"],
        "kappa_du": kb["kappa_du"], "kappa_da": kb["kappa_da"],
    }


def boronToNitrogen(p, epsilon: float) -> NitrogenParams:
    """boron `Params` -> nitrogen parameter vector at sequestration epsilon.

    the inverse of `nitrogenToBoron` in the mu = 1 gauge.  `p` may be a
    `proteostasis.model.Params` or anything with the same attribute names, so
    this module never has to import the boron package.
    """
    kb = bindingConstants(epsilon, p.c_tot, p.d_tot)
    return NitrogenParams(
        j_u=MU_GAUGE * p.j, j_fold=0.0, j_other=0.0,
        n_load=p.nu, o_load=0.0,
        c_u=kb["c_u"], c_a=kb["c_a"], d_u=kb["d_u"], d_a=kb["d_a"],
        c_tot=1.0, d_tot=1.0,
        ref_capacity=MU_GAUGE * p.c_tot,
        ref_k=p.kappa_ref,
        deg_u_capacity=MU_GAUGE * p.rho_U * p.d_tot,
        deg_a_capacity=MU_GAUGE * p.rho_A * p.d_tot,
        deg_u_k=p.kappa_u, deg_a_k=p.kappa_a,
        nucleation=MU_GAUGE * p.alpha_n, m=p.m,
        growth=MU_GAUGE * p.alpha_g,
        disaggregation=MU_GAUGE * p.alpha_d * p.c_tot,
        disaggregation_k=p.kappa_dis,
    )


def sequestrationEpsilon(p) -> float:
    """epsilon read back off an arbitrary boron `Params`.

    the audit's definition: epsilon = max(c_tot/kappa_cu, d_tot/kappa_du,
    c_tot/kappa_ca, d_tot/kappa_da).  on the ladder all four agree; off it this
    returns the worst axis, which is the quantity that bounds su - 1.
    """
    return max(p.c_tot / p.kappa_cu, p.d_tot / p.kappa_du,
               p.c_tot / p.kappa_ca, p.d_tot / p.kappa_da)


def boundSubstrateFraction(q: NitrogenParams, epsilon: float) -> float:
    """1 - 1/su at the low-burden state, i.e. the fraction of the damaged
    monomer pool that boron holds machinery-bound and nitrogen sets to zero.

    evaluated at u = a = 0, where cf = c_tot/(1+nu) and df = d_tot, so
    su - 1 = epsilon * (1/(1+nu) + 1).  this is the physical size of the
    difference between the two model forms, and it is what the ladder scans.
    """
    sigma0 = 1.0 / (1.0 + q.n_load + q.o_load) + 1.0
    su = 1.0 + epsilon * sigma0
    return 1.0 - 1.0 / su
