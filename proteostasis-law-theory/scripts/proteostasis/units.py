"""dimensional algebra and the dimensional (unit-carrying) shadow model.

phase 1 analysis runs on a nondimensional model, but a nondimensionalization is
only trustworthy if the dimensional model it came from is consistent. this
module supplies:

  1. `Dim` -- an exponent-tracking dimension type. addition of mismatched
     dimensions raises `DimensionError`, so a wrong rate-constant dimension is a
     hard failure rather than a silent rescaling.
  2. `Quantity` -- value + dimension, used to evaluate the rate laws
     symbolically-in-units.
  3. `DimensionalParams` and `dimensionalRhs` -- the model written with explicit
     units (amount volume^-1, time).
  4. `toNondimensional` -- the exact map onto `model.Params`, whose correctness
     is checked by a trajectory-rescaling test rather than by inspection.

base dimensions: amount (e.g. molecules or moles), volume, time.
concentration is amount * volume^-1; a flux is concentration * time^-1.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict


class DimensionError(ValueError):
    """raised when quantities with incompatible dimensions are combined."""


_BASE = ("amount", "volume", "time")


@dataclass(frozen=True)
class Dim:
    """a physical dimension held as integer/rational exponents of base units."""

    amount: float = 0.0
    volume: float = 0.0
    time: float = 0.0

    def asDict(self) -> Dict[str, float]:
        return {b: getattr(self, b) for b in _BASE}

    def __mul__(self, other: "Dim") -> "Dim":
        return Dim(*(getattr(self, b) + getattr(other, b) for b in _BASE))

    def __truediv__(self, other: "Dim") -> "Dim":
        return Dim(*(getattr(self, b) - getattr(other, b) for b in _BASE))

    def __pow__(self, k: float) -> "Dim":
        return Dim(*(getattr(self, b) * k for b in _BASE))

    def isDimensionless(self) -> bool:
        return all(abs(getattr(self, b)) < 1e-12 for b in _BASE)

    def matches(self, other: "Dim") -> bool:
        return all(abs(getattr(self, b) - getattr(other, b)) < 1e-12 for b in _BASE)

    def __str__(self) -> str:
        parts = [f"{b}^{getattr(self, b):g}" for b in _BASE if abs(getattr(self, b)) > 1e-12]
        return " ".join(parts) if parts else "dimensionless"


DIMENSIONLESS = Dim()
AMOUNT = Dim(amount=1.0)
VOLUME = Dim(volume=1.0)
TIME = Dim(time=1.0)
CONC = AMOUNT / VOLUME
FLUX = CONC / TIME
RATE = DIMENSIONLESS / TIME                 # first-order rate constant, time^-1
SECOND_ORDER = RATE / CONC                  # second-order rate constant


@dataclass(frozen=True)
class Quantity:
    """a numeric value tagged with a dimension.

    arithmetic mirrors physical rules: multiplication/division combine
    exponents, addition and subtraction demand an exact dimension match.
    """

    value: float
    dim: Dim = field(default=DIMENSIONLESS)

    def __add__(self, other: "Quantity") -> "Quantity":
        if not self.dim.matches(other.dim):
            raise DimensionError(f"cannot add [{self.dim}] to [{other.dim}]")
        return Quantity(self.value + other.value, self.dim)

    def __sub__(self, other: "Quantity") -> "Quantity":
        if not self.dim.matches(other.dim):
            raise DimensionError(f"cannot subtract [{other.dim}] from [{self.dim}]")
        return Quantity(self.value - other.value, self.dim)

    def __mul__(self, other) -> "Quantity":
        if isinstance(other, (int, float)):
            return Quantity(self.value * other, self.dim)
        return Quantity(self.value * other.value, self.dim * other.dim)

    __rmul__ = __mul__

    def __truediv__(self, other) -> "Quantity":
        if isinstance(other, (int, float)):
            return Quantity(self.value / other, self.dim)
        return Quantity(self.value / other.value, self.dim / other.dim)

    def __pow__(self, k: float) -> "Quantity":
        return Quantity(self.value ** k, self.dim ** k)

    def __neg__(self) -> "Quantity":
        return Quantity(-self.value, self.dim)

    def require(self, dim: Dim, what: str = "quantity") -> "Quantity":
        """assert a dimension; used to pin the expected units of each term."""
        if not self.dim.matches(dim):
            raise DimensionError(f"{what} has [{self.dim}], expected [{dim}]")
        return self


def one() -> Quantity:
    """the dimensionless unit, needed for `1 + x/K` style denominators."""
    return Quantity(1.0, DIMENSIONLESS)


# ---------------------------------------------------------------------------
# dimensional shadow model
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DimensionalParams:
    """the two-state conserved-resource model with explicit units.

    concentrations are amount volume^-1; `A` is carried in monomer-equivalent
    concentration so that transfer between `U` and `A` conserves mass 1:1.
    """

    # loads
    J: float = 1.0e-3          # total abnormal-monomer influx, conc time^-1
    N_free: float = 5.0e-1     # free ordinary nascent-chain concentration, conc
    K_N: float = 1.0           # nascent-chain / chaperone dissociation constant, conc

    # conserved rescue pools, conc
    C_tot: float = 6.0e-1
    D_tot: float = 4.0e-1

    # kinetics
    k_ref: float = 1.0e-2      # per-chaperone refolding turnover, time^-1
    k_U: float = 1.0e-2        # per-protease soluble degradation turnover, time^-1
    k_A: float = 5.0e-3        # per-protease aggregate clearance turnover, time^-1
    k_dis: float = 3.0e-3      # per-chaperone disaggregation turnover, time^-1
    k_n: float = 5.0e-3        # nucleation coefficient, conc^(1-m) time^-1
    k_g: float = 1.0e-2        # aggregate growth coefficient, conc^-1 time^-1
    m: float = 2.0             # nucleation order, dimensionless, m > 1

    # michaelis / dissociation constants, conc
    K_ref: float = 5.0e-1
    K_U: float = 5.0e-1
    K_A: float = 5.0e-1
    K_dis: float = 5.0e-1
    K_CU: float = 3.0e-1       # chaperone-U dissociation constant
    K_CA: float = 3.0e-1       # chaperone-aggregate dissociation constant
    K_DU: float = 5.0e-1       # protease-U dissociation constant
    K_DA: float = 5.0e-1       # protease-aggregate dissociation constant


def dimensionalRateLaws(U_f: Quantity, A_f: Quantity, C_f: Quantity,
                        D_f: Quantity, dp: DimensionalParams) -> Dict[str, Quantity]:
    """evaluate every rate law as a `Quantity`, pinning each to flux units.

    dimensions of the rate constants are declared here, not assumed. if any
    declaration is wrong the `.require(FLUX, ...)` call raises.
    """
    k_ref = Quantity(dp.k_ref, RATE)
    k_U = Quantity(dp.k_U, RATE)
    k_A = Quantity(dp.k_A, RATE)
    k_dis = Quantity(dp.k_dis, RATE)
    # nucleation coefficient must carry conc^(1-m) time^-1 for k_n * U^m to be a flux
    k_n = Quantity(dp.k_n, (CONC ** (1.0 - dp.m)) / TIME)
    k_g = Quantity(dp.k_g, SECOND_ORDER)

    K_ref = Quantity(dp.K_ref, CONC)
    K_U = Quantity(dp.K_U, CONC)
    K_A = Quantity(dp.K_A, CONC)
    K_dis = Quantity(dp.K_dis, CONC)

    laws = {
        "J": Quantity(dp.J, FLUX),
        "v_ref": k_ref * C_f * (U_f / (K_ref + U_f)),
        "v_degU": k_U * D_f * (U_f / (K_U + U_f)),
        "v_degA": k_A * D_f * (A_f / (K_A + A_f)),
        "v_dis": k_dis * C_f * (A_f / (K_dis + A_f)),
        "nucleation": k_n * (U_f ** dp.m),
        "growth": k_g * U_f * A_f,
    }
    for name, q in laws.items():
        q.require(FLUX, f"rate law '{name}'")
    return laws


def dimensionalRhs(U_f: Quantity, A_f: Quantity, C_f: Quantity,
                   D_f: Quantity, dp: DimensionalParams):
    """assemble dU/dt and dA/dt from unit-carrying terms.

    every `+`/`-` below is a dimensional check: a term with the wrong units
    raises `DimensionError` instead of quietly producing a number.
    """
    r = dimensionalRateLaws(U_f, A_f, C_f, D_f, dp)
    dU = r["J"] - r["v_ref"] - r["v_degU"] - r["nucleation"] - r["growth"] + r["v_dis"]
    dA = r["nucleation"] + r["growth"] - r["v_dis"] - r["v_degA"]
    dU.require(FLUX, "dU/dt")
    dA.require(FLUX, "dA/dt")
    return dU, dA


def dimensionalConservationCheck(U_f: Quantity, A_f: Quantity, C_f: Quantity,
                                 D_f: Quantity, dp: DimensionalParams) -> Dict[str, Quantity]:
    """rapid-equilibrium pool balances, written only in free concentrations.

    per DECISIONS.md D004 a total substrate pool must never be substituted into
    a free-resource binding formula; these expressions use free substrate only,
    and the joint solve in `model.solveFreePools` supplies it.
    """
    K_N = Quantity(dp.K_N, CONC)
    K_CU = Quantity(dp.K_CU, CONC)
    K_CA = Quantity(dp.K_CA, CONC)
    K_DU = Quantity(dp.K_DU, CONC)
    K_DA = Quantity(dp.K_DA, CONC)
    N_f = Quantity(dp.N_free, CONC)

    chaperone_occupancy = one() + N_f / K_N + U_f / K_CU + A_f / K_CA
    protease_occupancy = one() + U_f / K_DU + A_f / K_DA
    chaperone_occupancy.require(DIMENSIONLESS, "chaperone occupancy factor")
    protease_occupancy.require(DIMENSIONLESS, "protease occupancy factor")

    return {
        "C_tot": C_f * chaperone_occupancy,
        "D_tot": D_f * protease_occupancy,
        "U_tot": U_f * (one() + C_f / K_CU + D_f / K_DU),
        "A_tot": A_f * (one() + C_f / K_CA + D_f / K_DA),
    }


def dimensionalFreePools(U_tot: float, A_tot: float, dp: DimensionalParams):
    """solve the dimensional pool balances independently of `model.py`.

    deliberately a DIFFERENT formulation and algorithm from
    `model.solveFreePools`: all four free concentrations are unknowns and the
    system is solved by bounded least squares rather than by 2-d newton on the
    eliminated system. agreement between the two is therefore evidence about
    the algebra, not about one implementation agreeing with itself.
    """
    import numpy as np
    from scipy.optimize import least_squares

    def residual(z):
        U_f, A_f, C_f, D_f = z
        return [
            U_f * (1.0 + C_f / dp.K_CU + D_f / dp.K_DU) - U_tot,
            A_f * (1.0 + C_f / dp.K_CA + D_f / dp.K_DA) - A_tot,
            C_f * (1.0 + dp.N_free / dp.K_N + U_f / dp.K_CU + A_f / dp.K_CA) - dp.C_tot,
            D_f * (1.0 + U_f / dp.K_DU + A_f / dp.K_DA) - dp.D_tot,
        ]

    z0 = [min(U_tot, dp.K_CU), min(A_tot, dp.K_CA),
          dp.C_tot / (1.0 + dp.N_free / dp.K_N), dp.D_tot]
    sol = least_squares(residual, z0, bounds=(0.0, np.inf), xtol=1e-15,
                        ftol=1e-15, gtol=1e-15)
    return tuple(float(v) for v in sol.x)


def dimensionalRhsNumeric(U_tot: float, A_tot: float, dp: DimensionalParams):
    """numeric dU/dt, dA/dt in dimensional units, via the unit-carrying rhs."""
    U_f, A_f, C_f, D_f = dimensionalFreePools(U_tot, A_tot, dp)
    dU, dA = dimensionalRhs(Quantity(U_f, CONC), Quantity(A_f, CONC),
                            Quantity(C_f, CONC), Quantity(D_f, CONC), dp)
    return dU.value, dA.value


def toNondimensional(dp: DimensionalParams):
    """map dimensional parameters onto nondimensional `model.Params`.

    concentration scale  phi = C_tot + D_tot  (total rescue machinery)
    time scale           tau = 1 / k_ref      (per-chaperone refolding turnover)

    returns (Params, phi, tau). the scaling exponents are verified by
    `tests/test_units.py::testNondimensionalizationIsExactRescaling`, which
    integrates both models and compares trajectories.
    """
    from .model import Params  # local import avoids a module cycle

    phi = dp.C_tot + dp.D_tot
    tau = 1.0 / dp.k_ref

    p = Params(
        j=dp.J * tau / phi,
        nu=dp.N_free / dp.K_N,
        c_tot=dp.C_tot / phi,
        d_tot=dp.D_tot / phi,
        rho_U=dp.k_U / dp.k_ref,
        rho_A=dp.k_A / dp.k_ref,
        alpha_d=dp.k_dis / dp.k_ref,
        alpha_n=dp.k_n * (phi ** (dp.m - 1.0)) / dp.k_ref,
        alpha_g=dp.k_g * phi / dp.k_ref,
        m=dp.m,
        kappa_ref=dp.K_ref / phi,
        kappa_u=dp.K_U / phi,
        kappa_a=dp.K_A / phi,
        kappa_dis=dp.K_dis / phi,
        kappa_cu=dp.K_CU / phi,
        kappa_ca=dp.K_CA / phi,
        kappa_du=dp.K_DU / phi,
        kappa_da=dp.K_DA / phi,
    )
    return p, phi, tau
