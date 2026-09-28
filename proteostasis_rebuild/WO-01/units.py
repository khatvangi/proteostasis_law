"""
minimal dimensional checker for sympy expressions.

a unit is (scale, dims): the physical quantity equals numeric_value * scale *
base, with base dimensions conc (molar, M) and time (s), plus two count
dimensions, codon and protein, so that "per codon" and "per protein" cannot be
silently exchanged. uM therefore has scale 1e-6.

rules enforced by unit_of():
  - Add: every term must have identical dims AND identical scale. identical dims
    with different scale (uM + M) is a numeric bug and raises.
  - Mul/Pow: dims add / multiply, scales multiply / power.
  - functions (sqrt, exp, log): exp/log require a dimensionless, scale-1 argument.
  - bare numbers are dimensionless with scale 1. unit conversion must therefore
    be written with an explicit conversion symbol (e.g. UM_TO_M), never as a
    naked 1e-6, so the checker can see it.
"""
from dataclasses import dataclass
import math

import sympy as sp

BASE = ("conc", "time", "codon", "protein")


class DimensionError(Exception):
    pass


@dataclass(frozen=True)
class Unit:
    scale: float = 1.0
    dims: tuple = (0, 0, 0, 0)

    def __mul__(self, o):
        return Unit(self.scale * o.scale, tuple(a + b for a, b in zip(self.dims, o.dims)))

    def __truediv__(self, o):
        return Unit(self.scale / o.scale, tuple(a - b for a, b in zip(self.dims, o.dims)))

    def __pow__(self, k):
        return Unit(self.scale ** float(k), tuple(a * k for a in self.dims))

    def same(self, o):
        return self.dims == o.dims and math.isclose(self.scale, o.scale, rel_tol=1e-12)

    @property
    def dimensionless(self):
        return all(d == 0 for d in self.dims) and math.isclose(self.scale, 1.0)

    def __str__(self):
        parts = [f"{n}^{d}" if d != 1 else n for n, d in zip(BASE, self.dims) if d != 0]
        s = " ".join(parts) if parts else "1"
        return s if math.isclose(self.scale, 1.0) else f"{self.scale:g}*{s}"


def U(scale=1.0, conc=0, time=0, codon=0, protein=0):
    return Unit(scale, (conc, time, codon, protein))


# named units used across the rebuild
ONE = U()
PER_S = U(time=-1)
S = U(time=1)
M = U(conc=1)
UM = U(1e-6, conc=1)
UM_PER_S = U(1e-6, conc=1, time=-1)
PER_UM_PER_S = U(1e6, conc=-1, time=-1)
PER_M_PER_S = U(conc=-1, time=-1)
PER_CODON = U(codon=-1)
CODON_PER_PROTEIN = U(codon=1, protein=-1)
M_PER_UM = U(1e6)          # conversion factor: numeric value 1e-6, unit M/uM


def unit_of(expr, units):
    """return the Unit of a sympy expression; raise DimensionError if inconsistent."""
    if isinstance(expr, sp.Symbol):
        if expr.name not in units:
            raise DimensionError(f"no unit declared for symbol {expr.name}")
        return units[expr.name]
    if expr.is_Number:
        return ONE
    if isinstance(expr, sp.Add):
        us = [unit_of(a, units) for a in expr.args]
        for u in us[1:]:
            if not u.same(us[0]):
                raise DimensionError(f"Add mixes {us[0]} and {u} in {expr}")
        return us[0]
    if isinstance(expr, sp.Mul):
        out = ONE
        for a in expr.args:
            out = out * unit_of(a, units)
        return out
    if isinstance(expr, sp.Pow):
        base, ex = expr.args
        if not ex.is_Number:
            if not unit_of(ex, units).dimensionless:
                raise DimensionError(f"exponent not dimensionless in {expr}")
            if not unit_of(base, units).dimensionless:
                raise DimensionError(f"symbolic power of dimensional base in {expr}")
            return ONE
        return unit_of(base, units) ** sp.Rational(ex)
    if isinstance(expr, (sp.exp, sp.log)):
        u = unit_of(expr.args[0], units)
        if not u.dimensionless:
            raise DimensionError(f"{expr.func.__name__} of dimensional argument {u}")
        return ONE
    raise DimensionError(f"unsupported node {type(expr).__name__}: {expr}")


def check_equation(lhs_unit, rhs, units):
    """rhs must reduce to lhs_unit; returns the list of (term, unit) for reporting."""
    report = [(t, unit_of(t, units)) for t in (rhs.args if isinstance(rhs, sp.Add) else [rhs])]
    u = unit_of(rhs, units)
    if not u.same(lhs_unit):
        raise DimensionError(f"rhs has unit {u}, expected {lhs_unit}")
    return report
