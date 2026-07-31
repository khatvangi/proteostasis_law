"""phase 1 computational falsification package for the proteostasis law.

these modules test the INTERNAL mathematical robustness of the minimal
conserved-resource model defined in `theory/DYNAMICAL_SYSTEM.md`. nothing here
calibrates the model to an organism, and nothing here validates the theory
against data. the outputs are statements about the model, conditional on its
kinetic forms and parameter ranges -- 'model consequences' in the claim
taxonomy of `theory/SCOPE_AND_NONCLAIMS.md`.
"""

from .model import (Params, ModelError, rhs, rhsVector, jacobian, numericalJacobian,
                    fluxes, solveFreePools, solveFreePoolsCertified,
                    poolBalanceResiduals, massBalanceResidual, removalCeiling,
                    influxAdmitsNoBoundedState, allocationParams, paramsFromDict)
from .equilibria import (Equilibrium, findEquilibria, lowestStableEquilibrium,
                         solveEquilibriumFrom, traceBranch, findFold,
                         eigenvalueAlongBranch, FoldResult)
from .simulate import Trajectory, simulate, basinScan, defaultInitialConditions, recoveryTime

__all__ = [
    "Params", "ModelError", "rhs", "rhsVector", "jacobian", "numericalJacobian",
    "fluxes", "solveFreePools", "solveFreePoolsCertified", "poolBalanceResiduals",
    "massBalanceResidual", "removalCeiling", "influxAdmitsNoBoundedState",
    "allocationParams", "paramsFromDict", "Equilibrium", "findEquilibria",
    "lowestStableEquilibrium", "solveEquilibriumFrom", "traceBranch", "findFold",
    "eigenvalueAlongBranch", "FoldResult", "Trajectory", "simulate", "basinScan",
    "defaultInitialConditions", "recoveryTime",
]
