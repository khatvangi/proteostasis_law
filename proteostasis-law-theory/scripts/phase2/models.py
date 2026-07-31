"""one adapter interface over the two model forms.

the point of this module is that `protocols.py` must be blind to which model it
is solving.  if the boron arm and the free-limit arm went through different
protocol code, the 2x2 factorial would no longer isolate the model-form effect
from the solver effect -- which is the whole reason the factorial exists.

both adapters expose exactly:

    name, epsilon, jScale, params
    rhsVector(x)          -> ndarray(2,)
    jacobian(x)           -> ndarray(2,2)
    numericalJacobian(x)  -> ndarray(2,2)
    freeCoordinates(x)    -> (uf, af)
    residualScale()       -> float
    removalCeiling()      -> float

`freeCoordinates` is what lets the benchmark score BOTH the total-pool and the
free-pool admissibility criteria on every equilibrium from both model forms,
instead of silently comparing a total-burden threshold on one side against a
free-burden threshold on the other.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np

from . import nitrogen_limit
from .mapping import NitrogenParams, nitrogenToBoron


class ModelAdapter:
    """common surface; subclasses supply the model."""

    name = "abstract"

    @property
    def epsilon(self) -> float:
        raise NotImplementedError

    @property
    def jScale(self) -> float:
        raise NotImplementedError

    def rhsVector(self, x) -> np.ndarray:
        raise NotImplementedError

    def jacobian(self, x) -> np.ndarray:
        raise NotImplementedError

    def numericalJacobian(self, x) -> np.ndarray:
        raise NotImplementedError

    def freeCoordinates(self, x) -> Tuple[float, float]:
        raise NotImplementedError

    def residualScale(self) -> float:
        raise NotImplementedError

    def removalCeiling(self) -> float:
        raise NotImplementedError


class BoronAdapter(ModelAdapter):
    """the canonical total-substrate model at finite sequestration epsilon.

    this calls `scripts/proteostasis/model.py` directly -- not a transcription.
    that is deliberate: cell 1 of the factorial is only a code-equivalence test
    if one side really is the shipped boron code path.
    """

    name = "boron"

    def __init__(self, q: NitrogenParams, epsilon: float):
        from proteostasis import model as boron_model   # lazy: nitrogen never needs it
        from proteostasis import equilibria as boron_equilibria
        from . import boron_continuation

        self._model = boron_model
        self._equilibria = boron_equilibria
        self._cont = boron_continuation
        self._epsilon = float(epsilon)
        self.nitrogen_params = q
        self.params = boron_model.Params(**nitrogenToBoron(q, epsilon)).validate()
        #: how many evaluations left the shipped nonnegative code path.  the
        #: benchmark records this so a cell can never claim to be a pure boron
        #: code path when it was not.
        self.continuation_evaluations = 0

    @property
    def epsilon(self) -> float:
        return self._epsilon

    @property
    def jScale(self) -> float:
        return self.params.j

    def _count(self, u: float, a: float) -> None:
        if self._cont.continuationUsed(u, a):
            self.continuation_evaluations += 1

    def rhsVector(self, x) -> np.ndarray:
        u, a = float(x[0]), float(x[1])
        self._count(u, a)
        return self._cont.rhsVectorExtended((u, a), self.params)

    def jacobian(self, x) -> np.ndarray:
        u, a = float(x[0]), float(x[1])
        self._count(u, a)
        return self._cont.jacobianExtended(u, a, self.params)

    def numericalJacobian(self, x) -> np.ndarray:
        return self._model.numericalJacobian(float(x[0]), float(x[1]), self.params)

    def freeCoordinates(self, x) -> Tuple[float, float]:
        uf, af, _, _ = self._cont.solveFreePoolsExtended(
            float(x[0]), float(x[1]), self.params)
        return float(uf), float(af)

    def residualScale(self) -> float:
        return self._equilibria.residualScale(self.params)

    def removalCeiling(self) -> float:
        return self._model.removalCeiling(self.params)


class FreeLimitAdapter(ModelAdapter):
    """the epsilon = 0 face: nitrogen's free-substrate model.

    `epsilon` is still carried, because the free-limit model is evaluated at the
    SAME binding constants as its boron counterpart at each rung of the ladder.
    those constants still act (they titrate the resource pools); what is absent
    is the substrate sequestration factor.  carrying epsilon here is therefore
    not a contradiction -- it labels which boron cell this cell is matched to.
    """

    name = "free"

    def __init__(self, q: NitrogenParams, epsilon: float):
        self._epsilon = float(epsilon)
        self.nitrogen_params = q
        self.params = q

    @property
    def epsilon(self) -> float:
        return self._epsilon

    @property
    def jScale(self) -> float:
        return self.params.j

    def rhsVector(self, x) -> np.ndarray:
        return nitrogen_limit.rhsVector(x, self.params)

    def jacobian(self, x) -> np.ndarray:
        return nitrogen_limit.jacobian(float(x[0]), float(x[1]), self.params)

    def numericalJacobian(self, x) -> np.ndarray:
        return nitrogen_limit.numericalJacobian(float(x[0]), float(x[1]), self.params)

    def freeCoordinates(self, x) -> Tuple[float, float]:
        return float(x[0]), float(x[1])

    def residualScale(self) -> float:
        return self.params.residualScale()

    def removalCeiling(self) -> float:
        return self.params.removalCeiling()


def makeAdapter(model_form: str, q: NitrogenParams, epsilon: float) -> ModelAdapter:
    """`model_form` is 'boron' or 'free'.

    the free-limit arm is built from the SAME nitrogen parameter vector `q`, with
    the epsilon-dependent binding coefficients already written into `q` by
    `lhs.parametersForEpsilon`.  so both arms see identical binding constants and
    identical kinetics; the only difference is the sequestration factor.
    """
    if model_form == "boron":
        return BoronAdapter(q, epsilon)
    if model_form == "free":
        return FreeLimitAdapter(q, epsilon)
    raise ValueError(f"unknown model form {model_form!r}; expected 'boron' or 'free'")
