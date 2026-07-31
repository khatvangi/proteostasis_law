"""the deterministic 7-D latin hypercube, matched to nitrogen's sweep.

the sampler and the seven inverse transforms are transcribed from
nitrogen:.../src/lhs_sweep.py:33-48.  they are reproduced rather than imported
so the benchmark runs on either host, and the reproduction is tested against a
pinned sample-matrix hash (`SAMPLE_HASHES`) that was verified to be BIT-IDENTICAL
on boron (scipy 1.14.0 / numpy 2.2.6) and nitrogen (scipy 1.15.2 / numpy 1.26.4).
that cross-host check matters: if the two hosts disagreed on the sample matrix,
every cross-host cell comparison would be comparing different parameter draws.

the other sixteen coordinates are pinned -- see `mapping.PINNED_NITROGEN` for the
nitrogen-side values and `mapping.D_TOT_GAUGE` for the one gauge freedom that the
seven sampled groups do not fix.
"""

from __future__ import annotations

import hashlib
from typing import Dict, Sequence

import numpy as np
from scipy.stats import qmc

from .mapping import (D_TOT_GAUGE, LHS_COORDINATES, PINNED_NITROGEN,
                      NitrogenParams, bindingConstants)

DEFAULT_SEED = 20260731
DEFAULT_N = 2000

#: sha256 of the float64 sample matrix, verified identical on boron and nitrogen.
SAMPLE_HASHES: Dict[int, str] = {
    8: "1b21fdccaccd72796a324afb53366d5de31b563ebddc56933353b8213b302589",
    2000: "a937f4bc68faa9bb404bfa81c190a168f19a866d5402954a691e9b2d650d85fa",
}


def sampleMatrix(n: int = DEFAULT_N, seed: int = DEFAULT_SEED) -> np.ndarray:
    """the raw unit-hypercube sample.  identical call to nitrogen's."""
    return qmc.LatinHypercube(d=7, seed=seed).random(n)


def sampleHash(matrix: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(matrix, dtype=np.float64).tobytes()).hexdigest()


def _logmap(q: float, lo: float, hi: float) -> float:
    return float(np.exp(np.log(lo) + q * (np.log(hi) - np.log(lo))))


def sampledCoordinates(z: Sequence[float]) -> Dict[str, float]:
    """one hypercube row -> the seven nitrogen groups it controls."""
    out: Dict[str, float] = {}
    for k, (name, kind, lo, hi) in enumerate(LHS_COORDINATES):
        q = float(z[k])
        out[name] = _logmap(q, lo, hi) if kind == "log" else float(lo + q * (hi - lo))
    return out


def parametersForEpsilon(z: Sequence[float], epsilon: float) -> NitrogenParams:
    """the matched nitrogen parameter vector at one hypercube row and one epsilon.

    both arms of the factorial are built from THIS object: the free-limit arm
    uses it directly, and the boron arm passes it through
    `mapping.nitrogenToBoron`.  the binding coefficients written here and the
    kappas that `nitrogenToBoron` derives are exact reciprocals by construction,
    which `tests/phase2/test_mapping.py` asserts rather than assumes.
    """
    coords = sampledCoordinates(z)
    c_tot_boron = coords["ref_capacity"] * PINNED_NITROGEN["c_tot"]   # mu = 1
    kb = bindingConstants(epsilon, c_tot_boron, D_TOT_GAUGE)
    kw = dict(PINNED_NITROGEN)
    kw.update(coords)
    kw.update({"c_u": kb["c_u"], "c_a": kb["c_a"], "d_u": kb["d_u"], "d_a": kb["d_a"]})
    return NitrogenParams(**kw)


def describeDesign(n: int, seed: int) -> Dict:
    """the design record written into every benchmark manifest."""
    m = sampleMatrix(n, seed)
    h = sampleHash(m)
    return {
        "n_samples": int(n),
        "seed": int(seed),
        "sample_matrix_sha256": h,
        "sample_matrix_hash_pinned_match": SAMPLE_HASHES.get(n) == h if n in SAMPLE_HASHES else None,
        "sampled_coordinates": [
            {"name": nm, "transform": kind, "lo": lo, "hi": hi}
            for nm, kind, lo, hi in LHS_COORDINATES
        ],
        "pinned_nitrogen": dict(PINNED_NITROGEN),
        "d_tot_gauge": D_TOT_GAUGE,
    }
