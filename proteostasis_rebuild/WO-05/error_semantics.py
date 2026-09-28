"""
WO-05 error-to-burden semantics.

three quantities, kept apart by type (all dimensionless, per codon translated):

  RAW   e_raw   probability that decoding a codon incorporates a non-cognate
                tRNA, synonymous outcomes included. not measured by MS.
  SUB   f_sub   probability that the incorporated amino acid differs from the
                encoded one, summed over all destination amino acids.
                under the legacy uniform-S assumption f_sub = e_raw (1 - S).
  MS    r_MS    Landerer 2024 "error detection rate": per codon and dataset,
                PSMs covering the codon position that carry a substitution
                there / all PSMs covering it (Methods, "Calculating Error
                Detection Rates", PMC10939442). it is amino-acid level (MS
                cannot see synonymous errors) and summed over destinations,
                so it estimates f_sub, from below: I/L substitutions are
                invisible (identical mass) and rare events are missed.
                Data_S2 "mean" is the mean of per-dataset rates over the
                datasets in which the codon had >= 1 detected substitution
                (inferred: n = (sd/se)^2 is an integer 3..68 for every codon,
                never the 80 datasets).

burden flux into the misfolded pool (legacy units, fraction/s):
  J = f_sub * N_prot * p_misfold / T_gen
    = e_raw * N_prot * (1 - S) * p_misfold / T_gen
"""
import math
from dataclasses import dataclass

RAW, SUB, MS = "raw_decoding", "aa_substitution", "ms_detection"


@dataclass(frozen=True)
class ErrorRate:
    value: float
    kind: str               # RAW | SUB | MS

    def __post_init__(self):
        if self.kind not in (RAW, SUB, MS):
            raise ValueError(f"unknown error kind {self.kind!r}")
        if not (0.0 <= self.value <= 1.0):
            raise ValueError("an error probability per codon must lie in [0, 1]")


def substitution_from_raw(e: ErrorRate, S: float) -> ErrorRate:
    """the one place (1 - S) is allowed: raw decoding error -> substitution."""
    if e.kind != RAW:
        raise TypeError(f"(1-S) applies only to raw decoding error, got {e.kind}")
    return ErrorRate(e.value * (1.0 - S), SUB)


def substitution_from_ms(r: ErrorRate, detectability: float = 1.0) -> ErrorRate:
    """MS detection rate -> substitution rate. detectability is the fraction of
    substitutions MS can see (<= 1, unmeasured); 1.0 gives a LOWER BOUND on
    f_sub. no (1 - S) here: MS already sees only amino-acid changes."""
    if r.kind != MS:
        raise TypeError(f"expected an MS detection rate, got {r.kind}")
    if not (0.0 < detectability <= 1.0):
        raise ValueError("detectability must lie in (0, 1]")
    return ErrorRate(min(r.value / detectability, 1.0), SUB)


def as_substitution(e: ErrorRate, S: float) -> ErrorRate:
    if e.kind == RAW:
        return substitution_from_raw(e, S)
    if e.kind == MS:
        return substitution_from_ms(e)
    return e


def flux(e: ErrorRate, N_prot: float, p_misfold: float, T_gen: float, S: float) -> float:
    """corrected mapping: every kind is first converted to f_sub exactly once."""
    return as_substitution(e, S).value * N_prot * p_misfold / T_gen


def flux_legacy(f: float, N_prot: float, p_misfold: float, T_gen: float, S: float) -> float:
    """the legacy mapping (two_pool_ode.py:513 and scripts 09/11/12): always
    multiplies by (1 - S). correct only if f is a raw decoding error."""
    return f * N_prot * (1.0 - S) * p_misfold / T_gen


def flux_checked_legacy(e: ErrorRate, N_prot, p_misfold, T_gen, S):
    """the legacy formula, guarded: refuses a substitution-level input."""
    if e.kind != RAW:
        raise TypeError(f"legacy (1-S) flux applied to {e.kind} input: double discount")
    return flux_legacy(e.value, N_prot, p_misfold, T_gen, S)


# ---------------------------------------------------------------- thresholds
def threshold_raw(N, P_correct, S, p_misfold):
    """exact per-protein arithmetic threshold on e_raw: (1-P^(1/N)) / ((1-S) p_m)."""
    return (1.0 - P_correct ** (1.0 / N)) / ((1.0 - S) * p_misfold)


def threshold_sub(N, P_correct, p_misfold):
    """same threshold stated on f_sub (no synonymous filter)."""
    return (1.0 - P_correct ** (1.0 / N)) / p_misfold


def threshold_legacy_quoted(N, P_correct):
    """what the manuscript quotes as 1.19e-3: -ln(P)/N, i.e. the factor
    (1-S) p_m forced to 1.0 (arithmetic_stress_test.py:202-205)."""
    return -math.log(P_correct) / N
