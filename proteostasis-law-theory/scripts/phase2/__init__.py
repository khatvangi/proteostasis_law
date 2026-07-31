"""phase 2A -- matched boron/nitrogen equivalence benchmark.

phase 1 asked whether the boron model is internally consistent.  phase 2A asks a
different and prior question: are the boron and nitrogen phase-1 implementations
even models of the same vector field, and if not, how far apart are they?

the read-only cross-implementation audit answered the first half: nitrogen's
free-substrate model is exactly boron's total-substrate model at sequestration
strength epsilon = 0.  this package makes that claim executable.

  `mapping`          the parameter map and the epsilon ladder
  `nitrogen_limit`   locally transcribed epsilon = 0 right-hand side
  `models`           one adapter interface over both model forms
  `protocols`        the two root-search + attraction protocols, model agnostic
  `lhs`              the deterministic 7-D latin hypercube, matched to nitrogen
  `t0_equivalence`   T0: the epsilon -> 0 identity and O(epsilon) scaling test
  `run_matched_benchmark`  the 2x2 factorial benchmark driver

nothing in this package validates the proteostasis law.  every output is a
statement about two pieces of code.
"""
