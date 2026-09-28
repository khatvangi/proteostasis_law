# WO-01 — variables, units, conservation

## gates (fixed before analysis)

- G1.1 every rebuild state/parameter has a declared unit; no fraction state.
- G1.2 automated dimensional check of every legacy term, with its unit reported.
- G1.3 conservation laws stated with every source and sink.
- G1.4 each legacy term classified: dimensionally consistent? conservative flux?

## what was done

- `units.py`: sympy dimension walker. A unit is (scale, dims) over
  conc/time/codon/protein; addition requires equal dims AND equal scale.
- `variables.py` + `VARIABLES.md`: the rebuild's 5 states (all uM) and 13
  parameters, and the two conservation laws.
- `legacy_units_audit.py`: the legacy two-pool pieces (two_pool_ode.py
  lines 17–20, 128–159, 262) transcribed as sympy, with the legacy's literal
  `1e-6` written as an explicit conversion symbol so a uM/M slip would show.
  Output: `legacy_units_audit.json`.

## results

**Dimensions (G1.2).** All 11 legacy pieces are dimensionally consistent,
including the uM→M conversions in `v_agg` and `drain`. A deliberate uM + M
negative control is rejected, so the pass is not vacuous.

**Bookkeeping (G1.4).** Dimensionally consistent is not the same as conservative:

| legacy flux | donor → receiver | conservative? |
|---|---|---|
| `J_bare` | synthesis → P | yes, as an external source (synthesis not modelled) |
| `J_bare(Phi−1)` | **none** → P | **no**: inflow scales with `v_agg`, no pool debited, no process named |
| `k_deg P` | P → degraded | yes |
| `v_fold P` | P → untracked native pool | yes (harmless while native is not a state) |
| `drain(1−A_sat)` | P → A | yes |
| `k_clear A` | A → removed | yes as written; biologically disaggregation returns A to P or N |
| dilution `mu P`, `mu A` | — | **absent** |

**Two findings beyond the prior audit.**

1. *Synthesis-rate factor.* `J_bare = f N (1−S) p / T_gen` uses `1/T_gen` as
   the synthesis rate per existing protein. In balanced exponential growth that
   rate is `mu = ln2 / T_gen`. The legacy mapping therefore overstates inflow by
   `1/ln2 = 1.4427` at any T_gen. (Integrated over a generation, `1/T_gen` is
   synthesis relative to the *starting* proteome; an instantaneous ODE needs the
   rate relative to the *current* proteome, which is ln2/T_gen.)
2. *Dilution is missing.* At T_gen = 3600 s, `mu = 1.93e-4 /s` — comparable to
   `k_deg = 3e-4` and `k_clear = 4e-4`. Neither P nor A is diluted by growth. If
   `k_deg` was meant to include dilution, A still has none.

Both enter WO-02 as corrections: synthesis `s_P = mu P_T` (balanced growth, no
turnover term assumed) and `−mu x` on every state.

## tests

`python -m unittest test_wo01 -v` → 13 tests, OK (checker positives/negatives,
declarations, legacy audit, Phi non-conservation, ln2 factor).

## adversarial self-review

- *Is the legacy dimension pass meaningful?* Only partly. P and A are declared
  dimensionless, which makes the fraction-based equations easy to pass. The real
  defect of fraction states is not dimensional: it is that the denominator
  (total protein) is never a state, so nothing forces `P + A ≤ 1`. That is why
  G1.4 exists, and it catches the non-conservation that G1.2 cannot.
- *Could `1/T_gen` be intended?* It would be right for "fraction of the
  initial proteome synthesized per generation", but it is used as an
  instantaneous rate inside `dP/dt`. The factor is a correction, not a matter of
  taste. Its direction (overstating inflow) makes the legacy headroom
  *conservative*, so correcting it helps rather than hurts the old headline;
  that is recorded, not suppressed.
- *Did I assign units the legacy never stated?* Yes, `f_codon` as proteins
  affected per codon and `N_prot` as codons per protein. These are the only
  readings under which line 262 is dimensionally consistent; they are my
  reading, labelled as such.
- *Are the rebuild parameters complete?* For the core model yes; WO-03 adds
  nascent-chain species and must declare them with the same checker.

VERDICT: PASS — G1.1, G1.2, G1.3, G1.4 all met.
