# WO-01 — variables, units, conservation

Machine-readable source of truth: `variables.py` (imported by later WOs).
Units: concentrations uM, time s. Checker: `units.py` (tracks scale, so uM+M fails).

## states (all uM)

| symbol | meaning |
|---|---|
| N | native, functional protein (monomer units) |
| U | free non-native monomer (misfolded or error-bearing, chaperone-free) |
| B | chaperone–client complex (one client + one chaperone) |
| A | aggregated protein, monomer-equivalents |
| C | free chaperone (effective folding-competent pool) |

Derived totals: `P_T = N + U + B + A` (protein), `C_T = C + B` (chaperone).
Fractions are reported only as explicit ratios (`A/P_T`, `U/P_T`), never as states.

## parameters

| symbol | unit | meaning |
|---|---|---|
| s_P | uM/s | total protein synthesis flux |
| eps | 1 | probability a newly made chain enters U |
| s_C | uM/s | chaperone synthesis flux |
| mu | 1/s | growth = dilution rate, ln2 / doubling time |
| k_on | 1/(uM s) | chaperone–client association |
| k_off | 1/s | unproductive release |
| k_cat | 1/s | completion of one chaperone cycle |
| phi | 1 | P(completed cycle yields native) |
| k_d | 1/s | degradation of U |
| k_a | 1/(uM s) | aggregation; monomer flux `k_a U²` |
| k_dis | 1/s | disaggregation A → U |
| k_dA | 1/s | degradation of A |
| k_mis | 1/s | spontaneous unfolding N → U |

## conservation laws (proved symbolically in WO-02)

```
dP_T/dt = s_P − k_d U − k_dA A − mu P_T      sources: synthesis
                                             sinks: degradation of U, A; dilution
dC_T/dt = s_C − mu C_T                       source: chaperone synthesis; sink: dilution
```

Every flux inside the network (binding, release, folding, aggregation,
disaggregation, unfolding) must appear with opposite signs in exactly two pool
equations, so it cancels in the totals.
