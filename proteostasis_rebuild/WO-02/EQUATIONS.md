# WO-02 — conservative equations

Units uM and s (declared in `WO-01/variables.py`). Code: `model.py`.

```
dN/dt = (1−eps) s_P + phi k_cat B − k_mis N                                 − mu N
dU/dt =  eps s_P + k_mis N + ((1−phi) k_cat + k_off) B − k_on C U
         − k_d U − k_a U² + k_dis A                                          − mu U
dB/dt =  k_on C U − (k_off + k_cat) B                                        − mu B
dC/dt =  s_C + (k_off + k_cat) B − k_on C U                                  − mu C
dA/dt =  k_a U² − (k_dis + k_dA) A                                           − mu A
```

Totals: `P_T = N+U+B+A`, `C_T = C+B`.

```
dP_T/dt = s_P − k_d U − k_dA A − mu P_T
dC_T/dt = s_C − mu C_T          ⇒ C_T → s_C/mu, closed form
```

Balanced growth: `s_P = mu P_T*`, `s_C = mu C_T*`.

## legacy → rebuild mapping

| legacy term (two_pool_ode.py) | rebuild | status |
|---|---|---|
| P, A as proteome fractions (:5–6) | U, A in uM; fractions = U/P_T, A/P_T | modified: denominator is now a state |
| `J_bare` (:14) | `eps s_P`, with `s_P = mu P_T` | modified: ln2 synthesis rate (WO-01) |
| `J_bare(Phi−1)` (:13–14) | none | **rejected**: no donor pool (G2.5) |
| `C_free = C_tot/(1+M/K_d)` (:10) | explicit C, B with mass action | replaced; alternatives in WO-03 |
| `v_fold P` (:11,15) | `phi k_cat B` → N | modified: rescue returns mass to N |
| `k_deg P` (:15) | `k_d U` | retained |
| `drain(1−A_sat)` (:16–19) | `k_a U²` | modified: saturation by A dropped (no mechanism given) |
| `k_clear A` (:20) | `k_dis A` → U and `k_dA A` → degraded | split: disaggregation returns protein |
| none | `−mu x` on every state | added (WO-01 finding) |
| none | `k_mis N` → U | added, default 0 |
| `A_max` gate (:92) | none in dynamics | removed from dynamics; any viability threshold must be an external, measured criterion |
