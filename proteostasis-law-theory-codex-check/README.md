# Independent Phase 1 falsification workspace

This directory is an independent numerical replication of the minimal finite-resource
proteostasis model described in the canonical theory. It does not import or modify
canonical code and it is not a validation claim.

The nondimensional states `U` and `A` are free, binding-competent soluble and aggregate
burdens. Conserved free resources are

```
C_f = C_tot / (1 + N/K_N + U/K_CU + A/K_CA)
D_f = D_tot / (1 + U/K_DU + A/K_DA)
```

so ordinary nascent occupancy `N/K_N` is separate from abnormal influx `J`. The rate
laws and signs match the canonical two-state baseline. At both coordinate boundaries
the vector field points into the nonnegative quadrant, giving an analytic positivity
check; numerical integrations also terminate and fail if a materially negative state
is observed.

Run locally with:

```bash
python3 -m unittest discover -s tests -v
python3 sweep.py --config config/phase1.json --output results/phase1
```

Sweep labels are operational numerical falsification categories, not biological
validation. In particular, `no_bounded_attracting_state` means no stable equilibrium
was found by the specified deterministic search and all prescribed trajectories met
the escape criterion. It is not a global proof excluding cycles or remote attractors.

