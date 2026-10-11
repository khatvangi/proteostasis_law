# Positive invariance and uniqueness (fixed production, constant growth)

Model as in MINIMAL_DYNAMICS.md with τ = 0, β0 > 0, K > 0, h strictly between 0 and 1 for finite
positive F.

## The box

```
B = { 0 ≤ U ≤ Umax = L/(μ+δ),   Cmin = β0/μ ≤ C ≤ Cmax = (β0+βm)/μ }
```

On each edge the vector field points strictly inward:

| Edge | Sign | Why |
|---|---|---|
| U = 0 | dU/dt = L > 0 | no client to serve or dilute |
| U = Umax | dU/dt = −kCU/(K+U) < 0 | C ≥ Cmin > 0 |
| C = Cmin | dC/dt = βm·h(F) > 0 | h > 0 for finite F |
| C = Cmax | dC/dt = βm·(h(F) − 1) < 0 | h < 1 because F > 0 |

Tested on a 21-point grid along every edge for 100 random placeholder draws.

## Consequences

1. Every trajectory starting in B stays in B; every trajectory starting outside with U ≥ 0 enters
   the C-band (dC/dt has the right sign outside it) and then B. Pools and client stay bounded.
2. The vector field on B is inward and has no zeros on the boundary, so by Poincaré–Hopf the indices
   of the equilibria inside sum to +1. Every equilibrium has det J > 0 (MINIMAL_DYNAMICS §1), i.e.
   index +1. Hence **exactly one equilibrium, and it is locally stable.** Confirmed numerically: one
   root in 200/200 draws.
3. The divergence is negative everywhere, not just at equilibria:
   ∂(dU/dt)/∂U + ∂(dC/dt)/∂C = −kCK/(K+U)² − μ − δ + βm·h′(F)·K/(K+U) − μ < 0, since h′ ≤ 0.
   By Bendixson's criterion there is no closed orbit in B; B is compact and holds one equilibrium,
   which is a stable node or focus (so no homoclinic loop). Poincaré–Bendixson then gives
   **global convergence within B when τ = 0.** Tested: numeric divergence < 0 at random points of B.

## What it does not give

- **Not survival.** Boundedness of U ≤ L/(μ+δ) is bookkeeping: if L is large, U* is large and
  "bounded" means a large standing pool of unfolded protein. Viability needs U* below a tolerance,
  which the model does not supply.
- **Not robust to delay.** With τ > 0 the box is still invariant for the undelayed C-band argument
  (h stays in (0,1) whatever F was), but the equilibrium can lose stability (MINIMAL_DYNAMICS §3);
  trajectories then oscillate inside B.
- **Not robust to load-dependent growth.** If μ falls with U, Umax = L/(μ+δ) is no longer a bound
  (dilution flux can fall as U rises), the box argument fails, and multiple equilibria appear
  (MINIMAL_DYNAMICS §5; tested).
