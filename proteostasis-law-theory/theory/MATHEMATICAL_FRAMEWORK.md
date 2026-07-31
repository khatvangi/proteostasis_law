# Mathematical Framework

## Translation strategy and objectives

Let `u in U(e)` denote a translation strategy in environment `e`. Map it to

`Phi(u,e) = (T, F, Y, V, -G, -Q)`,

where `T` is productive throughput, `F` decoding fidelity, `Y` folding yield, `V` functional/chemical diversity supported, `G` energetic cost, and `Q` quality-control demand. Orienting costs negatively lets all coordinates be maximized. The Pareto set contains strategies for which no attainable alternative improves one coordinate without worsening another.

## Site-resolved influx

For protein `p`, site `i`, cognate residue `a(p,i)`, encoded codon `c(p,i)`, translation flux `E_p(u,e)` (molecules volume^-1 time^-1), site criticality `kappa_{p,i}`, and delivered residue `r`, define

`J_U(u,e) = sum_p E_p(u,e) sum_i kappa_{p,i} sum_{r != a(p,i)} Pr(r | c(p,i), u, e) d(r,p,i,e)`.

Here `d` maps a specific substitution at a specific site into expected entry into the modeled damaged/misfolded pool per translated chain. If `kappa` and `d` overlap conceptually, one must be a dimensionless modifier and the other a calibrated burden contribution; their definitions must prevent double counting. Correctly translated chains that fail to fold may be included through a separate influx `J_fold(u,e)`. The total abnormal monomer input is `J = J_U + J_fold + J_other`, while ordinary productive nascent-chain occupancy is handled separately in resource conservation.

The global mean mistranslation rate

`epsilon_bar = [sum_p E_p sum_i sum_{r != a} Pr(r|c,u,e)] / [sum_p E_p L_p]`

is an explicitly **coarse approximation**. It cannot generally determine `J_U`, because identical values can conceal different expression weights, substitution identities, site criticalities, and damage functions.

## Dynamics and feasibility

Let `x=(U,A,...)` and `dx/dt=F(x;J,N,theta)`, where `N` denotes ordinary nascent-chain load. For thresholds `H`, define a safe region `D_H`. A strategy is feasible if its loads and allocations yield a bounded attracting invariant state inside `D_H` for the relevant initial conditions.

## Dimension and uncertainty rules

- Every flux has amount volume^-1 time^-1; every pool has amount volume^-1.
- Terms may be added only when units and biological meanings match.
- Probabilities, criticality weights, and damage conversions must be identified separately.
- Parameter and measurement uncertainty propagate to a probability or interval of feasibility; a point classification is not assumed exact.

## Logical status

Pareto membership and feasibility are definitions. Stability results become theorem-level only after assumptions on `F` are stated and proved. Parameter-dependent behavior is a model consequence. Correspondence between the model and living systems is an empirical hypothesis.

