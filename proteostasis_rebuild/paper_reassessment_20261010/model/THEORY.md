# A finite-capacity route scaffold

A client class k has input flux d_k≥0. A qualified route r has flux x_kr≥0, one terminal outcome (productive rescue, degradation, or sequestration), machine service requirement a_mkr≥0, and ATP cost e_kr≥0. A serial route consumes **every** listed machine; alternative routes are distinct columns. A route such as TF→GroEL/GroES and a route DnaK/DnaJ/GrpE→GroEL/GroES share their downstream resource. A TF-only route cannot replace either for an obligate client.

Let H map routes to classes, A map routes to machine work, and P_k select productive routes of class k. The physical constraints are

    Hx + u = d,  u≥0;     Ax≤C,  x≥0;
    P_k x≥q_k            for each essential client class k.

A client exits exactly once; degradation does not count as rescue. Each capacity C_m = measured pool_m × unknown effective turnover_m is in machine-work/time. Each a_mkr is machine work/client; their sums require a declared common work convention. Rescaling a machine's work unit rescales its C_m and changes a naive scalar total: pooling is not a unit-free physical operation. The experiments fix a convention and vary unknown costs and turnover explicitly. Effective service requirements can absorb expected repeated cycles; the benchmark does not infer their values.

Sequestration is not an infinite steady-state sink: its capacity is an explicit horizon-normalized remaining storage budget. A continuing steady state needs storage release or growth dilution, neither modeled here. The LP assumes a fixed environment, complete route enumeration, divisible stationary fluxes, and no kinetic timing, finite intermediate pools, or spontaneous adaptation. Route order identifies a qualified itinerary, but A alone cannot distinguish two itineraries with identical service requirements in different order. This is necessary resource bookkeeping, not a kinetic proof of serial order.

## Exact dual and the meaning of a cut

For nonnegative route reward w, with no essential lower bounds, maximize wᵀx subject to Ax≤C, Hx≤d, x≥0. Its exact dual is

    minimize Cᵀp + dᵀz
    subject to Aᵀp + Hᵀz ≥ w,   p,z ≥ 0.

Weak duality follows by multiplying both primal inequalities by p,z; finite-dimensional LP strong duality gives equality since x=0 is feasible and H bounds every route. The right side is a fractional resource/demand cover: each productive route must be covered by its machine prices plus its class-demand price. Complementarity gives p_m(C_m−(Ax)_m)=0 and x_kr((Aᵀp+Hᵀz)_kr−w_kr)=0. At differentiable capacity points, ∂V/∂C_m=p_m; at a kink optimal dual prices are supergradients of the concave value function V(C). Prices can be nonunique. Prices are conditional on the objective, demand, and any imposed quotas, not intrinsic machine importance.

For essential quotas P x≥q, introduce v≥0. The dual becomes

    minimize Cᵀp + dᵀz − qᵀv
    subject to Aᵀp + Hᵀz − Pᵀv ≥ w,  p,z,v≥0.

Equivalently, essential feasibility has this exact Farkas certificate: for every p,z,v≥0 satisfying Aᵀp+Hᵀz−Pᵀv≥0, one must have Cᵀp+dᵀz−qᵀv≥0. A negative value certifies infeasibility. Full service Hx=d is not needed: unserved nonessential clients remain explicit.

For productive routes only, maximum common demand scale λ obeys Ax≤C and Hx=λd. Its exact dual is

    λ* = min Cᵀp
    subject to Aᵀp ≥ Hᵀy,  dᵀy≥1,  p≥0, y unrestricted.

The primal variable λ≥0 gives the normalization inequality. For positive finite optimum one can choose equality by scaling. This is the radial capacity used to set benchmark loads. It is an all-client folding margin, not a growth rate.

**Ordinary min-cut is a special case.** When every route uses exactly one eligible machine, with unit service per client and productive outcome, build source→class arcs of capacity d_k, class→eligible-machine arcs of infinite capacity, and machine→sink arcs of capacity C_m. Then

    max total rescue = min_{S subset of classes} [sum_{k not in S} d_k + sum_{m in N(S)} C_m].

All demand is feasible iff d(S)≤C(N(S)) for every S (capacitated Hall condition). General serial routes with arbitrary costs do not have this ordinary graph min-cut identity. Treating a serial route as a collection of interchangeable alternatives destroys the AND requirement; claiming ordinary integral min-cut for the general benchmark would be incorrect. The LP dual above is the exact generalization.

## Non-identification theorem for one pool

Fix the scalar observable S=Σ_m C_m, total demand, and even the vector of client demands. There is no function of these observables alone that recovers essential productive feasibility over the class of route models.

Proof by witness: two machines have capacities (1,0), and two classes have demands (1/2,1/2), the first essential. In architecture A, both classes have productive routes using machine 1 alone. Both demands are feasible. In architecture B, the essential class requires machine 2 and the other uses machine 1. Essential folding is impossible. Capacities, demands, and scalar observables coincide. The route architecture is therefore not identifiable from the pooled observables. Even with architecture B fixed, swapping capacities to (0,1) keeps S=1 but changes essential deliverability from zero to full (with the other class unserved). A scalar model must predict identical outcomes for these matched totals and hence cannot be correct in both. □

This is a non-identification theorem, not a claim that every scalar approximation is useless. Pooling is exact, for example, for unit-cost single-machine routes where every client can use every machine. Nor does the theorem imply any particular biological failure frequency or prove that GroEL is the dominant bottleneck.

## Three strictly defined nested models

Keep the same route columns, rewards, client demands, and costs; only change resource inequalities.

1. **Scalar:** (1ᵀA)x≤1ᵀC. Class demand balances remain bookkeeping, but machine identity has no effect at fixed total remaining capacity. Minimum route work can differ by class; this is already a more generous scalar comparator than a single average client.
2. **Independent classes:** retain the scalar inequality and impose A^(k)x≤C separately for each client class, where A^(k) zeros all columns belonging to other classes. Each class receives its own copy of every machine capacity. It retains within-client route qualification and serial requirements but omits competition between classes. This is an optimistic independent-client relaxation, not a physical allocation of the shared pool.
3. **Architecture:** Ax≤C, with one shared capacity per machine. This implies both previous constraint sets because all entries are nonnegative.

Thus F_architecture⊆F_independent⊆F_scalar, for every capacity perturbation. This nesting is in feasible sets, not in a claim of nested statistical fits. Equal physical C and equal scalar totals are used for all models; the copied budgets of the independent relaxation are explicit, never called measured extra capacity.

The benchmark lexicographically maximizes (i) essential productive flux, (ii) all productive flux, (iii) all handled flux, then (iv) negative ATP expenditure. Earlier attained optima are fixed as lower bounds in subsequent LPs. This avoids arbitrary utility weights and prevents cheaper disposal from masquerading as essential rescue. With multiple essential classes, first maximize a common quota fraction (or impose specified q) rather than sum their fluxes; the supplied benchmark has one aggregate essential class. Pooling those 13 proteins into one class can conceal within-class shortages, so its essential fraction is only a class-level proxy. Primary optima inherit nesting; secondary optima need not when primary optima differ. Tests and ensemble checks respect this distinction.

At finite dose, positive excess loss can arise from scalar thresholding alone. It is not sufficient evidence for serial cooperation, proteotoxic collapse, a dynamical bifurcation, growth inhibition, or death. Those require additional state equations and independent measurements.
