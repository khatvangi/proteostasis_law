# Finite-Resource Dynamical System

## State equations

Let `U` be damaged or misfolded soluble monomer and `A` aggregate burden, in compatible concentration units. One baseline model is

`dU/dt = J - v_ref(U,C_f) - v_degU(U,D_f) - n(U) - g(U,A) + v_dis(A,C_f)`

`dA/dt = n(U) + g(U,A) - v_dis(A,C_f) - v_degA(A,D_f)`.

Possible rate laws are

- refolding: `v_ref = k_ref C_f U/(K_ref+U)`;
- soluble degradation: `v_degU = k_U D_f U/(K_U+U)`;
- nucleation: `n = k_n U^m`, `m>1`;
- growth: `g = k_g U A`;
- disaggregation: `v_dis = k_dis C_f A/(K_dis+A)`;
- aggregate clearance: `v_degA = k_A D_f A/(K_A+A)`.

Other dimensionally valid kinetics may replace these forms. Disaggregation returns material to `U` in this baseline; if it directly yields foldable/native product, the stoichiometry must be changed in both equations.

## Conserved resource pools

Chaperone occupancy includes normal nascent chains `N`, soluble damaged substrate, aggregates, and optionally other clients:

`C_tot = C_f + C_N + C_U + C_A + C_O`,

`C_N = N_f C_f/K_N`, `C_U = U_f C_f/K_CU`, `C_A = A_f C_f/K_CA` under rapid-equilibrium binding.

If the listed substrate variables denote **free binding-competent** concentrations, then

`C_f = C_tot / (1 + N_f/K_N + U_f/K_CU + A_f/K_CA + O_f/K_CO)`.

If `U` or `A` is a total pool, free substrate and complexes must instead be solved jointly from the conservation equations. It is incorrect to substitute total substrate directly into a free-resource formula.

Similarly,

`D_tot = D_f + D_U + D_A + D_O`,

with rapid-equilibrium relations only in terms of free substrates. Degradation capacity can also require catalytic queue states; conservation must then include all occupied enzyme complexes.

Ordinary nascent-chain folding is not itself damage influx, but `C_N` consumes capacity and can alter the stability boundary.

## Stability criterion

At an equilibrium `x*`, local asymptotic stability requires all eigenvalues of the Jacobian `partial F/partial x` to have negative real parts. Viability additionally requires `x*` (or another bounded attractor) to remain below damage thresholds and the relevant initial state to lie in its basin. Existence alone is insufficient.

## Illustrative one-variable reduction

For intuition only, collapse aggregation and resources into

`dM/dt = j - a M/(K+M) + b M^2`, with `j,a,K,b>0`.

The vector field is convex because its second derivative is positive. Depending on parameters it has no positive fixed point, one double root at a fold, or two positive fixed points. With two roots, the lower is stable and the upper is unstable; the upper root is a threshold beyond which burden grows in this reduced model. At the fold the stable and unstable roots annihilate. This model does **not** contain a second stable high-burden attractor and does not imply hysteresis. Such phenomena require additional state structure and separate demonstration.

