# WO-04 — dynamics, stability, bifurcation

## Verdict

**VERDICT: PASS** — G4.1–G4.5 are met. Evidence was checked against the
independent source files and `wo04_results.json`, and the deterministic test
suite was rerun (10 tests, all passed).

## G4.1 — legacy scalar identities

`run_wo04.py:g41_scalar` reconstructs the cubic and equilibrium polynomial
with SymPy independently of the legacy audit. With `rho=4`, `chi=0.15`, the
cubic coefficients are `[-0.3, 0.4, 1.7, 5.0]`; the equilibrium polynomial
coefficients are `[-chi, 1-chi, -lambda+rho+1, -lambda]`. The computed
positive fold is `x=3.8907582151173496`, `lambda_fold=4.80218919587308`.
At `lambda=2`, all three polynomial roots have residual at most
`1.34e-15`. `g''=-2 chi-2 rho/(x+1)^3` is strictly negative for `x>-1`,
so the positive stationary point is the unique interior maximum; the fold is
non-degenerate (`g''=-0.3683850876`).

## G4.2 — operational collapse classification

At the published anchoring, operational collapse is the imposed aggregation
gate: `A=0.2501906286` at `P=0.0274023693`, matching the `A_max=0.25` gate
within the operational solver tolerance. The mathematical `J(P)` maximum is
at `P=0.0608124838`, with `A=0.6514810` and `J_fold/J_operational=1.19719`;
it occurs after the imposed gate. Removing the donorless `Phi` makes the
tested legacy clearance curve monotone for each of the six parameter
anchorings. This distinguishes the operational gate from a true fold.

## G4.3 — conservative-model steady-state counts

The analytic reduction in `bifurcation.py:analytic_v0` parameterizes the
steady-state curve by free non-native protein `U`. It proves the remaining
balances hold on that curve and writes the scalar balance derivative as a sum
of negative terms. Its boundary value is
`G(0)=s_P*(eps*mu+k_mis)/(k_mis+mu)>0` on the declared positive-`eps`
domain, while `G(U_max)<0`; strict decrease therefore gives exactly one
physical V0 steady state. V0 is also checked to reduce to WO-02 and conserve
both totals.

The deterministic parameter domain is declared in `bifurcation.py:sample`:
log-uniform `eps=[1e-4,0.5]`, `mu=[1e-5,1e-3]`, `P_T=[300,5000]`,
`C_T=[1,100]`, `k_on,k_off=[1e-2,10]`, `k_cat=[1e-3,1]`,
`k_d=[1e-5,1e-2]`, `k_a=[1e-6,1e-1]`, `k_dis=[1e-5,1e-2]`,
`k_dA=[1e-7,1e-3]`, `k_mis=[1e-7,1e-4]`; `phi` is uniform `[0.05,1]`.
V1 additionally samples `k_onA,k_offA=[1e-3,10]`; V2 samples
`k_dcat=[1e-4,1]` and fixes `k_dis=0`. With seed 100-series and 1,000
sets per variant, root counts are V0 `{1:1000}`, V1 `{1:946,3:54}`, V2
`{1:962,3:38}`. Maximum normalized full-system residuals are respectively
`4.06e-13`, `3.54e-13`, and `2.42e-13`; there are no zero-root cases or bad
residuals. A separate full six-dimensional Newton crosscheck over 150 V0
parameter sets found 467 physical converged solutions and no disagreement
with the reduced unique root.

## G4.4 — full-system stability

Every scanned root is reconstructed in the full six-state system and
classified from eigenvalues of its full Jacobian (`classify`), not from the
scalar reduction. There were zero single-root instability cases and zero
violations of the observed stable/unstable/stable ordering on multiroot
curves. The independent integration check perturbs both sides of each
representative root for three V1/V2 examples; all reported endpoint
classifications agree with the eigenvalue classifications. Deterministic
test `test_g44_pattern` also reproduces V1's three-root
`[stable, unstable, stable]` pattern.

## G4.5 — continuation and fold checks

Continuation of ten stored multistable examples each for V1 and V2 finds two
physical folds per example, two determinant sign changes, and three
equilibrium-curve crossings at the sampled epsilon, matching three roots.
At each located fold, `G/s_P` is within `1e-10` of zero, the smallest
absolute real Jacobian eigenvalue is below `1e-4` of the nearby reference,
and the determinant changes sign across the fold. The independent epsilon
sweep/fsolve locations agree with continuation (`G4.5_methods_agree=true`).
Ten V0 controls have zero folds and determinant sign changes and strictly
increasing epsilon along the curve. Thus the V0 mechanism missing a fold is
the nonlinear competition variant; V1/V2 restore saddle-nodes.

## Adversarial self-review

The 1,000-set scan is numerical evidence only on the explicitly sampled
domain; it does not establish global counts for V1/V2. V0's global uniqueness
claim is instead supported by the analytic monotonicity argument. Root scans
can miss tangent roots, but the reduction argument excludes them for V0, and
V1/V2 candidate folds are independently checked by continuation, `det(J)`
sign changes, and a vanishing eigenvalue. Full-system integration is a
crosscheck, not the stability proof. The legacy operational classification
is evaluated at the stated parameter anchors and does not establish a
biological collapse mechanism. None of these limits contradict a fixed gate.

## Test evidence

Command: `PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s WO-04 -p
'test_wo04.py' -v` — 10 tests passed in 49.319 s.

The supplied `run_1000.log` ends with `exit=0`; the independently parsed
`wo04_results.json` has `pass.G4.1` through `pass.G4.5` all `true` and
reports `n=1000` for each scan.
