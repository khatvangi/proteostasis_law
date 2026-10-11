# Minimal dynamics of the σ32 capacity loop

All results here are algebra on a structural model. **No parameter is calibrated.** The code
(`feedback.py`) and tests (`test_feedback.py`, 21 tests, all pass) use dimensionless placeholder
values only to check signs, identities and thresholds; none of those numbers is a statement about
E. coli. Corrections made after the independent review are listed in `CORRECTIONS.md`.

## Model (ASSUMED structure; each arrow is LIT, see EXPERIMENTAL_CONTROL.md)

```
dU/dt = L − k·C·U/(K+U) − (μ+δ)·U            client: influx, chaperone service, dilution, proteolysis
dC/dt = β0 + βm·h(F(t−τ)) − μ·C               chaperone: basal + σ32-driven synthesis, dilution
F     = C·K/(K+U)                             free chaperone (saturating occupancy)
h(F)  decreasing, 0<h<1                       σ32 activity repressed by free chaperone
```

`U` unfolded client needing the chaperone, `C` total σ32-regulated chaperone (one lumped machine),
`L` influx, `k` maximum turnover per chaperone, `K` half-saturating load, `μ` growth (= dilution),
`δ` proteolytic removal, `τ` delay from σ32 activity to working chaperone. σ32 itself is collapsed
into `h` plus the delay; the separate FtsH stability loop is folded into `h` (both lower σ32 output
when free chaperone is high).

What is assumed and what is observed:

| Element | Status |
|---|---|
| Free chaperone represses σ32; unfolded client titrates it; σ32 induces chaperones and proteases; TF outside the regulon | LIT (EXPERIMENTAL_CONTROL.md) |
| One lumped chaperone; one discrete delay τ acting on the whole free-pool signal | ASSUMED |
| Hill repression h(F); saturating service k·C·U/(K+U); free pool C·K/(K+U) with U not depleted by binding | ASSUMED |
| Constant μ (§1–§4, §7), or synthesis ∝ μ (§5 growth mode) | ASSUMED; the §5 conclusion depends on it |
| "Served" = chaperone turnover. No native-product state and no aggregate state | ASSUMED; served flux is **not** native protein output |
| G, τ, k, Cmax, β's | UNMEASURED |

One identity: **spare capacity = k·F.** The unused turnover of the pool, k·C − k·C·U/(K+U), equals
`k` times the free chaperone, which is the quantity σ32 senses. This holds **by construction** of
the assumed saturating-occupancy form, not because of data. Within that form, the loop senses spare
capacity directly, not load.

## 1. Without delay every balance point is stable

Linearising at an equilibrium (fixed production, constant μ):

```
J = [ −a   −b      ]    a = kCK/(K+U)² + μ + δ   (>0)   client removal slope
    [  c   −(μ+e)  ]    b = kU/(K+U)             (>0)   service gain per chaperone
                        c = βm|h′|·CK/(K+U)²     (>0)   more client → less free → more synthesis
                        e = βm|h′|·K/(K+U)       (>0)   more chaperone → more free → less synthesis
```

trace = −(a+μ+e) < 0, det = a(μ+e) + bc > 0, so every equilibrium is a stable node or focus for
any positive gains. With the invariant box and negative divergence (POSITIVE_INVARIANCE.md) the
equilibrium is unique and globally attracting **for τ = 0, fixed production and constant μ**.
Tested: analytic J equals finite-difference J and has this sign pattern in 200/200 random placeholder
draws.

## 2. What the feedback buys: error divided by (1 + loop gain)

Define the zero-frequency loop gain `G = (e·a + b·c)/(a·μ)`. Then det = aμ(1+G) exactly, and

```
dF*/dL = −[CK/(K+U)²] / [a·(1+G)]        (open loop: G = 0)
```

(the c- and e-terms cancel exactly in this derivative). A rise in load always lowers spare capacity
(dF*/dL < 0: proportional control, no perfect adaptation), and σ32 feedback divides that loss by
`1+G`. Tested against finite differences. "Feedback restores spare capacity" is therefore *partial*
restoration; it is complete only in the limit G → ∞. Caveat: `a` depends on the equilibrium C, so
comparing "open" with "closed" loop compares two different values of `a`.

## 3. Delay: stable at every delay if and only if G ≤ 1

This section is a **linear, local** result: it describes small perturbations about one equilibrium
of the fixed-production, constant-μ model. By the linearised-stability principle for retarded delay
equations it gives local asymptotic stability or instability of the nonlinear equilibrium. It says
nothing global for τ > 0.

With the σ32 arm delayed by τ the characteristic function is

```
f(s) = P(s) + Q(s)·e^(−sτ),    P(s) = (s+a)(s+μ),    Q(s) = e·s + e·a + b·c
```

**Theorem.** For a, b, c, e, μ > 0:
(i) if G ≤ 1, every root of f has negative real part for every τ ≥ 0;
(ii) if G > 1, there is a finite τ* > 0 such that the equilibrium is stable for τ < τ* and unstable
for every τ > τ*.

**Proof.**

1. *τ = 0.* f(s) = s² + (a+μ+e)s + aμ(1+G): both coefficients positive, so both roots are in the left
   half-plane.
2. *Roots move continuously in τ and enter the right half-plane only across the imaginary axis.*
   deg P = 2 > deg Q = 1 (retarded type). Roots that appear as τ leaves 0 come from Re s = −∞. This
   is the standard root-continuity property of retarded quasipolynomials.
3. *s = 0 is never a root.* f(0) = aμ + ea + bc = aμ(1+G) > 0.
4. *Crossing condition.* s = iω (ω > 0) requires |P(iω)|² = |Q(iω)|². With z = ω² this is

   ```
   F(z) = z² + A2·z + A0 = 0,   A2 = a² + μ² − e²,   A0 = (aμ)² − (ea+bc)² = (aμ)²(1 − G²)
   ```

5. *Case G < 1.* ea + bc < aμ and bc > 0, so ea < aμ, so **e < μ**. Then A2 = a² + (μ² − e²) > 0 and
   A0 > 0. A quadratic with positive coefficients has no positive root, so there is no crossing at any
   τ. With step 1 and step 2, this proves (i) for G < 1.
6. *Case G = 1 (equality).* A0 = 0. Again ea < aμ gives e < μ, so A2 > 0. The roots are z = 0 and
   z = −A2 < 0. z = 0 means s = 0, which step 3 excludes. So there is no crossing, and the loop is
   stable for every finite τ. The margin is not uniform: it vanishes as τ → ∞. The review measured
   rightmost real parts −0.65, −5.5×10⁻³ and −8.3×10⁻⁶ at τ = 1, 10 and 100 (a = μ = 1, e = 0.5,
   bc = 0.5; `../feedback_review/check_output.txt`). That proves (i).
7. *Case G > 1.* A0 < 0, so F has exactly one positive root z₊, and F′(z₊) = 2z₊ + A2 = √(A2² − 4A0) > 0.
   By the standard crossing-direction rule, sign Re(ds/dτ) = sign F′(z₊). So every crossing at
   ω₊ = √z₊ moves a pair of roots from left to right, and none ever returns. The crossings occur at
   τ_j = τ* + 2πj/ω₊, with

   ```
   τ* = [−arg(−P(iω₊)/Q(iω₊)) mod 2π] / ω₊   (> 0 because no root lies on the axis at τ = 0)
   ```

   So the loop is stable exactly for τ < τ*, has an imaginary pair at τ*, and is unstable for all
   τ > τ*. That proves (ii). ∎

**What the theorem rules out.** The earlier text said "G < 1 is necessary but not sufficient (a
crossing can still exist if e² > a²+μ²)". That was **wrong**. Step 5 shows G < 1 forces e < μ, so
e² > a²+μ² can never occur together with G < 1. Neither "e < μ" nor "e² < a²+μ²" is sufficient on its
own. The witness is a = b = 1, c = 2, e = 0.5, μ = 1: both conditions hold, yet G = 2.5 and
τ* = 1.267.

**Tests.**
- 5,000 random loops. Every G ≤ 1 draw has e < μ, A2 > 0 and no crossing. Every G > 1 draw has
  exactly one positive z, F′(z₊) > 0, and f(iω₊) = 0 at τ*.
- The equality case: no crossing.
- Direct integration of the delayed linear system from a **nonconstant** (oscillating) history:
  - G = 0.5 decays at τ = 0.5, 5 and 20;
  - G = 1 decays at τ = 1;
  - the G = 2.5 witness decays at 0.9τ* and grows at 1.1τ*, 2τ* and 3τ*.

  A nonconstant history is used so that stability is not checked from one special start.
- The independent review checked the root locations with a spectral method: stable at 0.97τ* and
  unstable at 1.03τ* in 40/40 loops.

**Not checked.** Beyond τ* the linear analysis only says the equilibrium is unstable. Whether the
nonlinear loop settles into a small stable oscillation or jumps to a large one (super- vs subcritical
Hopf) was not checked. Distributed delays and the growth-mode model with delay were not analysed.

**Plain meaning, with §2.** The steady-state loss of spare capacity is divided by 1+G. So feedback that
cuts that loss by half or less (G ≤ 1) is stable at any delay. Feedback that cuts it by more than half
(G > 1) necessarily has a finite delay limit. Tight homeostasis therefore comes with a delay limit.

**What is known for E. coli.** The real loop is sensitive to small perturbations: <1.5-fold changes in
DnaK/J matter (Tomoyasu 1998, LIT). The response is transient (Straus 1987, LIT). G has not been
measured, so whether it exceeds 1, and hence whether any delay limit exists, is unknown. τ is not
measured either.

## 4. Finite capacity: feedback cannot raise its own ceiling

Let Cmax = (β0+βm)/μ. Because h ≤ 1, dC/dt ≤ μ(Cmax − C), so

```
C(t) − Cmax ≤ (C(0) − Cmax)·e^(−μt)        (constant μ)
```

The band C ≤ Cmax is therefore **forward-invariant and attracting**. It is not an absolute bound:

- **Pool already below the ceiling.** If the pool starts at or below Cmax and growth is constant, it
  stays there. This holds for any delay too, since h ∈ (0,1) whatever the history.
- **Pool above the ceiling.** A pool that starts above Cmax relaxes toward it at least at rate μ. This
  covers pre-induced cells, and any change that lowers Cmax, such as a growth upshift.
  - *Example 1 (initial condition above the ceiling).* C(0) = 3Cmax, with L = 5, k = 2, μ = 0.5,
    β0 = 0.05, βm = 1, giving Cmax = 2.1. C stays above Cmax until t ≈ 3.3. Served flux reaches 6.3,
    which is above k·Cmax = 4.2 and above the influx L.
  - *Example 2 (growth upshift).* The steady state at μ = 0.25 has C above the new ceiling for
    μ = 0.5, and C stays above it for t > 1.
- **Growth arrest.** In fixed mode Cmax → ∞ as μ → 0, so there is no ceiling at all when growth
  stops. This is why the code refuses load-dependent growth in fixed mode.

Once C ≤ Cmax, served flux ≤ k·Cmax. So

```
if L > k·Cmax and C ≤ Cmax:   unserved flux = L − served ≥ L − k·Cmax > 0
```

for every gain, Hill coefficient and delay. The excess sits as unfolded client (diluted by growth) or
is degraded.

Tested at equilibria in 100 draws with L up to 10× k·Cmax. The trajectory exceptions above are tested
in `TestCeiling`.

"Served" here means chaperone throughput, not native yield. The model has no native-product or
aggregate pool (no B_agg term) and assumes one pass of service per client. So "unserved ≥ L − k·Cmax"
bounds chaperone throughput, not folded output.

Within its scope, this is the dynamic counterpart of the static inequality `B ≤ C`: feedback moves
the pool *within* [β0/μ, Cmax], but it cannot move the ceiling. Negative feedback fails in exactly two
ways here: the ceiling is too low (βm/μ small), or G > 1 and the delay is too long (τ > τ*).

## 5. Growth, dilution and energy

**The ceiling is an allocation, not a rate.** Suppose synthesis scales with growth: production =
μ·P·(φ0+φm·h), a Scott 2010-type allocation (ASSUMED). Then dC/dt ≤ μ(U)·(P(φ0+φm) − C), so
`Cmax = P·(φ0+φm)` for any μ, and the cap is a maximal proteome fraction.

- This holds for constant P and C(0) ≤ P(φ0+φm).
- A pool above the cap relaxes only at rate μ(U), which goes to zero if growth stops.
- Tested at equilibria for three μ values.
- Energy and ATP enter only through `k` and through what else that proteome fraction would have
  built. Neither is quantified here.

**Growth arrest can create a collapsed state, under a specific assumption.** If load slows growth
(μ(U) decreasing), dilution of unfolded client falls as load rises. Where d(μ(U)·U)/dU is negative
enough, the U-diagonal entry of the growth-mode Jacobian, ∂(dU/dt)/∂U, becomes positive. Then
det < 0 and a saddle appears.

- The saddle separates a served state from a collapsed state. In the collapsed state σ32 is fully on
  and C is already at its ceiling.
- Tested with a placeholder case: 3 equilibria, the middle one a saddle, with U-diagonal +0.0095.
  A small further load increase removes the served state.
- This is the known growth-dilution bistability (Klumpp 2009; Tan 2009), not new.

**This collapse depends on assuming chaperone synthesis is proportional to growth, so synthesis stops
when growth stops.** In the placeholder case, synthesis in the collapsed state is about 1000× below
the served state (tested). With growth-independent synthesis (fixed mode) the ceiling disappears
instead. Heat-shock protein synthesis is usually described as continuing or rising while growth
pauses, but no measurement of that is retrieved in this folder. So the collapsed state is
conditional on the assumption, not a finding about E. coli.

**Local stability of the served state does not imply survival.** Within the growth-mode model, a large
enough transient crosses into the collapsed basin.

## 6. Stoichiometry: one input cannot set two capacities

DnaK/J and GroEL/S both repress σ32 (LIT), but σ32 is one input. Split the model into two machines
whose synthesis is driven by the same σ32 activity (a dimension argument, not coded here). Then, as
σ32 activity varies, the reachable steady pools {C_DnaK, C_GroEL} trace a one-parameter curve.
Zhang 2014's ~6× vs ~4× is one point on that curve, in an engineered strain. Two consequences:

- A client mix can need a different GroEL:DnaK ratio than the curve provides. It is then
  *admissible* (enough total capacity exists for some allocation) yet *not reachable* by σ32 control.
- σ32 senses a combined free-chaperone signal, so free DnaK can keep σ32 repressed while GroEL is
  saturated. A single sensor channel cannot guarantee spare capacity in each machine separately.

Trigger factor is outside the regulon (LIT), so TF-dependent routes get no feedback at all.

## 7. Admissible allocation ≠ implemented allocation ≠ productive allocation

These are three different statements, often conflated:

1. **Admissible.** Some routing of client flux through machine capacities serves the load (the static,
   Hall-type condition in the route model). This is pure existence; it says nothing about how.
2. **Implemented.** The σ32 loop actually drives the pools to such an allocation, stably. Within this
   model that requires L < k·Cmax with C(0) ≤ Cmax (§4), τ < τ* if G > 1 (§3), the right
   stoichiometry (§6), and being in the served basin (§5).
3. **Productive.** The balance is closed by folding, not by destruction. The σ32 regulon also raises
   proteases: ClpX/P/A ~2× in Zhang 2014, and FtsH is σ32-induced (Herman 1995).

The exact bookkeeping at steady state is `L = served + δU* + μU*` (tested to 1e-9).

**Raising δ always lowers chaperone-served flux (constant μ, fixed mode).** Proof:
- Suppose δ rises and served flux does not fall.
- Then (μ+δ)U* = L − served does not rise. With δ larger, U* falls.
- Served = kCU/(K+U) not falling while U falls requires C to rise.
- At steady state, C = (β0+βm·h(F))/μ, so C rising requires h to rise, i.e. F to fall.
- But F = CK/(K+U) rises when C rises and U falls. Contradiction. ∎

Numerically, the review found 0 exceptions in 100 draws, and the test shows it on a ladder of δ.

So a cell can satisfy the capacity balance by degrading clients it needed: balanced, but harmful. Even
"served" is not "folded". The model has no native-product state, so served flux is an upper bound on
chaperone-assisted folding, not a measure of it.

The static inequality speaks only to (1). Feedback is the experimentally established mechanism for
(2), though in this model only under the conditions listed. Neither guarantees (3). That needs a
measured folded-vs-degraded fate per client class.
