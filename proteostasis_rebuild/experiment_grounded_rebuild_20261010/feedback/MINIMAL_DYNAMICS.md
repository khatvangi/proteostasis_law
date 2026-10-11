# Minimal dynamics of the σ32 capacity loop

All results here are algebra on a structural model. **No parameter is calibrated.** The code
(`feedback.py`) and tests (`test_feedback.py`, 12 tests, all pass) use dimensionless placeholder
values only to check signs, identities and thresholds; none of those numbers is a statement about
E. coli.

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

One useful identity: **spare capacity = k·F.** The unused turnover of the pool is exactly `k` times
the free chaperone, which is the quantity σ32 senses. The loop therefore senses spare capacity
directly, not load.

## 1. Without delay every balance point is stable

Linearising at an equilibrium:

```
J = [ −a   −b      ]    a = kCK/(K+U)² + μ + δ   (>0)   client removal slope
    [  c   −(μ+e)  ]    b = kU/(K+U)             (>0)   service gain per chaperone
                        c = βm|h′|·CK/(K+U)²     (>0)   more client → less free → more synthesis
                        e = βm|h′|·K/(K+U)       (>0)   more chaperone → more free → less synthesis
```

trace = −(a+μ+e) < 0, det = a(μ+e) + bc > 0, so every equilibrium is a stable node or focus for
any positive gains. Combined with the invariant box and negative divergence (POSITIVE_INVARIANCE.md) the equilibrium is
unique and globally attracting when τ = 0. Tested: analytic J equals finite-difference J and has this sign pattern in 200/200 random
placeholder draws.

## 2. What the feedback buys: error divided by (1 + loop gain)

Define the zero-frequency loop gain `G = (e·a + b·c)/(a·μ)`. Then det = aμ(1+G) exactly, and

```
dF*/dL = −[CK/(K+U)²] / [a·(1+G)]        (open loop: G = 0)
```

so a rise in load always lowers spare capacity (dF*/dL < 0 — proportional control, no perfect
adaptation), and σ32 feedback divides that loss by `1+G`. Tested against finite differences.
"Feedback restores spare capacity" is therefore *partial* restoration; it is complete only in the
limit G → ∞.

## 3. Delay: high gain plus delay oscillates

With the σ32 arm delayed by τ the characteristic equation is

```
(λ+a)(λ+μ) + (eλ + ea + bc)·e^(−λτ) = 0
```

A root reaches the imaginary axis at λ = iω where z = ω² solves
`z² + (a²+μ²−e²)z + a²μ² − (ea+bc)² = 0`, giving
`τ* = [−arg(−P(iω)/Q(iω)) mod 2π] / ω`.

- If **G > 1** the constant term is negative, so a positive root always exists: there is a finite
  critical delay τ* beyond which the loop oscillates (Hopf). Tested in all G>1 draws, and by direct
  integration of the delayed linear system (decays at 0.8τ*, grows at 1.25τ*).
- G < 1 is necessary but not sufficient for stability at all delays (a crossing can still exist if
  e² > a²+μ²).

Section 2 and 3 together are the trade-off: **tight homeostasis needs G well above 1, and G > 1
guarantees a delay limit.** Tomoyasu 1998 shows the real loop is very sensitive (high gain, LIT,
qualitative), and the response is observed to be transient (Straus 1987). Neither G nor τ is measured,
so whether E. coli sits near τ* is unknown.

## 4. Finite capacity: feedback cannot serve load above maximum service

The pool can never exceed `Cmax = (β0+βm)/μ` (h ≤ 1). At any state, served flux ≤ k·Cmax, so

```
if L > k·Cmax:   unserved flux = L − served ≥ L − k·Cmax > 0
```

for every gain, Hill coefficient and delay. The excess must sit as unfolded client (diluted only by
growth) or be degraded. Tested in 100 draws with L up to 10× k·Cmax. This is the dynamic counterpart
of the static inequality `B ≤ C`: feedback can move the pool *within* [β0/μ, Cmax]; it cannot move the
ceiling. Negative feedback fails in exactly two ways here: the ceiling is too low (βm/μ small) or the
delay is too long (τ > τ*).

## 5. Growth, dilution and energy

- **Ceiling is an allocation, not a rate.** If synthesis scales with growth (production =
  μ·P·(φ0+φm·h), Scott 2010 type allocation), then `Cmax = P·(φ0+φm)` for any μ: the cap is a maximal
  proteome fraction. Tested for three μ values. Energy/ATP enters only through `k` and through what
  else that proteome fraction would have built; neither is quantified here.
- **Growth arrest can create a collapsed state.** If load slows growth (μ(U) decreasing), dilution of
  unfolded client falls as load rises. Where d(μ(U)·U)/dU is negative enough, `a` turns negative,
  det < 0 and a saddle appears: one served state and one collapsed state with σ32 fully on and C
  already at its ceiling. Tested with a placeholder case (3 equilibria, middle one a saddle); a small
  further load increase removes the served state. This is the known growth-dilution bistability
  (Klumpp 2009; Tan 2009), not new. **Local stability of the served state therefore does not imply
  survival**: a large enough transient crosses into the collapsed basin.

## 6. Stoichiometry: one input cannot set two capacities

DnaK/J and GroEL/S both repress σ32 (LIT), but σ32 is one input. If the model is split into two
machines whose synthesis is driven by the same σ32 activity (a dimension argument, not coded here),
the reachable set of steady pools {C_DnaK, C_GroEL} as σ32 activity varies is a one-parameter curve
(Zhang 2014's ~6× vs ~4× is one point on it). Two consequences:

- A client mix that needs a different GroEL:DnaK ratio than the curve provides can be *admissible*
  (enough total capacity exists for some allocation) yet *not reachable* by σ32 control.
- Because σ32 senses a combined free-chaperone signal, free DnaK can keep σ32 repressed while GroEL
  is saturated. A single sensor channel cannot guarantee spare capacity in each machine separately.

Trigger factor is outside the regulon (LIT), so TF-dependent routes get no feedback at all.

## 7. Admissible allocation ≠ implemented allocation ≠ productive allocation

Three different statements, often conflated:

1. **Admissible**: there exists a routing of client flux through machine capacities that serves the
   load (the static / Hall-type condition in the route model). Pure existence; says nothing about how.
2. **Implemented**: the σ32 loop actually drives the pools to such an allocation, stably. Requires
   L < k·Cmax, τ < τ* (if G > 1), the right stoichiometry (§6), and being in the served basin (§5).
3. **Productive**: the balance is closed by folding, not by destruction. The σ32 regulon also raises
   proteases (ClpX/P/A ~2× in Zhang 2014; FtsH is σ32-induced, Herman 1995). Exact bookkeeping:
   `L = served + δU* + μU*` (tested to 1e-9). Raising δ lowers U* and *lowers chaperone-served flux*
   while the books still balance (tested). A cell can satisfy the capacity balance by degrading
   clients it needed — balanced, but harmful.

The static inequality speaks only to (1). Feedback is the experimentally established mechanism for
(2). Neither guarantees (3); that needs a measured folded-vs-degraded fate per client class.
