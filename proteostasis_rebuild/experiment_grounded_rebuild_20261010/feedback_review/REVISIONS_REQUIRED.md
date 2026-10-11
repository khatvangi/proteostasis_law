# Revisions required: feedback/MINIMAL_DYNAMICS.md (and related files)

Independent review, 2026-10-10. Files under review were **not edited**.
Every number below comes from `check_feedback_claims.py` (output saved in `check_output.txt`,
fixed seeds, ~8 min on boron). All parameters are dimensionless placeholders; nothing is fitted
to E. coli.

## Headline

**§3 is wrong in the direction that matters.** `G < 1` is not "necessary but not sufficient". It is
sufficient, and with the boundary it is exact: the delayed loop is stable at every delay
**if and only if `G ≤ 1`**. The doc's exception ("a crossing can still exist if `e² > a²+μ²`") can never
happen: `G < 1` forces `e < μ`. 0 cases in 200,000 random draws; a two-line proof is below.

## 1. Delay-independent stability, derived

Characteristic function (checked against `findLoopTerms`: `u' = −au − bx`,
`x' = c·u(t−τ) − e·x(t−τ) − μx`):

```
f(s) = P(s) + Q(s)·e^(−sτ),   P = (s+a)(s+μ),   Q = e(s+a) + bc      (CORRECT in the doc)
```

- At τ = 0 every root is in the left half-plane (trace < 0, det = aμ(1+G) > 0).
- deg P = 2 > deg Q = 1 (retarded type), so roots can only enter the right half-plane by crossing
  the imaginary axis. `s = 0` is never a root: f(0) = aμ(1+G) > 0.
- `s = iω` requires |P(iω)|² = |Q(iω)|², i.e. `z² + A2·z + A0 = 0`, z = ω²,
  `A2 = a² + μ² − e²`, `A0 = (aμ)² − (ea+bc)² = (aμ)²(1 − G²)`. (CORRECT in the doc.)

**Case G < 1.** ea + bc < aμ, and since bc > 0, ea < aμ, so **e < μ**. Hence A2 = a² + (μ² − e²) > 0 and
A0 > 0. Both coefficients positive → no positive root z → no crossing at any τ → stable for every τ ≥ 0.

**Case G = 1.** A0 = 0, A2 > 0: roots z = 0 and z = −A2 < 0. z = 0 would mean s = 0, which is not a
root. So stable for every finite τ, but with vanishing margin: rightmost real part −0.65, −5.5e−3,
−8.3e−6 at τ = 1, 10, 100 (a = μ = 1, e = 0.5, bc = 0.5).

**Case G > 1.** A0 < 0 → exactly one positive root z₊ (139,707/139,707 draws), and F′(z₊) > 0, so by
the standard crossing-direction rule every crossing at ω₊ = √z₊ is destabilising: stable iff τ < τ*.
Verified numerically only *near* τ*: in 40/40 loops a root sits on the imaginary axis at the code's τ*,
the loop is stable at 0.97τ* and unstable at 1.03τ* (spectral discretisation + Newton polish).
The `τ*` formula `[−arg(−P/Q) mod 2π]/ω` and `findCriticalDelay` are CORRECT.

Numerical support for G < 1: 40 random loops × τ ∈ {0, 0.1, 1, 10, 50}; largest rightmost real part
−5.7e−4 (always negative).

**Plain meaning.** Since steady-state sensitivity is divided by (1+G), the result says: a loop that
cuts the load-induced loss of spare capacity by at most half is stable at any delay; a loop that cuts it
by more than half has a finite delay limit. This is sharper than what the doc claims.

### Counterexamples to wrong inequalities
| Candidate condition | Verdict | Witness |
|---|---|---|
| "G < 1 but e² > a²+μ² can cross" (doc §3) | impossible | proof above; 0 hits in 2×10⁵ draws |
| e < μ suffices | false | a=b=1, c=2, e=0.5, μ=1: G = 2.5, τ* = 1.267 (Re of rightmost root −0.040 at 0.9τ*, +0.030 at 1.1τ*) |
| e² < a²+μ² suffices | false | same witness |
| e² > a²+μ² is a separate failure mode | not separate | it forces e > μ, hence G > 1 (e.g. a=μ=1, e=2: G = 2.0001, τ* = 1.209) |

## 2. Claim-by-claim verdicts

| # | Claim (file §) | Verdict | Replacement / note |
|---|---|---|---|
| 1 | Characteristic equation and z-polynomial (MD §3) | PASS | — |
| 2 | "If G > 1 … a finite critical delay τ* beyond which the loop oscillates" (MD §3) | PASS (linear) | Add: "Linear result only. Whether the nonlinear loop settles into a small stable oscillation or a large one (super- vs subcritical Hopf) was not checked." |
| 3 | "G < 1 is necessary but not sufficient for stability at all delays (a crossing can still exist if e² > a²+μ²)" (MD §3) | **FAIL** | "**The linearised loop is stable at every delay if and only if G ≤ 1.** G < 1 forces e < μ, so both coefficients of the crossing polynomial are positive and no crossing exists; at G = 1 the only candidate is ω = 0, which is not a root. For G > 1 there is exactly one crossing frequency and the loop is stable exactly for τ < τ*." |
| 4 | "tight homeostasis needs G well above 1, and G > 1 guarantees a delay limit" (MD §3) | PASS, can be sharpened | "Feedback that more than halves the steady-state loss of spare capacity (1+G > 2) necessarily has a finite delay limit; feedback that halves it or less is stable at any delay." |
| 5 | test `test_weak_loop_without_crossing…` comment "G < 1 and a²+μ² > e²" | PASS, redundant | second clause is implied by the first |
| 6 | J, trace, det, stable at τ = 0 (MD §1) | PASS | — |
| 7 | det = aμ(1+G); dF*/dL = −[CK/(K+U)²]/[a(1+G)] (MD §2) | PASS | derived independently: the c- and e-terms cancel exactly in dF*/dL. Note a itself depends on the equilibrium C, so "open vs closed loop" compares different a. |
| 8 | spare capacity = k·F (MD intro) | PASS, but by construction | it follows from the ASSUMED saturating-occupancy form F = CK/(K+U), not from data. Say so. |
| 9 | "The pool can never exceed Cmax" / "At any state, served flux ≤ k·Cmax" / "for every gain, Hill coefficient and delay" (MD §4) | **FAIL as stated** | see §3 below |
| 10 | Test of §4 (`test_load_above_max_service…`) | weaker than the claim | checks equilibria only, not trajectories. Doc should say "Tested at equilibria in 100 draws." |
| 11 | Box invariance, uniqueness, global convergence at τ = 0 (POSITIVE_INVARIANCE.md) | PASS | conditions (constant μ, fixed mode, β0 > 0) are already stated there. Note MD §4 contradicts its own Consequence 1 (trajectories *starting outside* the band enter it, i.e. C can start above Cmax). |
| 12 | Growth-coupled ceiling P(φ0+φm) "for any μ" (MD §5) | PASS with condition | holds for constant P and C(0) ≤ P(φ0+φm). |
| 13 | Collapsed state under growth arrest (MD §5) | PASS mathematically, physical caveat **required** | see §4 below |
| 14 | L = served + δU* + μU*; raising δ lowers served flux (MD §7) | PASS, and it is a theorem | Proof (constant μ, fixed mode): if served rose, (μ+δ)U would fall, so U falls; served = kCU/(K+U) rising with U falling needs C up; C up at steady state needs h up, i.e. F down; but F = CK/(K+U) goes up. Contradiction. Numerically 0/100 exceptions. Replace "(tested)" with "(always, at constant growth; proof in feedback_review)". |
| 15 | `findSigmaActivity` returns h′ = 0 at F = 0 (code) | **FAIL (minor)** | true h′(0) = −1/Ks for n = 1, 0 for n > 1, −∞ for n < 1. Harmless for every result shown, since F > 0 whenever C > 0 (β0 > 0). Fix: `dh = -n * F**(n-1) / Ks**n / (1+r)**2` for n ≥ 1, and forbid n < 1. |

## 3. "C never above Cmax": what is actually true

`dC/dt ≤ β0 + βm − μC` gives `C(t) − Cmax ≤ (C(0) − Cmax)·e^(−μt)`. So the band is **attracting and
forward-invariant**, not an absolute bound. It needs (i) C(0) ≤ Cmax and (ii) constant μ.
Counterexamples (L = 5, k = 2, K = 1, μ = 0.5, β0 = 0.05, βm = 1, Ks = 1, n = 2; Cmax = 2.1, kCmax = 4.2):

- **Initial condition.** C(0) = 3Cmax: C stays above Cmax until t = 3.27; served flux reaches 6.3,
  above kCmax = 4.2 and above the influx itself (1.26 L). The "excess load is never served" claim fails
  transiently.
- **Growth shift.** Steady state at μ = 0.25 has C = 2.844; when μ jumps to 0.5 the new ceiling is 2.1
  and C stays above it until t = 1.73. Any growth upshift (or any pre-induced history) puts the cell
  above its new ceiling.
- **Growth arrest, fixed mode.** Cmax = (β0+βm)/μ → ∞ as μ → 0 (2.1, 21, 210 for μ = 0.5, 0.05, 0.005).
  There is no ceiling in fixed mode without growth; the code already refuses `Ktox` in fixed mode for
  this reason.

Replacement for MD §4 first paragraph:
> "If the pool starts at or below Cmax = (β0+βm)/μ and growth is constant, it stays there, because
> h ≤ 1; a pool that starts above Cmax (pre-induced cells, or a growth upshift that lowers Cmax) relaxes
> toward it at rate μ. Once C ≤ Cmax, served flux ≤ k·Cmax, so if L > k·Cmax the unserved flux is at
> least L − k·Cmax for every gain, Hill coefficient and delay. Tested at equilibria in 100 draws; the
> transient exception is shown in feedback_review/."

## 4. Growth coupling: the collapsed state depends on a modelling choice

The §5 bistability exists only in "growth" mode, where synthesis = μ(U)·P·(φ0+φm·h). That choice makes
**chaperone synthesis go to zero when growth stops**, and makes the time to refill the pool scale with
1/μ. The opposite choice (fixed mode) has no ceiling at all when growth stops. Real cells are neither:
heat-shock protein synthesis is usually described as continuing or rising while growth pauses, though no
measurement of that is retrieved in this folder. So the collapsed state is a property of the assumed
coupling, not a finding.

Replacement sentence for MD §5, second bullet, appended:
> "This collapse depends on assuming chaperone synthesis is proportional to growth, so synthesis stops
> when growth stops. With growth-independent synthesis the ceiling disappears instead. Which is closer
> to E. coli under load is not measured here; the collapsed state is conditional on that assumption."

Also: the linear loop terms (a, b, c, e) in §1 to §3 are derived for fixed production and constant μ.
In §5, "a turns negative" uses the same symbol for a different Jacobian entry; say "the U-diagonal of the
growth-mode Jacobian turns negative".

## 5. Missing compartments (excess workload, native product)

- No native-product or aggregate state. "Served" is chaperone turnover, assumed to equal successful
  folding with one pass per client. The excess load in §4 can only be diluted or degraded; there is no
  aggregation term (B_agg) and no misfolded-then-rebound cycle. So "unserved flux ≥ L − kCmax" is a
  statement about chaperone throughput, not about native yield. MD §7 already says productive folding
  needs measured fates; §4 should say the same.
- F = CK/(K+U) treats U as not depleted by binding (Michaelis–Menten occupancy). It is ASSUMED.

## 6. Assumption vs experimental finding

| Element | Status |
|---|---|
| Free chaperone represses σ32 (DnaK/J, GroEL/S); FtsH degrades σ32; regulon includes chaperones and proteases; TF outside the regulon | experimental (LIT, per EXPERIMENTAL_CONTROL.md; not re-retrieved in this review) |
| Loop is sensitive to <1.5-fold chaperone change (Tomoyasu 1998) | experimental, **qualitative**. It is not evidence that G > 1 in the sense above. MD §3 "the real loop is very sensitive (high gain…)" should read: "the real loop is sensitive to small perturbations (Tomoyasu 1998), but G has not been measured, so whether it exceeds 1 (and hence whether any delay limit exists) is unknown." |
| Single lumped chaperone, single delay τ acting on both client and chaperone arms of F | ASSUMED |
| Hill repression h(F), saturating service kCU/(K+U), free pool CK/(K+U) | ASSUMED |
| Constant μ; or synthesis ∝ μ (growth mode) | ASSUMED, and the §5 conclusion depends on it |
| One pass of service = native protein; no aggregate pool | ASSUMED |
| G, τ, k, Cmax, β's | UNMEASURED |
| G ≤ 1 ⇔ delay-independent stability; G > 1 ⇒ finite τ* | mathematical result, linear, local to an equilibrium |

## Not checked (time bound)

Nonlinear behaviour beyond τ* (Hopf direction, amplitude); distributed rather than discrete delay;
growth-mode Jacobian with delay; behaviour beyond 1.03τ* for G > 1 (theory says it stays unstable;
numerics checked only near τ*).
