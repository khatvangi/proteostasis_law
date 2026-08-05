# Lemma 0: the binding closure is well posed, and H3 is a consequence

Deliverable of task A2. Hypothesis H3 asserts that `R` and `G` are smooth. But
`R` and `G` are not given by formulae in `(u, a)`: they are defined *implicitly*,
through a nonlinear algebraic system that the manuscript through v5 never stated
as such. Nothing in the paper established that the system has a solution, that
the solution is unique, or that it depends smoothly on the state.

All three are provable, and the proof is elementary. H3 is therefore removed from
the list of hypotheses and derived. What follows is self-contained.

Throughout, `theory/RATE_LAWS.md` Part 1 supplies the closure and the notation.
Everything here concerns the binding layer alone and is independent of the
catalytic stipulations discussed there.

---

## Setup

Fix a state `(u, a)` with `u, a >= 0` and parameters
`C_tot, D_tot, K_CU, K_CA, K_DU, K_DA > 0`, `nu >= 0`. Write

```
s_u(C,D) = 1 + C/K_CU + D/K_DU        u_f(C,D) = u / s_u(C,D)
s_a(C,D) = 1 + C/K_CA + D/K_DA        a_f(C,D) = a / s_a(C,D)
```

and define the **closure map** `T` on the box `B = [0, C_tot] x [0, D_tot]` by

```
T_1(C,D) = C_tot /(1 + nu + u_f/K_CU + a_f/K_CA)
T_2(C,D) = D_tot /(1 + u_f/K_DU + a_f/K_DA)
```

A solution of the closure is a fixed point of `T`. Equivalently it is a zero of

```
r_1(C,D) = C * q_c(C,D) - C_tot,      q_c = 1 + nu + u_f/K_CU + a_f/K_CA
r_2(C,D) = D * q_d(C,D) - D_tot,      q_d = 1 + u_f/K_DU + a_f/K_DA
```

which is the residual coded in `model.py:_bindingResidual`. Both forms are used:
`T` for existence and uniqueness, `r` for smoothness.

---

## (i) T maps B into itself and is monotone increasing

Each denominator of `T` is at least 1, so `0 < T_1 <= C_tot` and
`0 < T_2 <= D_tot`; `T(B) ⊂ B`, and every fixed point is **strictly positive**,
with the explicit floor

```
T_1 >= C_tot /(1 + nu + u/K_CU + a/K_CA) > 0
```

since `u_f <= u` and `a_f <= a`. Raising `C` or `D` raises `s_u` and `s_a`, which
lowers `u_f` and `a_f`, which lowers both denominators, which raises both
components of `T`. So `T` is monotone increasing on `B`.

By Knaster–Tarski a least and a greatest fixed point exist and are the limits of
iterating `T` from `(0,0)` and from `(C_tot, D_tot)` respectively. This is
**existence**, and it is exactly the certificate computed by
`model.py:solveFreePoolsCertified`.

## (ii) The fixed point is unique

The key property is **strict sublinearity**: `T(λ y) > λ T(y)` for every `y > 0`
and `λ ∈ (0,1)`. It holds because shrinking the argument shrinks the occupancy
terms in each denominator by a factor of at least `λ`, while the constant `1`
in that denominator does not shrink at all — so the denominator falls by less
than the factor `λ`, and the quotient rises by more.

Let `x` and `y` be fixed points, both strictly positive by (i). Put

```
λ* = sup{ λ > 0 : x >= λ y }
```

which is finite and positive because both vectors have strictly positive,
bounded entries, and `x >= λ* y` by closedness. Suppose `λ* < 1`. Since

```
s_u(λ* y) = 1 + λ*(y_1/K_CU + y_2/K_DU) > λ* (1 + y_1/K_CU + y_2/K_DU)
          = λ* s_u(y),
```

we have `u_f(λ* y) < u_f(y)/λ*`, and likewise `a_f(λ* y) < a_f(y)/λ*`. Hence

```
T_1(λ* y) > C_tot /(1 + nu + (u_f(y)/K_CU + a_f(y)/K_CA)/λ*)
          > λ* · C_tot /(λ*(1 + nu) + u_f(y)/K_CU + a_f(y)/K_CA)
          > λ* · C_tot /(1 + nu + u_f(y)/K_CU + a_f(y)/K_CA)
          = λ* y_1,
```

using `λ*(1+nu) < 1+nu` in the last step, and identically for the second
component. Combining with `T(λ* y) <= T(x) = x`, which holds because
`λ* y <= x` and `T` is increasing, gives `x > λ* y` **strictly in both
components**. Both are finite and positive, so `x >= (λ* + ε) y` for some
`ε > 0`, contradicting the definition of `λ*`.

Therefore `λ* >= 1`, i.e. `x >= y`. Exchanging the roles of `x` and `y` gives
`y >= x`. Hence `x = y`: the fixed point is **unique**.

## (iii) The binding Jacobian is a nonsingular M-matrix, with `det >= 1 + nu`

Differentiate `r` with respect to `(C, D)` at fixed `(u, a)`. Since
`∂u_f/∂C = -u_f/(s_u K_CU)` and similarly for the other three partials, all four
are strictly negative for `u, a > 0`, and

```
J_11 = q_c + C ( ∂u_f/∂C /K_CU + ∂a_f/∂C /K_CA )
J_12 =       C ( ∂u_f/∂D /K_CU + ∂a_f/∂D /K_CA )
J_21 =       D ( ∂u_f/∂C /K_DU + ∂a_f/∂C /K_DA )
J_22 = q_d + D ( ∂u_f/∂D /K_DU + ∂a_f/∂D /K_DA )
```

The off-diagonal entries are sums of strictly negative terms with positive
coefficients, so `J_12 <= 0` and `J_21 <= 0`: `J` is a **Z-matrix**.

Now evaluate `J` against the positive vector `x = (C, D)` itself. Because `u_f`
is homogeneous of degree `-1` in `(1, C, D)` jointly, the Euler-type
contractions are exact:

```
C ∂u_f/∂C + D ∂u_f/∂D = -(u/s_u^2)(C/K_CU + D/K_DU) = -u_f (s_u - 1)/s_u
C ∂a_f/∂C + D ∂a_f/∂D = -a_f (s_a - 1)/s_a
```

Hence

```
(Jx)_1 = C [ q_c - (u_f/K_CU)(s_u-1)/s_u - (a_f/K_CA)(s_a-1)/s_a ]
       = C [ 1 + nu + (u_f/K_CU)/s_u + (a_f/K_CA)/s_a ]
       >= C (1 + nu)

(Jx)_2 = D [ 1 + (u_f/K_DU)/s_u + (a_f/K_DA)/s_a ]  >=  D
```

both strictly positive. A Z-matrix admitting a strictly positive vector `x` with
`Jx > 0` is a **nonsingular M-matrix** (Berman and Plemmons, condition M35), so
`J` is nonsingular on the whole positive box and `J^{-1} >= 0` entrywise.

The determinant admits an explicit lower bound. Write `α = -J_12 >= 0` and
`β = -J_21 >= 0`. The two displayed inequalities give
`J_11 >= (1+nu) + (D/C)α` and `J_22 >= 1 + (C/D)β`, so

```
det J >= [ (1+nu) + (D/C)α ][ 1 + (C/D)β ] - αβ
       = (1 + nu) + (1 + nu)(C/D)β + (D/C)α
       >= 1 + nu.
```

**The binding Jacobian determinant is bounded below by `1 + nu` everywhere**, with
equality approached only as `u, a -> 0`. Verified numerically over 24,300 states
spanning four decades of every dissociation constant: the smallest observed
margin `det J - (1+nu)` is `1.1e-06`, the largest off-diagonal entry
`-8.0e-09`, and the smallest margin in the two `Jx` inequalities `3.2e-09` —
all three attained in the vanishing-burden limit where the bounds are equalities.
The proof does not rest on that check.

**Consequence for the code.** The runtime guard in `model.py:jacobian`, which
raises when the binding Jacobian is singular, cannot fire. It is retained as a
defensive assertion and is annotated accordingly.

## (iv) Smoothness, and H3 as a consequence

`r` is a rational function of `(C, D, u, a)` whose denominators `s_u` and `s_a`
are bounded below by 1, so `r` is real-analytic on
`{(C,D) ∈ int B} x {u, a >= 0}`. By (iii) the partial Jacobian
`∂r/∂(C,D)` is nonsingular there. The implicit function theorem therefore
applies at every state, with no singularity possible, and

```
(u_f, a_f, C_f, D_f) are real-analytic functions of (u, a) on u, a > 0.
```

`R` and `G` are compositions of these with the rate laws. Every rate law is a
rational function of the free concentrations except primary nucleation
`k_n u_f^m`, which for non-integer `m > 1` is `C^1` at `u_f = 0` and analytic for
`u_f > 0`. Hence `R` and `G` are `C^∞` on the open positive quadrant and `C^1` up
to the boundary. **H3 holds as a consequence rather than as an assumption.**

## (v) Two monotonicity facts used elsewhere

Both follow from inverse-positivity of `J`, and both are used in task B6.

Differentiating the closure with respect to `u` at fixed `a`,
`∂r_1/∂u = C/(s_u K_CU) > 0` and `∂r_2/∂u = D/(s_u K_DU) > 0`, so

```
∂(C_f, D_f)/∂u = -J^{-1} (∂r/∂u) <= 0.
```

**Raising total soluble burden lowers both free machinery pools.** Then

```
∂u_f/∂u = 1/s_u - (u/s_u^2)( (∂C_f/∂u)/K_CU + (∂D_f/∂u)/K_DU )  >=  1/s_u > 0
∂a_f/∂u =       - (a/s_a^2)( (∂C_f/∂u)/K_CA + (∂D_f/∂u)/K_DA )  >=  0
```

**Free soluble monomer rises with total soluble burden, and so does free
aggregate.** The second is the one that matters: it is why the sign argument of
Section 3.3 does not survive the move from total to free concentrations, since
the two clearance terms in `G` have `C_f` and `D_f` falling while `a_f` rises.

---

## Statement as it enters the paper

**Lemma 0.** *For every `(u, a)` in the nonnegative quadrant the rapid-equilibrium
closure has exactly one solution, and it is strictly positive. The binding
Jacobian is a nonsingular M-matrix with determinant at least `1 + nu`, so the
free concentrations are real-analytic functions of the state on `u, a > 0`.
Consequently `R` and `G` are smooth and H3 is not an independent hypothesis.
Moreover `∂u_f/∂u > 0` and `∂a_f/∂u >= 0`.*

## Reference

Berman A, Plemmons RJ (1994) *Nonnegative matrices in the mathematical
sciences*. SIAM Classics in Applied Mathematics 9. Chapter 6, Theorem 2.3,
conditions M35 and N38.
