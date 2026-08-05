# Theorem 2, proved

Deliverable of task B2. Through v5 Theorem 2 was **asserted**: the manuscript
stated the converse, listed four genericity conditions, and verified them
numerically over two ensembles. Numerical verification of (G1)–(G4) is evidence
that the conditions hold; it is not a proof that they imply the normal form.

The proof is half a page and needs no citation to Sotomayor, because the
nondegeneracy condition the standard theorem asks for turns out to be the same
object as (G3), sign included. That coincidence is the content of Step 4 and is
the reason the arclength formulation was worth introducing in Theorem 1.2.

Notation follows the manuscript: `F = (F_1, F_2)` with `F_1 = j − R − G` and
`F_2 = G`, so that

```
J = [ -R_u - G_u    -R_a - G_a ]
    [   G_u            G_a     ]
```

and at a constrained critical point `grad R = lambda grad G` for some scalar
`lambda`. Write `v` for a right null vector of `J`, `w` for a left one, `gamma`
for the arclength parametrisation of `{G = 0}` through the point, and
`r(s) = R(gamma(s))`.

---

## Statement

**Theorem 2 (converse).** *Let `F` be `C^2`, and let `x* = (u*, a*)` satisfy
`G(x*) = 0` and `det J(x*) = 0`. Assume*

* **(G1)** `grad G(x*) != 0`;
* **(G2)** `tr J(x*) != 0`;
* **(G3)** `r''(0) != 0`.

*Then `x*` is a generic saddle-node in the influx `j`, with
`j_crit = R(u*, a*)`.*

Parameter transversality is not assumed. It is automatic: see
`theory/G4_FROM_G1.md` (task B3), which shows `w · F_j = 1` identically under the
normalisation `w_1 = 1`.

---

## Step 1. The equilibrium count

By Theorem 1.2, `grad G != 0` makes `{G = 0}` a `C^2` curve near `x*` and the
equilibria on it are exactly the solutions of `r(s) = j`. Also by Theorem 1.2,
`r'(s) = det J/||grad G||`, so `det J(x*) = 0` is `r'(0) = 0`. Taylor:

```
r(s) = j* + (1/2) r''(0) s^2 + O(s^3),      j* = r(0) = R(u*, a*)
```

With `r''(0) != 0` this is a nondegenerate extremum, so for `j` on one side of
`j*` there are exactly two equilibria near `x*` and on the other side none. Two
equilibria collide and annihilate as `j` crosses `j*`. That is the saddle-node in
the equilibrium-count sense, and it already fixes the critical value.

The sign of `r''(0)` says which side. `r''(0) < 0` is a constrained maximum of
removal: solutions exist below `j*` and vanish above, so the pair is destroyed as
the influx rises — collapse-oriented. `r''(0) > 0` creates the pair as the influx
rises — birth-oriented, the lower turning point of a hysteresis loop. This is the
orientation classification of Section 3.4, and it is a corollary of Step 1 rather
than a separate observation.

## Step 2. The centre manifold is one-dimensional

In two dimensions `det J = 0` makes `0` an eigenvalue, and the other eigenvalue
is `tr J`. (G2) therefore does two things at once: it makes the zero eigenvalue
algebraically simple, and it puts the remaining eigenvalue off the imaginary
axis. The equilibrium is hyperbolic in the complementary direction, and the
centre manifold is one-dimensional.

(For `n >= 3` neither consequence follows from `tr J != 0`, which is why
Corollary 2's converse needs a different condition — task B4.)

## Step 3. The centre direction *is* the nullcline tangent

Let `v` span `ker J`. The second row of `Jv = 0` reads

```
G_u v_1 + G_a v_2 = 0,   i.e.   grad G . v = 0
```

so `v` is orthogonal to `grad G`, hence tangent to `{G = 0}`. The unit tangent
`gamma'(0) = (-G_a, G_u)/||grad G||` is orthogonal to `grad G` too, and in two
dimensions that orthogonal complement is one-dimensional, so

```
v  ∥  gamma'(0).
```

**The centre manifold and the arclength curve are the same object.** This is why
Step 4 can compare a centre-manifold quantity with a derivative of `r` at all:
they are computed along the same direction. Normalise `v = gamma'(0)`.

## Step 4. (G3) is Sotomayor's nondegeneracy condition, sign included

The left null vector. With `grad R = lambda grad G`, the vector `w = (1, 1+lambda)`
satisfies

```
wJ = ( -R_u - G_u + (1+lambda)G_u ,  -R_a - G_a + (1+lambda)G_a )
   = ( -R_u + lambda G_u ,  -R_a + lambda G_a )
   = -( grad R - lambda grad G ) = 0.
```

The second derivative of `F` against `w`. Since `F_1 = j - R - G` and `F_2 = G`,

```
w . D^2F(v,v) = D^2F_1(v,v) + (1+lambda) D^2F_2(v,v)
              = -D^2R(v,v) - D^2G(v,v) + (1+lambda) D^2G(v,v)
              = -D^2R(v,v) + lambda D^2G(v,v).                        (A)
```

The second derivative of `r`. Differentiating `r(s) = R(gamma(s))` twice,

```
r''(0) = D^2R(v,v) + grad R . gamma''(0).
```

The curve lies in `{G = 0}`, so `G(gamma(s)) = 0` identically; differentiating
that twice gives `D^2G(v,v) + grad G . gamma''(0) = 0`, i.e.
`grad G . gamma''(0) = -D^2G(v,v)`. Using `grad R = lambda grad G`,

```
grad R . gamma''(0) = lambda grad G . gamma''(0) = -lambda D^2G(v,v)

r''(0) = D^2R(v,v) - lambda D^2G(v,v).                                (B)
```

Comparing (A) and (B):

```
r''(0) = - w . D^2F(v,v).
```

Sotomayor's saddle-node theorem asks for `w . F_j != 0` and
`w . D^2F(v,v) != 0`. The first is automatic (task B3). The second is exactly
`r''(0) != 0`, which is (G3) — not merely equivalent to it but equal to it up to
sign. No black-box citation is required, and the sign convention is fixed rather
than inherited: a **negative** `w . D^2F(v,v)` is a constrained **maximum** of
removal and a collapse-oriented fold.

## Step 5. Conclusion

Steps 1–4 give a one-dimensional centre manifold on which the reduced equation
is, in the arclength coordinate,

```
ds/dt  ~  c ( j - j* ) - (1/2) |r''(0)| s^2 + O(s^3),   c != 0,
```

the standard saddle-node normal form, with the quadratic coefficient carrying the
sign of `-r''(0)`. Together with the equilibrium count of Step 1 and the
hyperbolicity of Step 2, `x*` is a generic saddle-node in `j` at
`j_crit = R(u*, a*)`. ∎

---

## What the numbers are now for

Table 1 previously carried the weight of the converse. It no longer does. What it
reports is whether the hypotheses **hold in the ensembles studied**, which is a
question about the model rather than about the theorem — and the two networks
that fail do so by failing (G1), where the constraint stops being a manifold and
`{G = 0}` has no tangent to differentiate along. Step 1 is the first thing that
breaks there, not the last.
