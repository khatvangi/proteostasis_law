# Methods

## Inputs and non-destructive scope

The redraw is standalone and does not import or modify the legacy ODE. Baseline
values were transcribed read-only from `../proteostasis-P1/two_pool_ode.py`:
`C_tot=50 µM`, `K_d=1 µM`, `Prot_tot=300 µM`, `k_obs_max=0.01 s^-1` in the
source. The redraw sets `k_obs_max=1` for the normalized rate panel because a
common multiplicative rate constant does not change the closure comparison.
The ligand sweep is `M_T` from 0 to 300 µM, inclusive, which is nonnegative and
does not exceed the source total-proteome concentration. It is a mathematical
feasibility sweep, not a claim that every point is a measured physiological
state.

## Scalar model

Use

`F(x)=lambda-x-rho*x/(1+x)+chi*x^2`, `g(x)=x+rho*x/(1+x)-chi*x^2`,

with `rho=4`, `chi=0.15`, and `F=lambda-g`. Stationary points are roots of
`g'(x)=1+rho/(1+x)^2-2 chi x`; the script obtains them from the expanded cubic
after multiplying by `(1+x)^2`, then filters real roots with `x>0`.
Equilibria at `g(x)=lambda` are obtained by multiplying by `(1+x)`:
`-chi*x^3 + (1-chi)*x^2 + (1+rho-lambda)*x - lambda=0`.
For `rho=4`, `chi=0.15`, and `lambda=2`, the positive roots are
`x=0.5808674541` and `x=7.9669676044`. Each retained root is independently
checked against the original rational function with `abs(g(x)-lambda)<1e-9`.
The third algebraic real root is `x=-2.8811683919` and is outside the
nonnegative physical domain; it also satisfies the cleared equation and the
original rational function.
Stability is evaluated from `F'(x)=-g'(x)`: an equilibrium is locally stable
when `g'(x)>0` (equivalently `F'<0`) and unstable when `g'(x)<0`.
The positive nonzero root for `lambda=0` is `x=9.2645937939`; it is obtained
from the same polynomial, not from an `x=10` marker.

## Exact finite-pool binding

Let total chaperone and ligand be `C_T` and `M_T`, with free species `C_f` and
`M_f`, and complex `C_b`. Conservation and mass action are

`C_T=C_f+C_b`,

`M_T=M_f+C_b`,

`K_d=C_f M_f/C_b`.

Eliminating free species gives

`C_b^2-(C_T+M_T+K_d)C_b+C_T M_T=0`.

The physical root is the smaller root. `redraw.py` evaluates it in a stable
form,

`C_b = 2 C_T M_T / (C_T+M_T+K_d + sqrt((C_T+M_T+K_d)^2-4 C_T M_T))`,

then returns `C_f=C_T-C_b`. The legacy approximation is evaluated as
`C_T/(1+M_T/K_d)`. Both are inserted into the same source-inspired functional
form `v_fold=k_obs_max C_f/(C_f+K_d)`.
The exact/approximate rate ratio is checked at representative loads; for
`M_T=0, 10, 50, 300 µM` it is approximately `1.00000, 1.19042, 1.75382,
1.16534`, respectively.

## Plots and files

`scalar_g_equilibria.png/svg` plots `g`, its positive stationary maximum, the
positive zero at `x=9.2646`, and the stable/unstable positive equilibria for
`lambda=2`; the horizontal `lambda=2` line makes the equilibrium construction
explicit. `finite_pool_vs_approx.png/svg` plots free chaperone and normalized
folding rate for exact and approximate closures over the concentration sweep.

## Reproduction

From this directory:

```bash
python redraw.py
pytest -q
```

Dependencies are NumPy, Matplotlib, and pytest. No fitting, literature-derived
parameter optimization, or experimental-data calibration is performed.
