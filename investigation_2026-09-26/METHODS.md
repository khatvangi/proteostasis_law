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
Equilibria at `lambda=2` are roots of
`-chi*x^3 + (1-chi)*x^2 + (rho+1)*x - lambda=0`. Stability is evaluated from
`F'(x)=-g'(x)`: an equilibrium is locally stable when `F'<0`.

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
