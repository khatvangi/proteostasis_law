# B1 resolved — ODE vs arithmetic bound direction

**Question (blocker B1).** `proteostasis-P1/README.md` states the ODE bound sits
"~30–100× **below**" the arithmetic bound (amplification closes the window early);
the latest `two_pool_summary.md` shows it ~10–100× **above** the [1e-4,1e-3]
window. Which is right, and can the two-pool result enter the combined paper?

## Answer: the README is stale. The ODE bound sits ABOVE (looser than) arithmetic.

Traced to the raw JSON (`two_pool_results.json`, `arithmetic_results.json`), not
the prose. Three reference points, low → high:

| Quantity | Value (/codon) | Source |
|---|---|---|
| Observed E. coli mistranslation rate | ~1×10⁻⁴ – 1×10⁻³ | literature consensus |
| Arithmetic tolerable bound | 1.19×10⁻³ point; 8.5×10⁻³ MC median | `arithmetic_results.json` A_reproduction, D_mc |
| ODE operational bound | 1.00×10⁻² deterministic; 1.70×10⁻² MC median | `two_pool_results.json` B_compare, D_monte_carlo |

**The two bounds are statistically consistent, with the ODE looser:**
- ratio arith/ODE median = 0.50 (ODE ~2× above arithmetic median)
- paired MC P(arith < ODE) = 0.654
- 95% CI overlap = 100% of the arithmetic CI width
- ODE / arithmetic-window: ×10 (upper 1e-3) to ×100 (lower 1e-4)

## Why the README says the opposite

The README describes a **monomer-runaway** regime in which chaperone-titration
amplification Φ(P) collapses folding *before* the per-protein arithmetic limit,
which would push the ODE bound *below* arithmetic. That is not the operative
mechanism in the current model. At baseline and in 99.94% of MC draws
(4,997 / 5,000) the bound is set by **aggregation_death** — the aggregated
fraction A reaching A_max = 0.25 — not by monomer runaway (3 / 5,000 draws).
Aggregation-death gives a *looser* bound than runaway would, so the ODE bound
ends up above arithmetic. The README narrative predates aggregation-death
becoming the dominant mechanism.

## The empirical-N update did NOT cause the flip

The deterministic two-pool bound is **identical** (1.000×10⁻²/codon,
aggregation_death) in both `two_pool_summary.backup_uniformN.md` and the current
`two_pool_summary.md`. Empirical-N bootstrap only widened the Monte Carlo
envelope (MC median moved, CIs widened); it did not move the deterministic bound
or change the direction. The README was already stale before empirical-N.

## A second confusion in the summary text

`two_pool_summary.md`'s "gap diagnosis" labels [1e-4, 1e-3] the "arithmetic bound
window." That interval is actually the **observed E. coli rate** window, not the
arithmetic bound (whose MC 95% CI is [1.24×10⁻³, 7.89×10⁻²], median 8.5×10⁻³).
Conflating the observed-rate window with the arithmetic bound inflates the
apparent gap. State the three quantities separately.

## The physically meaningful result (Part E)

At the observed E. coli rate f = 1×10⁻⁴/codon, the system sits on the low-burden
stable branch with **×158 headroom** to the collapse threshold P_dagger
(P* = 1.7×10⁻⁴ vs P_dagger = 2.7×10⁻²) and ×11,091 headroom on the aggregated
pool. This is the defensible, citable statement: the observed error rate sits
roughly two orders of magnitude inside the proteostasis envelope, and all three
independent estimates (observed, arithmetic, ODE) agree to within ~2 orders of
magnitude.

## Actions

1. **README correction (needs rw; `proteostasis-P1` is read-only here).** Replace
   "order-of-magnitude consistent with the arithmetic bound but sits ~30–100×
   below it because the amplification factor Φ(P) closes the feasible window"
   with: "order-of-magnitude consistent with the arithmetic bound and sits ~2×
   above it at the median (fully overlapping 95% CIs); the operative closure
   mechanism is aggregation-death (A→A_max), not monomer runaway."
2. **Combined draft:** B1 is resolved → Result 5b (two-pool consistency) may now
   be added with the corrected direction and the ×158 headroom framing. Added
   below.
