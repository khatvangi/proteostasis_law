# A minimal kinetic model of proteostasis maintenance, and the test it fails

Built and run on 10 October 2026, Boron. Every number below comes from `outputs/summary.json`, made by
`python run.py` (seed 20261010, 200 one-client and 200 two-client configurations, about 30 s on 60
cores, rerun byte-identical). Tests: `python -m unittest test_model -v`, 14/14 pass.

## Result first

**The decisive test I could run against existing data goes against the model.** In the model, deleting
HtpG in a ΔdnaKJ background raises steady-state aggregate in 200/200 one-client configurations and
199/200 two-client configurations. Fauvet et al. 2021 measured the opposite at 30 °C: insoluble protein
was 20.5% of total in ΔdnaKJ and 18.9% in ΔdnaKJ ΔhtpG (WT 18.5%)
([DOI](https://doi.org/10.3389/fmolb.2021.653073), full text retrieved from PubMed Central PMC8082187,
saved here as `fauvet2021_PMC8082187_fulltext.txt`). The size is the problem too. The largest increase
any configuration predicts is +1.46 percentage points (median +1×10⁻⁶). Counting HtpG-bound client as
insoluble can flip the model's sign, but measured HtpG abundance caps that effect at 0.635 µM of client,
0.02% of total protein (0.14% even with the seven-fold HtpG induction Fauvet report in ΔdnaKJ), against
an observed drop of 1.6 points.

**Caveat, in the same breath.** The text gives these three values with no spread. The per-replicate data
(five biological replicates per strain) are in a supplement I did not retrieve, so I cannot say whether
18.9 against 20.5 is outside the noise. "Insoluble" also includes membrane proteins, which make up most
of the 18.5% baseline. And my model holds growth rate fixed. If ΔdnaKJ ΔhtpG grows faster at 30 °C,
faster dilution would lower aggregate, and the model cannot express that. Fauvet report faster growth at
37 and 39 °C. I did not find a 30 °C growth comparison. So the result is a **point-estimate
contradiction whose significance is unresolved**, not a refutation.

## What the model is

The model tracks clients, not machines as abstract capacity. Per client (all in µM): free non-native `U`;
held by TF `CT`, by DnaK/J/GrpE `CK`, inside GroEL `CG`, held by HtpG `CH`, engaged by HslUV `CV`;
aggregate `A`; being extracted by ClpB `CB`; native `NAT`; and two cumulative sinks, degraded and diluted.
All flows are mass-action with non-negative constants: synthesis, spontaneous folding `kf·U`,
second-order aggregation `ka·U²`, capture by TF, DnaK, GroEL and HtpG at `kon·F·U` (`F` = free sites),
TF/DnaK→GroEL handoff, release back to `U`, productive exit to `NAT` ("rescue"), HtpG→HslUV transfer and
degradation, ClpB extraction back to `U`, and dilution at the measured growth rate µ on every state.
Every machine has its own finite pool, `F_x = M_x − Σ occupancy`, and no site can serve another
machine's route. The direct GroEL route is present. The model has no viability threshold. Its outputs
are native yield `Y` (native formation flux per synthesized chain) and aggregate `A`. Client mass is
conserved exactly (`d/dt Σ states = s`), and machine pools are conserved by construction.

**Measured inputs (all Schmidt 2016, BW25113 glucose exponential, matched in WO-07):** TF 14.58,
DnaK 10.88, GroEL14 0.950, HtpG2 0.613, ClpB6 0.0117 µM; µ = 1.61×10⁻⁴ s⁻¹; synthesis 0.478 µM/s;
total protein 2965 µM. HslUV is 0.0228 µM, limited by HslV12. It comes from the same Schmidt tables via
`paper_reassessment_20261010/evidence/machine_abundance.tsv`, because WO-07 has no HslUV row. GroES7
(1.24 µM) exceeds GroEL14, so GroEL is taken as limiting. **Everything kinetic is unknown.** It is
sampled log-uniformly in units of µ over ranges listed in `model.RANGES` (e.g. capture hazard on an empty
pool, 1–10⁴ µ). Client share of synthesis is sampled between 0.1% and 10%. Nothing is fitted.

## The one-client analytic core (exact)

With free pools fixed, each held state has a committed probability `a_x` of reaching native without
returning to `U` and `b_x` of returning. Steady state is exact:

- `s = U (Q + k′U)`, with `Q = kf + µ + Σ θ_x(1−b_x)`, `k′ = ka(1−b_A)`, `θ_x = kon_x F_x`
- `U* = 2s / (Q + √(Q² + 4k′s))`
- `Y = P·U*/s`, with `P = kf + Σ θ_x a_x`
- `A* = ka U*² / (kb F_B + µ)`

Free and occupied capacity follow directly. For example `CT = θ_T U*/(r_T+f_T+kh F_G+µ)` and
`F_x = M_x − Σ C_x`, closed self-consistently. The closed form matches ODE integration to ~10⁻¹⁵ and
matches the flux-defined yield to 12 decimal places.

Two exact consequences:

1. **Route criterion.** Raising capture by route x raises yield if and only if `a_x/(1−b_x) > P/(Q + 2k′U*)`.
   In words: the route's chance of success once it commits must beat the marginal yield. Second-order
   aggregation pushes this threshold below Y itself. Checked numerically in 800/800 route evaluations.
2. **A pure holdase cannot raise steady-state yield.** If a chaperone only holds and releases back to the
   same state, it cancels out of the `U` balance. It can only cost yield, through dilution while held and
   through occupancy. Holding helps only when release goes somewhere better (rescue), or transiently.

## Classification of signs

| Perturbation (output) | Result | Class |
|---|---|---|
| faster spontaneous folding, kf ×1.5 (Y, A) | Y up 200/200, A down 200/200; proved at fixed pools | **universal** |
| route with no productive exit, more capture (Y) | never raises Y (proof) | **universal** |
| delete TF, DnaK, HtpG or HslUV (A) | A up 200/200 each, one client | **universal for one client**; HtpG 199/200 with two clients |
| GroEL ×5 in Δtig ΔdnaK (A) | A down 200/200 | **universal given the direct route** |
| GroEL ×5 in Δtig ΔdnaK, serial architecture (Y) | change exactly 0 | **universal**: rescue is inexpressible without the bypass |
| delete HtpG, no HtpG productive exit (Y) | Y up 100/100 | **encoded by assumption**, not a prediction |
| delete HtpG, HtpG productive exit allowed (Y) | up 91, down 9 of 100 | empirically constrained (Fauvet) |
| GroEL ×5 rescues Δtig ΔdnaK yield (Vorderwülbecke) | holds in 150/200 | empirically constrained |
| ΔdnaKJ: ΔhtpG / ΔhslV raise yield (Fauvet growth) | 184/200 and 200/200 | weak constraint; ΔhslV is automatic because HslUV acts only through HtpG |
| all four directional constraints together | 104/200 configs pass | — |
| GroEL ×2 in WT (Y) | up 125, down 75; after constraints 72 vs 32 | **not identified** |
| cooling 30→23 °C, unknown activation energies (Y) | WT 99/101, Δtig ΔdnaK 94/106; after constraints 51/53 and 46/58 | **not identified** |
| faster folding lowers GroEL dependence (Reingewertz direction) | 115 lower, 8 higher of 123 with any dependence | **not identified**: not universal, and no one has measured whether the synonymous change raises kf |
| faster folding lowers DnaK dependence | 139/139 | universal in the sweep (not proved) |

Even requiring the 23 °C direction (Δtig ΔdnaK better cold) leaves the WT cooling sign at 77 up, 17
down. Vorderwülbecke's overproduction rescue does not imply that more GroEL helps wild type.

## Two clients sharing finite pools: sign reversals without optimal routing

The two clients are sampled independently, and routing is pure mass action. Each effect on client 1's
yield was computed with all pools scaled ×0.1, 0.3, 1, 3 and 10.

| Effect on client 1 | always + | always − | reverses with capacity |
|---|---:|---:|---:|
| delete HtpG | 140 | 2 | 58 |
| double client 2's synthesis | 46 | 28 | 126 |
| double client 2's folding rate (a codon change in *another* gene) | 30 | 49 | 121 |
| GroEL ×2 | 83 | 51 | 66 |

The sharpest case: when client 1 has no HtpG productive exit, deleting HtpG never lowers its yield on its
own (0/100). With a competitor present, the deletion lowers it at some capacity in 16/100. The disposal
route had been clearing the competitor off shared TF, DnaK and GroEL. A competing client can also
*raise* client 1's yield by tying up HtpG. These are the reasons single-client signs do not carry over to
a proteome.

## Held-out direction predictions

- **Tested now (D1, above): fails at point estimate, significance unresolved.**
- **Tested now (D2):** ΔdnaKJ Δtig should aggregate more than ΔdnaKJ. The model gives 200/200. Calloni
  2012 (30 °C) report 1,087 against 474 aggregated proteins. Consistent, but those are identity counts,
  not mass. I took them from the reassessment report and did not re-read the source here.
- **Prospective (no data accessed):** (i) because HslUV acts only through HtpG in this model,
  ΔdnaKJ ΔhtpG ΔhslV should equal ΔdnaKJ ΔhtpG. Any extra effect of ΔhslV refutes that structure.
  (ii) If a synonymous mDHFR variant lowers GroEL dependence by speeding folding, it must also lower DnaK
  dependence (139/139). An opposite move would mean the codon acts on route-specific capture, not on
  folding.

## Growth (stated separately, uncalibrated)

`growth / growth_ref = (Y_essential / Y_ref)^β`, with β > 0 unknown (`model.growthIndex`). Only its sign
is used. Growth does not feed back on µ in this model.

## What this does not show

- No kinetic constant is measured. The sweep shows which signs the structure forces, not which values
  hold in cells. Nothing here is "validated", and passing directional constraints is not a fit.
- Structural simplifications: one client per machine unit; ClpB acts without DnaK; HslUV degrades only
  HtpG-delivered clients; no σ32 induction (Fauvet: HtpG ×7 in ΔdnaKJ); fixed growth rate. The
  experiments come from different temperatures (30, 37, 39 °C). The GroEL overproduction fold (×5) is
  assumed, because it is not in our records.
- The prior art FoldEco (2012), Santra et al. PNAS 2017 and FoldME (2017) contain richer versions of this
  network. I did not re-read them, and this model makes no claim of novelty in structure.
- **Without newly accessed held-out datasets the full empirical test remains unresolved.** The
  experiments that would settle it are the Fauvet per-replicate insoluble fractions and 30 °C growth
  rates, a GroEL titration crossed with ΔhtpG and client-resolved proteomics (for the two-client
  reversals), and mDHFR codon variants in ΔdnaK. None of these was retrieved or run here.
