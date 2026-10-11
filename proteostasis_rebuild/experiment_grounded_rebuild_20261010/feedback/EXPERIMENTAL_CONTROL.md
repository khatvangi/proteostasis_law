# How E. coli actually keeps chaperone capacity: the observed σ32 loop

Labels: **LIT** = stated in a primary paper whose record I retrieved from PubMed on 2026-10-10
(abstract, plus full text where marked); **UNMEASURED** = no measured value exists that I could
retrieve for exponential-phase wild-type E. coli; **ASSUMED** = a modelling choice made here.

## The route, step by step

| Step | Status | Evidence |
|---|---|---|
| 1. Stress or a misfolding-prone client raises unfolded substrate | LIT | Tomoyasu 1998: fewer than 1,500 molecules of unstable firefly luciferase raised heat-shock protein levels and prolonged shut-off |
| 2. Unfolded substrate titrates free DnaK/J and GroEL/S | LIT (model of mechanism, strongly supported) | Lim 2013 full text: "extensive support for the model that free chaperones directly inactivate σ32 and that these chaperones are titrated by unfolded proteins" |
| 3. Free chaperone inhibits σ32 activity, directly | LIT | DnaK/J: Tomoyasu 1998. GroEL/S: Guisbert 2004 — GroEL binds σ32 and lowers σ32-dependent transcription with purified proteins; GroEL/S overexpression lowers σ32 activity in vivo, reversed by co-overexpressing a GroEL substrate |
| 4. Free chaperone also destabilises σ32; FtsH degrades it | LIT | Herman 1995: HflB(FtsH) depletion accumulates σ32 and induces the response; σ32 half-life up to 12× longer with reduced HflB, 1.8× shorter with HflB overexpressed. Tomoyasu 1998: lowering FtsH raises σ32 *level* but the extra σ32 is inactive — FtsH sequestration alone cannot induce the response |
| 5. σ32 rises transiently after a temperature upshift | LIT | Straus 1987: σ32 amount rises transiently after 30→42 °C via changes in both synthesis and stability. Lim 2013 (full text): heat also melts an inhibitory rpoH mRNA structure — a **feed-forward** arm, not feedback |
| 6. σ32 regulon products (chaperones, proteases) accumulate after a delay | LIT (qualitative); delay UNMEASURED here | Zhang 2014 (full text): feedback-deficient σ32-I54N raises DnaK, DnaJ, GrpE, GroEL, GroES; Hsp70 set ~6×, Hsp60 set ~4×, ClpX/ClpP/ClpA ~2× after 1 h by SILAC |
| 7. More chaperone → more free chaperone → σ32 shut off | LIT | Tomoyasu 1998: a <1.5-fold rise in DnaK/J lowered σ32 level and activity and sped shut-off; a small fall did the reverse and impaired refolding and growth at high temperature. Lim 2013: DnaK/J/GrpE overexpression dropped σ32 specific activity to 0.3 of wild type |
| 8. Trigger factor is outside the loop | LIT | Zhang 2014: TF "a cotranslational chaperone not under σ32 regulation" was largely unchanged and used as loading control |
| 9. The sensor needs membrane localisation via SRP | LIT | Lim 2013: σ32 must be SRP-targeted to the inner membrane to be properly regulated; a Tn5 insertion in the *ftsY* promoter (FtsY down 3–4-fold) gave σ32 level 6.8× and activity 7.7× of wild type, and activity no longer fell on DnaK/J/GrpE overexpression. The I54N mutant used by Zhang 2014 lies in the "homeostatic control region" and is mislocalised in the same way — which is why it escapes feedback |

Two independent sensors (DnaK/J and GroEL/S) drive one transcription factor; two loops act on it
(activity and FtsH-mediated stability, Lim 2013 full text says "two feedback loops"). The
FtsH-encoding operon is itself σ32-induced (Herman 1995). The 2008 review by the Gross lab says the
central chaperone-mediated loop is still not fully understood (Guisbert 2008).

## What is not measured (and so is not calibrated anywhere in this folder)

- **Loop gain.** No retrieved paper gives the gain from free-chaperone change to steady σ32 output
  in the sense used in MINIMAL_DYNAMICS.md. Tomoyasu 1998 shows the loop is *sensitive* (<1.5-fold
  perturbations matter; <1,500 client molecules are detected). That is qualitative evidence of
  sensitivity, not a number. It is **not** evidence that G > 1 in the sense of MINIMAL_DYNAMICS §3,
  so whether the native loop has a finite delay limit at all is unknown.
- **Delay** from σ32 activation to functional, assembled DnaK/GroEL — UNMEASURED here.
- **Maximum regulon output** in exponential growth without heat. Zhang 2014's ~6×/~4× are with a
  plasmid-overexpressed, feedback-blind mutant in BL21(DE3), an engineered state; they are a lower
  bound on reachable fold change, not the native ceiling.
- **Per-chaperone service rate in vivo** (the `k` of the model) — UNMEASURED; WO-06/WO-07 already
  found the legacy kinetic constants were in-vitro peptide numbers.

## What the 2014 experiment does and does not show

It shows that raising the σ32 regulon without heat (37 °C, minimal growth effect, 89% of detected
proteins changed by <50%) increased soluble and functional yield of two aggregation-prone recombinant
proteins about 2–4-fold and native A25T-TTR tetramer by 39%. So **capacity raised through the σ32
program is usable capacity for a heavy client.** It does not show the native loop reaches that state
on its own: the mutant is chosen precisely because wild-type σ32 output is transient and limited by
feedback and degradation (stated in their introduction). The ~6× vs ~4× fold changes also mean the
"evolved stoichiometry" is approximately, not exactly, preserved.

## Sources (all retrieved via PubMed, 2026-10-10)

- Straus, Walter & Gross 1987 Nature 329:348 — [10.1038/329348a0](https://doi.org/10.1038/329348a0)
- Herman, Thévenet, D'Ari & Bouloc 1995 PNAS 92:3516 — [10.1073/pnas.92.8.3516](https://doi.org/10.1073/pnas.92.8.3516)
- Tomoyasu, Ogura, Tatsuta & Bukau 1998 Mol Microbiol 30:567 — [10.1046/j.1365-2958.1998.01090.x](https://doi.org/10.1046/j.1365-2958.1998.01090.x)
- Guisbert, Herman, Lu & Gross 2004 Genes Dev 18:2812 — [10.1101/gad.1219204](https://doi.org/10.1101/gad.1219204) (abstract only; PMC returned no body text)
- Guisbert, Yura, Rhodius & Gross 2008 MMBR 72:545 — [10.1128/MMBR.00007-08](https://doi.org/10.1128/MMBR.00007-08) (abstract only)
- Lim et al. 2013 PLoS Biol 11:e1001735 — [10.1371/journal.pbio.1001735](https://doi.org/10.1371/journal.pbio.1001735) (full text)
- Zhang et al. 2014 ACS Chem Biol 9:1945 — [10.1021/cb5004477](https://doi.org/10.1021/cb5004477) (full text)
