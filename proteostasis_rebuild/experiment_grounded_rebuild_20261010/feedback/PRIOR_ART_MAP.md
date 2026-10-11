# Prior art, what is new, and how this feeds the codon question

## Already established (do not present as ours)

| Result | Where it is already |
|---|---|
| Free DnaK/J and GroEL/S inhibit σ32; unfolded clients titrate them; σ32 induces more chaperone — a chaperone-mediated negative feedback loop | Tomoyasu 1998; Guisbert 2004; reviewed Guisbert 2008 |
| FtsH degrades σ32, and is itself σ32-induced; second (stability) loop | Herman 1995; Tomoyasu 1998 |
| Heat feed-forward via rpoH mRNA thermosensor; σ32 rise is transient | Straus 1987; summarised in Lim 2013 |
| SRP-dependent membrane localisation is needed for σ32 to be regulated | Lim 2013 |
| Mathematical models of the E. coli heat shock loop decomposed into feedforward / feedback modules with robustness–performance trade-offs | El-Samad, Kurata, Doyle, Gross & Khammash 2005 PNAS 102:2736 [10.1073/pnas.0403510102](https://doi.org/10.1073/pnas.0403510102); Kurata et al. 2006 PLoS Comput Biol 2:e59 [10.1371/journal.pcbi.0020059](https://doi.org/10.1371/journal.pcbi.0020059) |
| Raising the σ32 regulon (feedback-blind I54N) raises usable folding capacity for heavy clients; TF not in regulon | Zhang 2014 |
| Growth-rate dilution can turn a circuit bistable when expression slows growth | Klumpp, Zhang & Hwa 2009 Cell 139:1366 [10.1016/j.cell.2009.12.001](https://doi.org/10.1016/j.cell.2009.12.001); Tan, Marguet & You 2009 Nat Chem Biol 5:842 [10.1038/nchembio.218](https://doi.org/10.1038/nchembio.218) |
| Protein allocation as a bounded proteome fraction under growth | Scott et al. 2010 Science 330:1099 [10.1126/science.1192588](https://doi.org/10.1126/science.1192588) |
| Finite-capacity chaperone networks with folding/degradation/aggregation competition, growth coupling and collapse | FoldEco, Santra/Farrell/Dill 2017, FoldME — listed by the sibling `theory/PRIOR_ART.md`; not re-retrieved here |
| Proportional feedback reduces error by 1+loop gain; high gain plus delay gives Hopf oscillation; actuator saturation caps what any controller can deliver | textbook control theory. The exact two-variable statement (stable at every delay iff G ≤ 1) follows from the standard crossing-frequency method for retarded delay equations; I have not checked whether the El-Samad/Kurata models state it |

I did not read the El-Samad/Kurata model equations, only their abstracts, so I cannot say which of the
algebra in MINIMAL_DYNAMICS.md they already contain. Assume most of §1–§3 is there in some form.

## What is (modestly) new here

Little, and none of it is new biology. It is framing that the project needs:

1. **σ32 senses spare capacity directly** (spare capacity = k × free chaperone). Simple identity, but it
   is the bridge between the static inequality and the dynamic controller.
2. **An explicit three-way split** — admissible vs implemented vs productive — with the bookkeeping
   identity L = served + degraded + diluted, and a proof (constant growth, fixed production) that raising
   degradation always lowers chaperone-served flux while the balance still closes. "Served" is chaperone
   turnover, not folded output; the model has no native-product state. This is what the earlier static `B_total ≤ C_buffer` framing
   could not distinguish.
3. **The one-input / two-machine rank limit**: σ32 can move pools only along one curve, so an
   admissible allocation needing a different GroEL:DnaK ratio is unreachable by this controller, and
   TF-dependent routes are uncontrolled. I have not seen this stated for σ32 specifically; I have not
   searched exhaustively either.

None of these is a quantitative prediction about E. coli. No gain, delay, ceiling or service rate is
measured, so no statement of the form "E. coli sits X-fold inside its limit" follows from this folder.

## How this becomes the foundation for the codon-conditioned demand question

The downstream question is whether codon choice changes *which machine* and *how much* each client
needs (the reassessment's surviving evidence, as cited in `paper_reassessment_20261010/SCIENTIFIC_REPORT.md` and not re-retrieved here:
Reingewertz 2024 synonymous DHFR variants change GroEL dependence; serine AGC decoding difference). This folder fixes three things that question must respect:

1. **Demand must be machine-resolved, not summed.** Because σ32 regulates DnaK and GroEL along a single
   curve, and TF not at all, a codon change that shifts a client from DnaK/TF to GroEL changes viability
   differently from one that raises total load by the same amount. So the codon effect must be
   expressed as a change in a demand *vector* (per machine), and tested against the σ32 reachable curve.
2. **The relevant readout is σ32 output, which is measurable.** If codon changes alter chaperone demand,
   free chaperone falls and σ32 activity rises — a σ32-reporter (or DnaK/GroEL level) shift in isogenic
   synonymous variants is the direct, parameter-free test. Same-total-load variants that shift
   machine demand should move σ32 differently depending on whether GroEL or DnaK is the titrated sensor.
3. **Fate, not balance, is the endpoint.** Since the loop also raises proteases, a codon effect can be
   absorbed by degradation. A codon claim therefore needs folded and degraded fractions of the client
   measured, not just solubility or growth.

What this folder does **not** do: it does not show codon choice affects chaperone demand in vivo for
any native E. coli protein, and it supplies no number for the codon-demand coupling.
