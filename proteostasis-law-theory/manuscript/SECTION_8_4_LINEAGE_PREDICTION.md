# §8.4 — The lineage prediction, and what would close it

*Drop-in text. §8.4 does not yet exist in a draft; this is the paragraph for it.
Numbers recomputable by `scripts/phase3/asymmetric_division.py`; see D031-D033.*

---

A dividing cell that partitions aggregate unevenly does not need two attractors
to show a lineage difference. The old-pole cell inherits an inclusion body at
every division, which is a renewed perturbation rather than a basin, and a
monostable system under a renewed perturbation settles to a stationary offset.
Symmetric partitioning in half the volume leaves concentration unchanged, which
is what the dilution term already encodes, so the no-asymmetry control returns an
effect of exactly zero rather than of small magnitude. The difference is therefore
attributable to the asymmetry alone, and no separatrix appears anywhere in the
construction.

Scoring the magnitude requires one more quantity than the model supplies. The
state variable is total aggregate, while the object that segregates indivisibly is
the visible focus. Writing `beta` for the focus share of total aggregate, a
fraction `beta` passes entire to the old-pole daughter and the remainder splits
evenly, so the effective partition is `f = (1 + beta)/2`. This is the scalar rule
reparameterised, not a second mechanism, and it removes `f` from the list of free
quantities: `f` follows from `beta`, and `beta` is physical.

Two further quantities cancel. The growth-burden slope enters as
`(B_old - B_new)/k_mu`, and because `k_mu` is the arrest burden divided by the
quality-control proteome share, this equals the measured slope times the burden
difference expressed as a proteome fraction, identically. The proteome share
cancels and the basal growth rate never appears. What remains is a single
requirement on the aggregate load of the old-pole cell, indexed by `beta`:

| `beta` | old-pole aggregate required, as a fraction of proteome |
|---|---|
| 1.00 | 0.047 % to 0.080 % |
| 0.75 | 0.061 % to 0.104 % |
| 0.50 | 0.091 % to 0.157 % |
| 0.25 | 0.182 % to 0.314 % |
| 0.145 | 0.315 % to 0.541 % |

The requirement scales as `1/beta`, so a smaller focus share demands a larger
aggregate load to produce the same lineage difference. We do not claim this
prediction is parameter-free. No published measurement bounds `beta` under the
unstressed growth the observation concerns. Winkler et al. report, under heat
stress, that 17,500 to 33,000 molecules aggregate per cell and that individual
aggregates hold 2,400 to 16,500 molecules; combining these bounds `beta` between
0.145 and 1, which excludes only the smallest focus shares. That combination is
ours rather than theirs, and the condition is not the one at issue.

The interval is nonetheless informative, because its direction is fixed. Every
departure from a wholly focal aggregate raises the required load and moves it
toward what can be detected. Against the only reported aggregate fraction in
E. coli, 5 % to 10 % of total protein in cells lacking the heat-shock sigma factor
at 30 degrees, the requirement sits 62 to 214 times lower when `beta` is 1 and 9
to 32 times lower at the weakest defensible focus share. Wild-type aggregation at
30 degrees has been reported only as undetected, which is a bound and not a value.

One measurement closes this. Quantitative fluorescence of the polar focus against
total aggregated protein recovered by sedimentation, in the same unstressed cells,
fixes `beta` and collapses the table to a single interval. A second measurement
then decides the prediction: the aggregate content of an old-pole cell, as a
fraction of total protein, to a precision near 0.01 %. We state the requirement
rather than a result, and note that the quantity has not been measured at that
precision in either direction.
