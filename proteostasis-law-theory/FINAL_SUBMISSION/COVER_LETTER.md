# Cover letter

To the Editors, *Bulletin of Mathematical Biology*

Dear Editors,

I submit **"An Exact Fold Condition for Mass-Balanced Models of Protein Quality
Control"** for consideration as a research article.

A cell clearing misfolded protein with finite chaperone and protease pools
tolerates damage influx only up to a threshold. Where that threshold sits is
normally located by numerical continuation, one parameter set at a time. The
paper shows it is algebraic. For any clearance network in which the damage
influx is state-independent and mass balance accounts for all outflow, the
Jacobian determinant factors as `det J = −(∇R × ∇G)`, with `R` the total removal
flux and `G` the aggregate nullcline field. A fold is therefore exactly a
constrained critical point of removal on the aggregate nullcline, with critical
influx `j_crit = R(u*, a*)`. The proof needs mass balance and one
determinant-preserving row operation, and nothing else; the same row operation
gives the `n`-state form, so a model in this class can be extended without
re-deriving its fold condition.

Two things bound the result and are stated where the claim is made. The converse
holds under three genericity conditions — a nonvanishing constraint gradient, a
simple zero eigenvalue, and a nondegenerate second derivative of removal along
the nullcline — together with a fourth that is proved from the first rather than
assumed; the nondegeneracy condition turns out to be Sotomayor's, up to sign,
which makes the classification a corollary rather than a separate appeal.
Second, solving the condition returns a candidate *set*, not a candidate:
identifying which candidate terminates the accessible low-burden branch requires
branch information in 9.1% of the kinetic ensemble, and the paper reports that
rather than treating solver convergence as identification.

Beyond the identity, the paper establishes that it survives growth dilution and
clearance machinery damaged by the influx it clears; quantifies how far below
the capacity ceiling the fold sits, with the closure-dependence of that factor
carried at the claim; and classifies the bifurcation structure preceding the
fold, including a set of networks that lose stability and regain it under
monotone loading, with the interior of that band integrated rather than
interpolated. Three codimension-two boundaries are located, each approached by a
different network.

The work is theoretical and computational; no new experimental data were
generated. All code, configurations, and the test suite asserting every reported
quantity are public under the MIT licence and archived on Zenodo
(doi:10.5281/zenodo.21794565). One limitation of the deposit — the generating
run's working tree was uncommitted at launch — is stated in the availability
section together with its scope.

The manuscript is not under consideration elsewhere, and a preprint is deposited
on bioRxiv. I declare no competing interests.

Yours sincerely,

Kiran Boggavarapu  
ORCID 0000-0003-0751-6459  
Department of Chemistry and Physics, McNeese State University  
Lake Charles, LA 70609, USA  
kiran@mcneese.edu
