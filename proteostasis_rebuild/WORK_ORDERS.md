# WORK ORDERS — proteostasis rebuild

Bounded sequential loop, single agent session, started 2026-09-27.
Scope: proteostasis only. Code origin and synonymous-codon evolution are out of
scope (other chats own them).

## loop rules

1. WOs run strictly in order WO-00 → WO-10. A WO may start only when every
   earlier WO is PASS.
2. Deliverables and gates below were written BEFORE any WO analysis and are not
   edited after results are seen. If a gate turns out to be ill-posed, that is
   recorded as a finding in the WO report, not fixed by rewording the gate.
3. PASS requires every listed gate. Any unmet gate → FAIL. A gate that cannot be
   evaluated for lack of an input that this run cannot produce → BLOCKED.
   FAIL or BLOCKED stops the loop.
4. Only new files under `proteostasis_rebuild/`. Legacy files are read-only;
   their sha256 is recorded in WO-00 and re-verified at the end of the run.
   Python runs with `PYTHONDONTWRITEBYTECODE=1` so importing legacy modules
   writes nothing into legacy directories.
5. Simulations are not experimental evidence. No number is called "validated"
   without data. Old headline numbers are not preserved if corrections break
   them.
6. Every WO report ends with an adversarial self-review that tries to falsify
   its own result, and a verdict line `VERDICT: PASS|FAIL|BLOCKED`.

Status vocabulary for claims (CLAIM_REGISTER.tsv):
`REPRODUCED` (legacy code gives the number; says nothing about truth),
`CORRECTED`, `REJECTED`, `CONDITIONAL`, `UNVERIFIED`, `OPEN`.

---

## WO-00 provenance / claim map

Deliverables: `WO-00/REPORT.md`, `CLAIM_REGISTER.tsv` (top level, populated),
`WO-00/legacy_hashes.tsv`, `WO-00/check_provenance.py`, `WO-00/reproduce_legacy.py`.

Gates:
- G0.1 every register row cites an existing legacy file and line, and an
  automated check confirms the quoted anchor text occurs on that line.
- G0.2 every row has a status from the vocabulary and an owning WO.
- G0.3 sha256 manifest of every legacy file read in this run is written.
- G0.4 the legacy headline model numbers (x24.8 headroom at mu = 6.33e-4;
  mechanism at baseline; 1.19e-3 vs 5.66e-3 arithmetic threshold) are re-run
  from legacy code, read-only, and match the stored outputs to rel. tol 1e-6.
  (Reproduction establishes what the legacy computes, not that it is right.)

## WO-01 variables, units, conservation

Deliverables: `WO-01/VARIABLES.md`, `WO-01/units.py` (dimension checker),
`WO-01/legacy_units_audit.py`, tests.

Gates:
- G1.1 every state variable and parameter of the rebuild has a declared unit;
  concentrations in uM, time in s. No "fraction of proteome" state without an
  explicit denominator pool.
- G1.2 an automated dimensional check passes for every term of every equation
  of the legacy two-pool model and reports each term's unit.
- G1.3 conservation laws the rebuild must satisfy are stated as equations
  (total protein, total chaperone) with every source and sink listed.
- G1.4 each legacy term is classified: dimensionally consistent or not;
  conservative flux or not (has a donor and a receiver pool).

## WO-02 source–sink conservative equations

Deliverables: `WO-02/model.py`, `WO-02/EQUATIONS.md`, tests.

Gates:
- G2.1 symbolic (sympy) proof that d(P_T)/dt = s_P − degradation − mu·P_T and
  d(C_T)/dt = s_C − mu·C_T hold identically.
- G2.2 every RHS term passes the WO-01 dimension checker.
- G2.3 nonnegative orthant is forward invariant: symbolic check that each
  dx_i/dt ≥ 0 on the face x_i = 0 for nonnegative parameters.
- G2.4 numerical integration from >= 200 random ICs/parameter sets conserves the
  totals (rel. err < 1e-6 against the analytical total ODE) with no negative
  state below −1e-9.
- G2.5 the legacy Phi inflow term is tested for a donor pool; the phantom mass
  creation rate at the legacy operating point is quantified.

## WO-03 chaperone allocation / kinetics alternatives

Deliverables: `WO-03/chaperone.py`, `WO-03/REPORT.md`, tests.

Gates:
- G3.1 exact finite-pool equilibrium binding recomputed independently; mass
  balance residual < 1e-10 over the sweep; audit values at M_T = 50 uM
  reproduced (C_f = 6.58872).
- G3.2 an ATP-driven cycle model (DnaK-like) is written with explicit
  K_M = (k_off + k_cat + mu)/k_on and shown analytically and numerically to
  differ from the equilibrium K_d benchmark when k_cat is not ≪ k_off.
- G3.3 nascent-chain competition model conserves chaperone and reduces
  exactly to the no-competition model when nascent load is zero.
- G3.4 report separates: equilibrium benchmark / kinetic cycle / competition,
  and names which real machine (DnaK, GroEL, ClpB) each can and cannot stand for.
- G3.5 QSS reductions are checked against the full ODE (rel. err stated).

## WO-04 dynamics, stability, bifurcation

Deliverables: `WO-04/bifurcation.py`, `WO-04/REPORT.md`, tests.

Gates:
- G4.1 the legacy scalar identities (stationary cubic, equilibria, lambda_fold)
  are recomputed symbolically, independently of the audit script.
- G4.2 the legacy operational "collapse" mechanism is classified at the
  evaluation point (true fold vs imposed A_max gate).
- G4.3 for the conservative model, the number of steady states is determined
  by an analytical argument AND by a numerical scan over a declared parameter
  domain; the two agree.
- G4.4 stability of every steady state found is classified by Jacobian
  eigenvalues of the full (non-reduced) system.
- G4.5 any saddle-node found is located by continuation and verified by
  det(J) = 0 with a sign change; if none is found, the report states the
  mechanism that is missing and which variant, if any, restores one.

## WO-05 error-to-burden semantics

Deliverables: `WO-05/error_semantics.py`, `WO-05/REPORT.md`, tests.

Gates:
- G5.1 raw decoding error, amino-acid substitution rate, and MS-detected
  per-substitution rate are defined separately with units and semantics.
- G5.2 every legacy error→flux mapping is classified (correct / double-discount
  / ambiguous) with file:line.
- G5.3 a corrected mapping is implemented; a unit test fails if a
  substitution-level input is passed through (1−S).
- G5.4 the arithmetic threshold identity is recomputed exactly (5.658e-3 at the
  stated parameters) and the legacy 1.19e-3 is explained.
- G5.5 the ×25 headroom is recomputed with the double discount removed.

## WO-06 literature parameter audit

Deliverables: `WO-06/parameter_audit.tsv`, `WO-06/REPORT.md`, checks.

Gates:
- G6.1 every parameter used by the legacy model and the rebuild has: value,
  unit, organism, condition (phase/medium/temperature), source, what was
  actually measured, verification status.
- G6.2 every citation is checked against a bibliographic record (PubMed or
  equivalent); ones that cannot be matched are flagged UNVERIFIED, mismatched
  ones MISCITED. Nothing is marked VERIFIED from memory.
- G6.3 an automated check fails if any row lacks organism/condition/status.

## WO-07 matched-condition bundles

Deliverables: `WO-07/bundles.tsv`, `WO-07/REPORT.md`, checks.

Gates:
- G7.1 separate exponential and stationary bundles; every entry carries a
  condition tag and a match flag (MATCHED / MISMATCHED / UNMEASURED).
- G7.2 no bundle mixes phases silently: an automated check fails if a bundle
  draws a value tagged with the other phase without a MISMATCHED flag.
- G7.3 the effect of each mismatch on a model output is quantified as a range,
  not hidden in a point value.

## WO-08 identifiability and sensitivity

Deliverables: `WO-08/identifiability.py`, `WO-08/REPORT.md`, tests.

Gates:
- G8.1 structural non-identifiability is quantified BEFORE any Monte Carlo:
  rank of the observable sensitivity matrix, by two independent methods
  (symbolic where tractable, numerical SVD at two step sizes) that agree.
- G8.2 each null direction is verified: moving along it changes observables
  by < 1e-6 relative.
- G8.3 Monte Carlo is run only over identifiable combinations, with the
  sampled ranges taken from WO-06/WO-07 and their status carried through.

## WO-09 competing models / discriminating predictions

Deliverables: `WO-09/compare.py`, `WO-09/REPORT.md`, tests.

Gates:
- G9.1 three variants implemented from shared code: legacy (Phi + A_max gate),
  conservative (WO-02/03), adaptive (stress-induced chaperone synthesis).
- G9.2 at least three observables computed for every variant; each observable
  is marked discriminating or not by an explicit, pre-stated criterion.
- G9.3 no statement of validation; every prediction is labelled
  model-conditional.

## WO-10 decisive experiment design

Deliverables: `WO-10/EXPERIMENT_DESIGN.md`.

Gates:
- G10.1 factors, levels, conditions, measurements, and pre-registered
  falsifiers are specified and trace to WO-09 discriminating observables.
- G10.2 sample sizes derive from a stated effect size and noise assumption,
  with the assumption flagged as such.
- G10.3 wet-lab execution is marked HUMAN_REQUIRED; nothing is claimed as run.
