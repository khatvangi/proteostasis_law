# theory/ — experiment-grounded rebuild, 2026-10-10

Read in this order:

1. `SCIENTIFIC_THEORY.md` — the argument in plain scientific English
2. `EQUATIONS.md` — conservation equations, theorems, the solved example, proofs
3. `EXPERIMENT_CONSTRAINTS.tsv` — which experiment fixes which assumption; sign-only vs quantitative; held-out vs used
4. `PRIOR_ART.md` — FoldEco / Santra / FoldME scope, and what is left to test
5. `CLAIM_DECISIONS.md` — what was dropped, narrowed, refused
6. `FALSIFICATION_REGISTER.md` — three pre-registered predictions and the blockers
7. `CORRECTIONS.md` — errors found by the independent A6/A7 audit and what replaced them

Verify:  `python checks.py`  (22 checks, ~20 s, writes `checks.json`, exit 0 on pass; reads
`../independent_review/results_a6_a7.json` for one cross-check)

No parameter in this directory is an E. coli estimate. No manuscript is to be drafted until a
prediction in the register has been run. Nothing outside this directory was modified.
