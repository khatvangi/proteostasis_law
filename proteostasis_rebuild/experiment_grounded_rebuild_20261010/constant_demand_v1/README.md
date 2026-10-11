# constant_demand_v1

Steady-state model of keeping a fixed function going under a constant environmental demand.
Intended as the first section of a theoretical proteostasis paper. Not a manuscript.

| file | content |
|---|---|
| DERIVATION.md | the section itself, in prose (draft paper text) |
| EQUATIONS.md | units, equations E1–E18, proofs, ASSUMED vs DERIVED table |
| model.py | scalar model, routing counterexample, multi-client LP with Farkas certificates |
| test_model.py | 12 tests |
| example_output.txt | output of `python model.py` |
| test_output.txt | output of the test run |

Run (numpy, scipy ≥ 1.9 for HiGHS):

    python model.py
    python -m unittest test_model -v

Main results:
- D ≤ a·Pmax/λ is necessary and sufficient only when P is set independently of N on [0, Pmax],
  λ > 0, and there is no route feedback. With feedback it is necessary only.
- For ordinary catalysis the cost of meeting D is turnover, λD/a, not service.
- With the same chaperone counts, poor routing gives P_real = 43.7 against P_req = 60 and
  Pmax = 95, with no machine more than 7% occupied (hypothetical numbers).
- Many clients: A N ≥ D is a linear feasibility problem. Farkas certificates decide it exactly.
  Hall-type subset counting fails once routes are serial.

All parameter values are hypothetical. No E. coli number is claimed.
