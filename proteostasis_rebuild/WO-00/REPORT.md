# WO-00 — provenance / claim map

## deliverables and gates (copied from WORK_ORDERS.md, fixed before analysis)

- G0.1 every register row cites an existing legacy file+line; anchor text verified on that line.
- G0.2 every row has a vocabulary status and an owning WO.
- G0.3 sha256 manifest of every legacy file read.
- G0.4 legacy headline numbers re-run read-only, match stored outputs to rel. tol 1e-6.

## legacy sources inspected (read-only)

- `proteostasis-P1/two_pool_ode.py` (the model; vendored byte-identical at
  `envelope-paper/scripts/vendor/two_pool_ode.py`, confirmed with `diff`)
- `proteostasis-P1/LITERATURE_ANCHORS.md`, `arithmetic_*`, `paired_mc*`, `two_pool_*`
- `envelope-paper/manuscript/MANUSCRIPT.md`, scripts 06/09/11/12, computed JSONs
- `investigation_2026-09-26/{AUDIT,METHODS,EVIDENCE_AND_FALSIFICATION}.md`

## result

`CLAIM_REGISTER.tsv` holds 33 claims (C01–C33): 15 from the manuscript, 18 from
legacy code, anchors, and computed outputs. All start as `OPEN`; each is owned by
the WO that must settle it. Statuses are updated as later WOs close.

Reproduction (`reproduce_legacy.py` → `legacy_reproduction.json`):

| check | recomputed | stored | rel. err |
|---|---:|---:|---:|
| headroom_P at usage-weighted mu (6.334e-4) | 24.81728 | 24.81728 | 0 |
| headroom_A at usage-weighted mu | 274.6522 | 274.6522 | 0 |
| headroom_P at 1e-4 | 158.0557 | 158.0557 | 0 |
| two-pool f_codon_crit (baseline) | 1.000321e-2 | 1.000321e-2 | 0 |
| two-pool P_dagger | 2.740237e-2 | 2.740237e-2 | 0 |
| arithmetic, exact, stated params | 5.658143e-3 | 5.658143e-3 | 0 |
| arithmetic with (1−S)p_m forced to 1 | 1.188210e-3 | 1.188210e-3 | 0 |

Findings that WO-00 surfaces (settled later, not here):

1. At baseline the legacy "collapse" mechanism is `aggregation_death`: the
   threshold is where the uncapped quasi-steady A reaches the imposed
   `A_max = 0.25`, not a fold of the dynamics (C07; owner WO-04).
2. The legacy's own stored `A_reproduction` block has `"factor": 1.0`, i.e.
   the 1.19e-3 "reproduction" set `(1−S)·p_misfold = 1`, while the manuscript
   (line 104) states p_misfold = 0.3 and S = 0.3. Under the stated parameters the
   legacy code itself gives 5.658e-3 (its own `B_length_sweep`, N=300). (C05, C06)
3. Scripts 09/11/12 pass the MS-derived mu through `(1 − S)`
   (C22–C24). Whether that is a double discount is WO-05's question.
4. `supraadditivity_summary.json` labels the capacity perturbation
   `C_tot_uM / factor` while the code default and the same file's
   `capacity_knob` field say `k_obs_max` (C32). The docstring of
   `09_supraadditivity.py` (lines 109–113) also quotes operating-point values
   (M = 0.052 uM, 97.9%) that belong to f = 1e-4, while the summary at the
   evaluation point reports 0.331 uM and 97.4%. Stale text, not a numerical error.
5. The pre-existing repo state: `proteostasis-paper` submodule was already
   modified before this run (recorded in STATE.md).

## deterministic checks run

```
PYTHONDONTWRITEBYTECODE=1 python check_provenance.py --hash   # 30 files hashed
PYTHONDONTWRITEBYTECODE=1 python check_provenance.py          # 33 claims, 0 failures
PYTHONDONTWRITEBYTECODE=1 python reproduce_legacy.py          # 7/7 PASS at 1e-6
```

The vendored legacy `__pycache__/two_pool_ode.cpython-312.pyc` still carries its
Jul 30 timestamp after import, so the import wrote nothing into the legacy tree.

## adversarial self-review

- *Could the anchors pass while pointing at the wrong thing?* The check only
  proves the anchor string is on the cited line. Several anchors are short
  (`7.4%`, `12 of 36`, `(1.0 - p.S_avg)`); a short anchor on the right line
  cannot be the wrong claim, but a row's *interpretation* column is my reading
  and is not machine-checked. Those readings are what later WOs test.
- *Did I weaken a gate?* Yes, in the first draft: I compared against 4-figure
  printouts with a 5e-4 tolerance. That silently loosened G0.4. I replaced
  those references with the full-precision stored JSON values and restored the
  uniform 1e-6 tolerance; all seven pass at 0 error.
- *Is rel. err = 0 suspicious?* No: the headroom rows run the same code on the
  same inputs, so bit-identical output is expected. The arithmetic rows use a
  formula I wrote independently and still agree to machine precision. The
  reproduction proves determinism and that I have located the generators; it
  does NOT prove any of these numbers is scientifically right.
- *Is the claim list complete?* No, and it is not required to be. It covers the
  abstract's numeric claims and every legacy assumption named by the prior audit.
  Codon-axis claims (mu/nu clustering, z-scores) are out of scope (other chats).

VERDICT: PASS — G0.1, G0.2, G0.3, G0.4 all met.
