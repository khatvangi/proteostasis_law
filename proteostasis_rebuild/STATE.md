# STATE — running log

Append-only log of loop transitions. Machine-readable mirror: `loop_state.json`.

## 2026-09-27

- scaffold written: WORK_ORDERS.md (gates fixed before analysis), STATE.md,
  STATUS.md, CLAIM_REGISTER.tsv, loop_state.json.
- pre-existing repo state noted: `git status` in proteostasis_law showed the
  `proteostasis-paper` submodule already modified BEFORE this run began. Not
  touched by this run.
- WO-00 started.
- WO-00 PASS (33 claims, anchors verified; 30 legacy files hashed; 7/7 legacy
  numbers reproduced at rel. err 0). One self-caught gate loosening fixed before
  verdict (see WO-00/REPORT.md). WO-01 started.
- WO-01 PASS (legacy dimensionally consistent; non-conservative Phi inflow,
  missing dilution, and a 1/ln2 synthesis-rate overstatement found; C34, C35
  added). WO-02 started.
- WO-02 PASS (conservative model; totals proved symbolically; orthant invariant;
  200 integrations conserve totals; legacy Phi inflow has no donor and is 48%
  of inflow at the legacy threshold). Self-caught hard-coded G2.5 replaced by a
  computed test. WO-03 started.
- WO-03 PASS (exact finite-pool binding; driven cycle set by K_M not K_d;
  four-state cycle shows drive-dependent ultra-affinity with detailed-balance
  control; nascent competition makes theta an output; tQSSA err < 0.5%).
  WO-04 started.
