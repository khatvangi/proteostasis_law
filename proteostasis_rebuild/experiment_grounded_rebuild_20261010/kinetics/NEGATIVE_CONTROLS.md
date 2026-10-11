# Negative controls

Each control is a test in `test_model.py` (14/14 pass) and, where noted, a tally in
`outputs/summary.json` from `run.py`.

| Control | What would show an artifact | Result |
|---|---|---|
| **Single pooled-capacity comparator** | if pooling capacity gave the same signs, resolving machines would add nothing | Pooled and resolved models are identical at vanishing load (test). At measured pools they disagree in sign for HtpG deletion in 19/200, GroEL ×2 in 69/200, and GroEL rescue of Δtig ΔdnaK in 52/200 two-client configs |
| **Same-total machine swap** (0.5×GroEL14 sites moved from DnaK to GroEL) | a pooled model cannot see it | Resolved model: client-1 yield changes in 200/200. Pooled model: 0/200, identical to 13 decimals (test) |
| **No HtpG harm when route disabled** (HtpG capture = 0) | an HtpG effect without an HtpG route would be leakage | Deleting HtpG changes yield by exactly 0 (test, 5 configs × 2 clients) |
| **Direct GroEL bypass required** | if serial TF/DnaK→GroEL could rescue, the bypass would be unnecessary | Serial: GroEL ×5 changes Δtig ΔdnaK yield by exactly 0 (200/200, max abs 0.0). With bypass it moves (test) |
| **No ClpB aggregate quota** | the old model rescaled demand to the tiny ClpB pool | Client demand is set only by its synthesis share. Deleting ClpB leaves every parameter unchanged and gives Y > 0, A > 0 (test). Key signs (GroEL ×2, ΔhtpG, faster folding) are unchanged at ClpB ×≈0, ×1 and ×10 in 197/200 |
| **Mass conservation** | lost or created client | `Σ states = s·t` to 10⁻⁸ relative over 50/µ (test, 3 two-client configs) |
| **Positivity and capacity** | negative states or overfilled pools | ODE minimum ≥ −10⁻¹²; F ≥ 0; occupancy ≤ pool; 0 ≤ Y ≤ 1 (test) |
| **Exactness** | algebraic steady state wrong | Closed form = ODE steady state to ≤10⁻⁶ relative, = flux yield to 12 decimals (test) |
| **Start independence** | multiple steady states | Two very different initial conditions reach the same state to 10⁻⁸ (one config; not a proof of uniqueness) |
| **Encoded-sign check** | a constraint counted as a prediction | HtpG-deletion benefit with no HtpG productive exit is labelled *encoded* (100/100 by construction). ΔhslV benefit is automatic for the same reason |
| **Reproducibility** | hidden randomness | Fixed seeds; `summary.json` byte-identical on rerun |

Solver note: the first full run failed to converge on 4/200 two-client configs (residual up to 10⁻²).
Gauss–Seidel stalls on the slow HtpG→HslUV mode. I added an ODE-relaxation fallback before the Newton
polish, and all 400 configs now converge (residual < 10⁻⁹). Non-converged configs were not dropped.
