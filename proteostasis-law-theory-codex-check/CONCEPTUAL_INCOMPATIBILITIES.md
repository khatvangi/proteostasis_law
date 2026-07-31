# Conceptual incompatibilities and limits

- The canonical law permits equilibria or other bounded attractors. This Phase 1
  implementation searches equilibria and samples trajectories; it cannot globally
  exclude limit cycles, remote equilibria, or basin structure.
- The canonical framework is dimensional and site-resolved upstream. This workspace is
  nondimensional and treats `J` and `N` as independent sweep inputs; it does not derive
  them from codons, substitutions, sites, or expression.
- `U` and `A` are explicitly free binding-competent concentrations. They are not total
  pools, so the algebraic rapid-equilibrium resource formulas are internally consistent
  but cannot be compared directly to total-burden measurements.
- The fixed thresholds and broad parameter ranges are stress tests, not calibrated
  biological estimates. Point labels do not propagate measurement uncertainty.
- The two-state baseline has no regulation, division/dilution, compartmentation,
  toxicity feedback, or energetic cost. A stable subthreshold state is therefore only
  model-feasible, not a fitness or viability validation.
- `no_bounded_attracting_state` is deliberately an operational numerical label: no
  stable equilibrium found within the declared search plus escape of the declared
  initial conditions. It is not a theorem of nonexistence.
