"""
the rebuild's declared variables and parameters (gate G1.1).

every quantity is a concentration in uM, a rate in 1/s, or a pure number.
there is no "fraction of proteome" state: fractions are reported only as
explicit ratios such as A / P_T, with P_T itself a state sum.

imported by WO-02 onward so that one table is the single source of units.
"""
import units as u

# name -> (unit, meaning)
STATES = {
    "N": (u.UM, "native (folded, functional) protein, monomer units"),
    "U": (u.UM, "free non-native monomer: misfolded or error-bearing, chaperone-free"),
    "B": (u.UM, "chaperone-client complex; counts one client and one chaperone"),
    "A": (u.UM, "aggregated protein, in monomer-equivalents"),
    "C": (u.UM, "free chaperone (effective folding-competent pool)"),
}

PARAMS = {
    "s_P": (u.UM_PER_S, "total protein synthesis flux"),
    "eps": (u.ONE, "probability a newly synthesized chain enters U"),
    "s_C": (u.UM_PER_S, "chaperone synthesis flux"),
    "mu": (u.PER_S, "specific growth rate = dilution rate (ln2 / doubling time)"),
    "k_on": (u.PER_UM_PER_S, "chaperone-client association"),
    "k_off": (u.PER_S, "unproductive client release"),
    "k_cat": (u.PER_S, "completion of one chaperone cycle"),
    "phi": (u.ONE, "probability a completed cycle yields native protein"),
    "k_d": (u.PER_S, "degradation of free non-native monomer"),
    "k_a": (u.PER_UM_PER_S, "aggregation: monomer-equivalent flux k_a U^2"),
    "k_dis": (u.PER_S, "disaggregation of A back to U"),
    "k_dA": (u.PER_S, "degradation of aggregated protein"),
    "k_mis": (u.PER_S, "spontaneous unfolding of native protein"),
}

UNITS = {k: v[0] for k, v in {**STATES, **PARAMS}.items()}

# conservation laws the rebuild must satisfy (gate G1.3). written here as
# text; WO-02 proves them symbolically.
CONSERVATION = {
    "total protein P_T = N + U + B + A":
        "dP_T/dt = s_P - k_d*U - k_dA*A - mu*P_T  "
        "(source: synthesis; sinks: degradation of U and A, dilution)",
    "total chaperone C_T = C + B":
        "dC_T/dt = s_C - mu*C_T  "
        "(source: chaperone synthesis; sink: dilution; chaperone not degraded)",
}
