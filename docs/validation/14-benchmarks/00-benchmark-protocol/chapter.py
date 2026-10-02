# %% [markdown]
# # Benchmark protocol
#
# Every benchmark is a validation campaign with a declared population.
#
# Campaign sequence:
#
# reference definition -> CFDX setup -> mesh/solution verification -> physical comparison -> uncertainty assessment -> decision.
#
# Required record: exact revision; geometry; dimensions; properties; dimensionless parameters; boundary conditions; mesh; numerical schemes; solver criteria; monitored quantities; independent reference; uncertainty; metric; acceptance criterion; raw result; post-processing script.
#
# For a scalar quantity:
#
# \[
# e_Q=|Q_{CFDX}-Q_{ref}|/Q_{scale}.
# \]
#
# A status is evidence-backed: PASS means the pre-declared criterion is satisfied; FAIL means the executed case does not satisfy it; NOT QUALIFIED means evidence is insufficient; NOT RUN means no valid execution artifact exists.
#
# Criteria are fixed before interpretation and are never widened after seeing the CFDX result.
