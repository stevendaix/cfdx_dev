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
# \\[
# e_Q=|Q_{CFDX}-Q_{ref}|/Q_{scale}.
# \\]
#
# A status is evidence-backed: PASS means the pre-declared criterion is satisfied; FAIL means the executed case does not satisfy it; NOT QUALIFIED means evidence is insufficient; NOT RUN means no valid execution artifact exists.
#
# Criteria are fixed before interpretation and are never widened after seeing the CFDX result.
#
# ## Executable campaign map
#
# The repository already contains executable entry points for the canonical benchmark
# families. run_campaign.py is the common orchestration layer and writes one log plus
# one machine-readable results.json record per campaign.
#
# | Case | Short regression | Full campaign entry point | Evidence meaning |
# |---|---|---|---|
# | Couette | test_couette_quick | test_phase9_acceptance | coupled pressure-velocity regression; quantitative profile/QoI evidence required for qualification |
# | Poiseuille | test_poiseuille_quick | test_poiseuille_diagnostics | solver result plus refinement/order and flow-rate evidence |
# | Ghia | test_ghia_cavity_quick | test_ghia_cavity | centreline comparison plus continuity and mesh/reference evidence |
# | VMFL036 | test_vmfl036_axisymmetric | test_vmfl036_axisymmetric | physical axisymmetric sphere drag campaign |
# | BFS | test_bfs_quick | test_bfs_qualification | quick test is only a contract smoke; full campaign needs coarse/medium/fine meshes |
# | NACA0012 | test_naca0012_quick | test_naca0012_qualification | quick test is only geometry/BC coverage; force convergence is separate |
# | Thermal | test_thermal_vv, test_cht_validation | same executables | component/solver V&V; full CHT qualification remains distinct |
# | Radiation | test_radiation_vv, test_s2s_radiation_vv | same executables | radiation model V&V; coupled thermal-radiation qualification remains distinct |
#
# The distinction between quick and full is intentional. A successful contract test
# never becomes a qualification PASS by implication.
#
# ## Evidence status rules
#
# - VERIFIED_CONTRACT: the short executable completed successfully; this proves only
#   the declared contract covered by that test.
# - EXECUTED: the full executable completed successfully; the output still requires
#   quantitative review against the declared oracle and acceptance criterion.
# - PASS: may be assigned only after the retained result satisfies the pre-declared
#   physical/QoI/conservation/refinement criterion.
# - FAIL: an executed campaign violates its criterion or exits non-zero.
# - NOT_RUN: the required executable or refinement meshes are unavailable.
#
# run_campaign.py deliberately never emits PASS or QUALIFIED. Qualification is
# a V&V decision based on the retained evidence, not a process-exit-code alias.
#