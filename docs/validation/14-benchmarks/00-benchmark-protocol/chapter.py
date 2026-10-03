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

# %% [markdown]
# ## Cross-benchmark evidence contract
#
# Every campaign record should be sufficient for an independent reviewer to reconstruct what was actually tested.
#
# ### A. Frozen case definition
#
# Record software revision/build configuration; exact geometry and dimensionality; fluid/material properties and units; governing dimensionless groups; boundary and initial conditions; mesh source/revision/hash and quality metrics; spatial/temporal schemes; gradient/reconstruction and limiter; pressure-velocity algorithm; linear solver/preconditioner; relaxation/CFL/time-step policy; and convergence/stopping criteria.
#
# A benchmark name is not a sufficient case definition.
#
# ### B. Independent observables
#
# Separate solver convergence, conservation/balance, field error or reference-QoI error, mesh/time/model sensitivity, and runtime/performance diagnostics.
#
# A physical QoI must never be the only convergence indicator.
#
# ### C. Reference provenance
#
# Each oracle must be labelled analytical, manufactured, discrete-exact, published experimental/reference data, published CFD/reference-code result, or internal regression baseline.
#
# Published CFD values are comparison data, not proof of CFDX correctness. If a reference code is used, preserve its model, mesh assumptions, physical parameters and extraction definition.
#
# ### D. Error definitions
#
# \[
# e_Q=\frac{|Q_{CFDX}-Q_{ref}|}{Q_{scale}},
# \qquad
# \Delta Q=Q_{CFDX}-Q_{ref}.
# \]
#
# The scale must be declared and the raw CFDX quantity retained.
#
# For profiles, retain pointwise errors and norms where appropriate:
# \[
# L_2=\left(\frac1N\sum_i e_i^2\right)^{1/2},
# \qquad
# L_\infty=\max_i|e_i|.
# \]
#
# ### E. Refinement requirements
#
# A steady spatial campaign should normally contain at least three systematically related meshes when observed order is claimed. Temporal order requires at least three relevant time-step levels. Radiation transport claims additionally require angular/model-space refinement when that discretisation controls the error.
#
# Coarse/medium/fine must identify actual mesh/time-step identifiers, not labels alone.
#
# ### F. Force and flux extraction
#
# Force, heat-flux and radiation QoIs must be extracted from CFDX fields independently of any hard-coded expected result. Record pressure and viscous contributions separately whenever both exist.
#
# \[
# \mathbf F=\int_S(-p\mathbf n+\boldsymbol\tau\mathbf n)\,dA.
# \]
#
# The projected coefficient definition, reference area and sign convention must be frozen with the benchmark.
#
# ### G. Qualification boundary
#
# \`EXECUTED\` means the solver ran and produced an artifact; it is not a physical PASS.
#
# \`PASS\` requires the retained artifact to satisfy all declared gates for that benchmark population.
#
# \`QUALIFIED\` is a higher-level decision over the declared population. A single successful mesh or parameter point cannot qualify a broader model family.
#
# ### H. Next campaign priorities
#
# The benchmark suite should be completed in dependency order:
#
# **Couette/Poiseuille/Ghia resynchronisation → VMFL036 → NACA0012/BFS → thermal/CHT → radiation → broader turbulence/external-flow populations.**
#
# This follows the qualification authority in Issue #118 and the numerical-method dependencies in Issue #461. Open implementation blockers remain visible rather than becoming documentation PASS states.
