# %% [markdown]
# # N9-S8 — NACA 0012 at Re=1000
#
# This executable Jupytext chapter is the canonical scientific presentation
# and literature-comparison study for N9-S8. The population is **laminar,
# incompressible and steady first**: this is a pressure--velocity robustness
# study, not a turbulence-model validation case.
#
# %% [markdown]
# ## 1. Frozen physical problem
#
# - NACA 0012, chord (c=1)
# - (U_\infty=1), (ho=1)
# - (Re_c=1000), hence (
u=10^{-3})
# - (alpha=0^\circ)
# - incompressible Navier--Stokes, laminar
# - no-slip airfoil wall
# - uniform far-field velocity
#
# [
# Re_c = \frac{U_\infty c}{\nu}=1000.
# ]
#
# Aerodynamic coefficients:
#
# [
# C_D=\frac{D}{\tfrac12\rho U_\infty^2c},quad
# C_L=\frac{L}{\tfrac12\rho U_\infty^2c},quad
# C_M=\frac{M}{\tfrac12\rho U_\infty^2c^2}.
# ]
#
# At zero incidence, the symmetric geometry provides an independent sanity
# gate: converged lift must approach zero without imposing (C_L=0).
#
# %% [markdown]
# ## 2. Independent mesh population
#
# NASA/TMR NACA 0012 C-grid family:
#
# | level | grid | role |
# |---|---:|---|
# | coarse | 225 x 65 | debugging/refinement |
# | medium | 449 x 129 | primary result |
# | fine | 897 x 257 | mesh convergence |
#
# These are external validation assets, not CFDX-generated meshes. Materialization
# is controlled and offline; CI must never silently download them.
#
# %% [markdown]
# ## 3. N9 algorithm matrix and retained evidence
#
# All six N9 pressure--velocity families use the same mesh, physics, spatial
# discretization and convergence contract. Retain nonlinear history, true
# residuals, linear iterations, continuity, force balance, (C_D,C_L,C_M),
# pressure/viscous drag, (C_p(x/c)), (C_f(x/c)) where available, runtime
# and iteration count.
#
# %% [markdown]
# ## 4. Literature reference hierarchy
#
# Kurtulus (2015), *On the Unsteady Behavior of the Flow Around NACA 0012 Airfoil
# with Steady External Conditions at Re=1000*, International Journal of Micro Air
# Vehicles 7(3), 301--326, DOI 10.1260/1756-8293.7.3.301, is the primary
# independent numerical reference. At (alpha=0^\circ), it reports approximately
# (C_D=0.12); symmetry supplies (C_L=0).
#
# This is a **numerical literature reference**, not an analytic solution and not
# a hard-coded CFDX tolerance. Kunz & Kroo (2000) and Kouser et al. (2021) are
# independent low-Re cross-checks. Swanson & Langer (2016) is provenance only
# because its cited cases are compressible.
#
# %% [markdown]
# ## 5. Quantitative comparison
#
# [
# E_{C_D}=\frac{|C_{D,CFDX}-0.12|}{0.12},
# qquad E_{C_L}=|C_{L,CFDX}|.
# ]
#
# Where reference distributions exist, compare (C_p(x/c)) and (C_f(x/c))
# with documented interpolation and surface-coordinate conventions. Literature
# discrepancy never replaces convergence, conservation, symmetry or refinement.
#
# %%
from pathlib import Path
import json

RESULT = Path("results/n9-s8/comparison.json")
if RESULT.exists():
    report = json.loads(RESULT.read_text(encoding="utf-8"))
    for row in report["results"]:
        print(
            f"{row['case']:24s} {row['qoi']:8s} "
            f"computed={row['computed']:.8g} reference={row['reference']:.8g} "
            f"abs={row['absolute_error']:.3e} rel={row['relative_error']}"
        )
else:
    print("No retained Total Validation comparison artifact is available yet.")
    print("N9-S8 remains NOT QUALIFIED until the controlled campaign executes.")

# %% [markdown]
# ## 6. Verification and qualification gates
#
# Required evidence: nonlinear convergence, independent conservation, zero-lift
# symmetry without constraining (C_L), stable integrated and distributed QoIs
# under the three-level refinement, reproducibility across all six N9 methods,
# and documented comparison with independent low-Re literature.
#
# **Status: NOT QUALIFIED** until the controlled Total Validation campaign has
# executed all three meshes and all six N9 methods with retained evidence.
