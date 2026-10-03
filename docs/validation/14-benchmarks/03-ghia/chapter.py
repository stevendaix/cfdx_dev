# %% [markdown]
# # Ghia lid-driven cavity
#
# A square cavity has stationary walls except for a moving lid of velocity (U). The flow is incompressible and closed, so net mass flux is zero and pressure has a constant null space.
#
# [
# Re=UH/
u.
# ]
#
# The classical Ghia, Ghia and Shin reference is used through centreline velocity profiles and vortex information.
#
# ## Reference comparison
#
# Freeze the exact Reynolds number, cavity size, lid velocity, mesh sequence, convection scheme, pressure-velocity algorithm, interpolation used to extract centreline values and reference-data table revision.
#
# Quantities include:
#
# - vertical-centreline horizontal velocity;
# - horizontal-centreline vertical velocity;
# - primary-vortex position/strength;
# - secondary vortices where resolved;
# - mass imbalance;
# - independently reassembled momentum residual.
#
# Profile errors should be retained as both pointwise and norm-based metrics. Vortex locations are supplementary because they can be more mesh-sensitive than centreline velocities.
#
# ## Scope
#
# Passing one Reynolds number does not qualify the entire cavity population. The evidence should remain separated for Re=100, Re=400, Re=1000 and any future cases.
#
# A residual-only convergence flag is insufficient for reference validation.
