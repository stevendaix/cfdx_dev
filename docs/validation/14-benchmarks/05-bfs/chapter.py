# %% [markdown]
# # Backward-facing step
#
# The backward-facing step is a separated-flow benchmark. Its value is not a single reattachment-length number: the separated shear layer, pressure recovery, wall shear and velocity field must be considered together.
#
# ## Frozen definition
#
# Record the exact step height \(h\), channel dimensions, inlet profile, bulk/centreline velocity definition, viscosity, Reynolds-number convention, outlet condition and wall treatment.
#
# \[
# Re=\frac{U_bh}{\nu}
# \]
# is only meaningful when \(U_b\) and \(h\) follow the declared reference convention.
#
# ## Quantities
#
# Primary quantities may include:
#
# - reattachment length \(x_r/h\);
# - wall-pressure distribution;
# - wall-shear sign and separation/reattachment locations;
# - velocity profiles at prescribed stations;
# - mass conservation;
# - turbulence quantities when a turbulent population is selected.
#
# Reattachment length should be extracted from the sign change of wall shear or an explicitly frozen reference definition, not inferred visually from a contour.
#
# ## Numerical sensitivity
#
# The campaign must record step-region mesh resolution, streamwise expansion ratio, near-wall treatment, convection scheme, pressure-velocity coupling and inlet profile. Coarse meshes can produce apparently plausible recirculation while moving \(x_r/h\) materially.
#
# ## Qualification boundary
#
# Reattachment-length agreement alone does not qualify a turbulence model or pressure-velocity algorithm. Field, balance and convergence evidence remain required.
