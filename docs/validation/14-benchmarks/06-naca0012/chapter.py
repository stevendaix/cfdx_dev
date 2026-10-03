# %% [markdown]
# # NACA0012
#
# NACA0012 is a curved-wall external-aerodynamics benchmark. It exercises geometry, surface reconstruction, boundary layers, pressure distribution and, depending on the population, transition/turbulence modelling.
#
# ## Geometry
#
# The symmetric four-digit thickness law is
#
# [
# y_t=5tcleft[0.2969sqrt{x/c}-0.1260(x/c)-0.3516(x/c)^2+0.2843(x/c)^3-0.1015(x/c)^4ight],
# ]
#
# with (t=0.12). The exact trailing-edge convention must be frozen because the coefficient affects the geometry.
#
# ## Aerodynamic coefficients
#
# [
# C_L=rac{L}{	frac12ho U_infty^2c},
# qquad
# C_D=rac{D}{	frac12ho U_infty^2c},
# qquad
# C_M=rac{M}{	frac12ho U_infty^2c^2}.
# ]
#
# Also retain (C_p(x/c)), surface pressure and, where relevant, skin-friction/separation data.
#
# ## Frozen reference parameters
#
# Each population must explicitly state angle of attack, Reynolds number, Mach/compressibility treatment, turbulence/transition model, far-field location, wall treatment and force reference point.
#
# ## Mesh and extraction
#
# Use at least three systematically refined curved-wall meshes for a mesh-convergence claim. Record first-cell height, wall-normal growth, surface spacing, wake refinement and mesh-quality metrics. Pressure and viscous force contributions must be extracted independently.
#
# One coefficient on one mesh is insufficient qualification.
