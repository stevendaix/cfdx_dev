# %% [markdown]
# # NACA0012
#
# NACA0012 is a curved-wall external-aerodynamics benchmark. It exercises geometry, surface gradients, boundary layers, pressure distribution and, when used, turbulence/transition modelling.
#
# The symmetric NACA four-digit thickness law is
#
# \[
# y_t=5tc[0.2969\sqrt{x/c}-0.1260(x/c)-0.3516(x/c)^2+0.2843(x/c)^3-0.1015(x/c)^4].
# \]
#
# For NACA0012, t=0.12.
#
# Aerodynamic coefficients:
#
# \[
# C_L=L/(0.5\rho U_\infty^2c),\quad
# C_D=D/(0.5\rho U_\infty^2c).
# \]
#
# Depending on the reference, also compare C_M, surface C_p and separation.
#
# Exact angle of attack, Reynolds number, Mach/compressibility treatment, transition/turbulence model and far-field boundary conditions must match the reference.
#
# Curved-wall mesh convergence is mandatory. One coefficient on one mesh is insufficient qualification.
