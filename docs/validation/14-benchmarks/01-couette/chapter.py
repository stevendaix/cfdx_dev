# %% [markdown]
# # Couette flow
#
# Plane Couette flow isolates viscous momentum transport. Two infinite plates are separated by (H); wall velocities are (U_0) and (U_H); the flow is steady, incompressible and has zero streamwise pressure gradient.
#
# [
# murac{d^2u}{dy^2}=0,
# qquad
# u(y)=U_0+(U_H-U_0)y/H.
# ]
#
# The analytical wall shear is
# [
# 	au_{xy}=mu(U_H-U_0)/H.
# ]
#
# The benchmark should explicitly verify velocity boundary-condition signs, gradient accuracy, viscous fluxes, pressure neutrality, continuity and wall-force extraction.
#
# ## Evidence
#
# Retain L1/L2/Linf velocity errors, wall-shear error, mass imbalance, independently reassembled momentum residual, nonlinear/linear iteration history and mesh-quality metrics. Use at least three systematically related meshes when an observed order is claimed.
#
# The iterative error must be demonstrably smaller than the spatial discretisation signal before interpreting mesh convergence.
#
# ## Scope
#
# Couette validates viscous transport and the selected solver path for this configuration. It does not by itself qualify convection schemes, turbulence models, difficult meshes or all pressure-velocity algorithms.
