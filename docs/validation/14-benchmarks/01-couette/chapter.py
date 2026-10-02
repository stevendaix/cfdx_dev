# %% [markdown]
# # Couette flow
#
# Plane Couette flow isolates viscous momentum transport. Two infinite plates are separated by H; wall velocities are U0 and UH; the flow is steady, incompressible and has zero streamwise pressure gradient.
#
# \[
# \mu d^2u/dy^2=0,
# \qquad
# u(y)=U_0+(U_H-U_0)y/H.
# \]
#
# The analytical wall shear is
#
# \[
# \tau_{xy}=\mu(U_H-U_0)/H.
# \]
#
# This case verifies velocity Dirichlet treatment, gradients, viscous fluxes, pressure neutrality and conservation.
#
# Campaign: at least three systematic meshes; L2 velocity error; wall-shear error; observed order; mass balance; iterative error smaller than the spatial signal.
#
# Acceptance: pre-declare expected order and tolerances. A solver convergence flag alone is not acceptance evidence.
