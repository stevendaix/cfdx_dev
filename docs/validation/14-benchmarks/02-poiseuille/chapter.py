# %% [markdown]
# # Plane Poiseuille flow
#
# Fully developed pressure-driven channel flow gives an analytical profile and flow rate.
#
# \[
# \mu d^2u/dy^2=dp/dx,
# \qquad
# u(y)=\frac{1}{2\mu}\frac{dp}{dx}y(y-H).
# \]
#
# The exact flow rate per unit span is
#
# \[
# Q=-H^3(dp/dx)/(12\mu).
# \]
#
# Verify pressure-gradient/source treatment, diffusion, wall conditions, conservation and convergence.
#
# Campaign: systematic mesh refinement; velocity L2 error; flow-rate error; wall shear; mass imbalance; observed order. The N=8,16,32,64 structure may be retained as a campaign, but current values must be regenerated from the current revision.
