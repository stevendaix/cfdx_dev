# %% [markdown]
# # Plane Poiseuille flow
#
# Fully developed pressure-driven channel flow gives an analytical velocity profile and flow rate:
#
# \[
# \mu\frac{d^2u}{dy^2}=\frac{dp}{dx},
# \qquad
# u(y)=\frac{1}{2\mu}\frac{dp}{dx}y(y-H).
# \]
#
# The exact flow rate per unit span is
# \[
# Q=-\frac{H^3}{12\mu}\frac{dp}{dx}.
# \]
#
# ## Evidence
#
# Verify pressure-gradient/source treatment, diffusion, wall conditions, conservation and matrix/solver behaviour independently. Retain velocity profile error, flow-rate error, wall shear, mass imbalance, independently reconstructed fluxes and true residual.
#
# A systematic mesh campaign should contain at least three levels for observed-order claims. The historical N=8,16,32,64 sequence may be retained as provenance, but quantitative values must be regenerated from the current revision before being reused as current evidence.
#
# ## Scope
#
# This is primarily a scalar/viscous FV verification target. It must not be used as sole evidence that the complete SIMPLE/PISO/PIMPLE/COUPLED pressure-velocity architecture is qualified.
