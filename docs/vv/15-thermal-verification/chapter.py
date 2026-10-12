# %% [markdown]
# # 15 — Thermal Verification
#
# Heat transfer verification covers conduction, convection, volumetric sources, radiation coupling and conjugate interfaces.
#
# ## Governing balance
#
# \[
# \rho c_p\frac{\partial T}{\partial t}
# +\nabla\cdot(\rho c_p\mathbf uT)
# =
# \nabla\cdot(k\nabla T)+S_T.
# \]
#
# ## Baselines
#
# Use 1-D conduction, transient diffusion, prescribed heat flux, convection boundaries and controlled source terms.
#
# ## CHT
#
# Verify perfect thermal contact and finite contact conductance independently. At an interface:
#
# \[
# q_1=q_2
# \]
#
# with the temperature jump determined by the declared contact resistance when present.
#
# ## Dimensionless checks
#
# Use Fourier, Peclet, Prandtl, Nusselt and Biot numbers to document the regime and ensure the benchmark is interpreted consistently.
