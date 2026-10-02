# %% [markdown]
# # Thermal benchmark family
#
# ## One-dimensional conduction
#
# \[
# d^2T/dx^2=0,
# \qquad
# T(x)=T_0+(T_L-T_0)x/L,
# \]
#
# and
#
# \[
# q_x=-k(T_L-T_0)/L.
# \]
#
# This isolates thermal diffusion and temperature boundary conditions.
#
# ## Conjugate heat transfer
#
# Perfect contact requires
#
# \[
# q_s=q_f,\qquad T_s=T_f.
# \]
#
# With contact resistance Rc,
#
# \[
# q=(T_1-T_2)/R_c.
# \]
#
# Verify temperature error, heat-flux error, interface conservation and coupling iteration error. Thermal validation is separated into conduction, convection and CHT populations rather than treated as one generic pass.
