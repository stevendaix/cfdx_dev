# %% [markdown]
# # 04 — Method of Manufactured Solutions
#
# For a differential operator
#
# \[
# L(u)=f,
# \]
#
# choose a smooth exact field \(u_e\) and construct
#
# \[
# f_{MMS}=L(u_e).
# \]
#
# The numerical error can then be measured directly:
#
# \[
# E_2=
# \left(
# \frac{\sum_P V_P|u_P-u_e|^2}
# {\sum_P V_P|u_e|^2}
# \right)^{1/2}.
# \]
#
# ## Campaign
#
# Select a smooth field, derive the source independently, impose compatible boundaries, refine the mesh, control iterative error and calculate observed order.
#
# ## Scope
#
# MMS should be used for gradients, diffusion, convection, source terms, temporal operators and coupled equations. Use multiple mesh families where geometry quality may influence the result.
#
# ## Diagnostics
#
# Inspect local error, boundary error, conditioning and conservation in addition to the global norm. A global norm can hide local defects.
