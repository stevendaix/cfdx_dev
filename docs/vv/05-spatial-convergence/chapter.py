# %% [markdown]
# # 05 — Spatial Convergence
#
# Assume
#
# \[
# E_h=Ch^p.
# \]
#
# For two levels,
#
# \[
# p_{obs}=\frac{\ln(E_H/E_h)}{\ln(H/h)}.
# \]
#
# ## Campaign requirements
#
# Use at least three levels when practical. Record the characteristic mesh-size definition, cell family, geometry perturbations and all numerical settings.
#
# ## Asymptotic regime
#
# Formal order should only be claimed when the sequence shows a defensible asymptotic trend. A single slope is insufficient when boundary reconstruction, geometry or conditioning can dominate.
#
# ## N2 application
#
# Gradient verification must distinguish polyhedral and tetrahedral meshes, boundary-neighbour reconstruction, conditioning and the separation between gradient calculation and face-value reconstruction.
