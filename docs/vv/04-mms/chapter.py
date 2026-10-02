# %% [markdown]
# # 04 — Manufactured Solutions
#
# Given
#
# \[
# L(u)=f,
# \]
#
# select a smooth exact field u_e and construct
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
# \frac{\sum_PV_P|u_P-u_e|^2}
# {\sum_PV_P|u_e|^2}
# \right)^{1/2}.
# \]
#
# ## Campaign
#
# 1. choose a smooth solution;
# 2. derive the source independently;
# 3. use multiple mesh levels;
# 4. control iterative error;
# 5. compute observed order;
# 6. inspect local error as well as global norms.
#
# MMS is useful for gradients, diffusion, convection, source linearisation and coupled operators.
#
# A wrong order is evidence to investigate, not a reason to change the criterion after execution.
