# %% [markdown]
# # 16 — Radiation Verification
#
# Radiation is nonlinear and must be verified independently before multiphysics coupling.
#
# ## Blackbody oracle
#
# \[
# E_b=\sigma T^4.
# \]
#
# ## Surface exchange
#
# Verify view-factor reciprocity and closure where applicable:
#
# \[
# A_iF_{ij}=A_jF_{ji},
# \qquad
# \sum_jF_{ij}=1
# \]
#
# for a closed enclosure under the assumptions of the method.
#
# ## Models
#
# Verify S2S, P1 and DOM/RTE components independently where implemented. Distinguish geometric/view-factor error from nonlinear radiation iteration error.
#
# ## Coupling
#
# A coupled thermal-radiation pass requires independent evidence for both radiation and energy conservation.
