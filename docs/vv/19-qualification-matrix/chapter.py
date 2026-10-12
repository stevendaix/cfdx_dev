# %% [markdown]
# # 19 — Qualification Matrix
#
# Qualification is a bounded population statement, not a global “CFDX is validated” label.
#
# ## Matrix
#
# Each row should contain:
#
# | Field | Meaning |
# |---|---|
# | Capability | equation/model/algorithm |
# | Population | geometry/regime/mesh/method range |
# | Requirement | measurable claim |
# | Test | reproducible campaign |
# | Criterion | pre-declared gate |
# | Evidence | retained artifact |
# | Status | implemented/verified/validated/qualified |
#
# ## Gate logic
#
# A capability is qualified only when every required population member and gate is satisfied. Optional or failed cases remain visible and must not be silently removed from the population.
#
# ## N2 example
#
# Gradient methods require separate evidence for Green–Gauss, LS/WLS, vertex methods, boundary reconstruction, conditioning and observed order. A single smooth benchmark cannot qualify the whole gradient capability.
