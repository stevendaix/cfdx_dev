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
# ## Requirements
#
# Define the characteristic \(h\), mesh family, refinement ratio, norm, quantity of interest and acceptance order before execution.
#
# Use at least three levels where practical. Formal order should only be claimed when the sequence is plausibly asymptotic.
#
# ## Mesh-family controls
#
# Separate effects of cell size, skewness, non-orthogonality, boundary reconstruction and topology. Tetrahedral, hexahedral and general polyhedral sequences are not interchangeable evidence.
#
# ## Richardson/GCI
#
# For a sufficiently regular sequence, Richardson extrapolation may estimate the limiting value. GCI should be reported with its assumptions and should not be used to manufacture an asymptotic regime.
