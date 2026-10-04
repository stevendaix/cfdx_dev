# %% [markdown]
# # 11 — Gradient Verification
#
# The gradient contract includes both the cell gradient and the independent face-value reconstruction.
#
# ## Green–Gauss
#
# \[
# \nabla\phi_P\approx\frac{1}{V_P}\sum_f\phi_f\mathbf S_f.
# \]
#
# ## Least squares
#
# \[
# (A^TWA)g=A^TWb.
# \]
#
# Verify rank and conditioning of \(A^TWA\); do not hide rank deficiency with arbitrary regularisation.
#
# ## Required population
#
# Test constant and linear fields, polyhedral and tetrahedral meshes, boundary-neighbour policies, skewness, non-orthogonality, conditioning and refinement order.
#
# ## Current qualification rule
#
# Formal second order must be demonstrated by an actual CFDX refinement campaign. A reference implementation or a theoretical order statement is not qualification evidence.
#
# Gradient limiting and face reconstruction must have separate tests so that one cannot mask an error in the other.
