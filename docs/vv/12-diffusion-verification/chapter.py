# %% [markdown]
# # 12 — Diffusion Verification
#
# Diffusion begins from
#
# \[
# \nabla\cdot(\Gamma\nabla\phi).
# \]
#
# ## Baseline tests
#
# Constant field → zero flux. Linear field on an appropriate mesh → exact/controlled gradient. Smooth MMS → observed order.
#
# ## Non-orthogonality
#
# Separate the orthogonal contribution from the correction:
#
# \[
# \mathbf d_f\cdot\nabla\phi
# =
# \text{orthogonal part}+\text{non\!\!-orthogonal correction}.
# \]
#
# Test each contribution and their sum.
#
# ## Boundaries
#
# Verify Dirichlet, Neumann and mixed/Robin conditions independently, including flux sign conventions and source linearisation.
#
# ## Acceptance
#
# Report flux error, field error, conservation defect and observed order. Do not qualify diffusion from a single benchmark.
