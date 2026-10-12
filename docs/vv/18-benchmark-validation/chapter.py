# %% [markdown]
# # 18 — Benchmark Validation
#
# Validation compares CFDX with an independent physical or reference basis. It is different from code verification.
#
# ## Benchmark record
#
# Define geometry, operating conditions, fluid properties, boundary conditions, numerical method, quantity of interest, reference source, uncertainty and acceptance criterion.
#
# ## Comparison
#
# For a quantity Q:
#
# \[
# e_Q=\frac{|Q_{CFDX}-Q_{ref}|}{Q_{ref}}.
# \]
#
# The metric must be appropriate to the reference uncertainty; pointwise comparison is not always meaningful.
#
# ## Existing CFDX campaigns
#
# Couette, Poiseuille, Ghia, thermal/radiation and the more difficult VMFL036/BFS/NACA families must be reported from actual current repository evidence. Do not infer current status from historical discussion.
#
# ## Independence
#
# A benchmark source, analytical reference or independent code should be identified. A value copied from CFDX documentation is not an independent validation oracle.
