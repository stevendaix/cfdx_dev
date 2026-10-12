# %% [markdown]
# # 00 — Development Philosophy
#
# CFDX is scientific software. Code quality and numerical correctness are coupled: an implementation without a mathematical contract is incomplete, while a passing benchmark without traceability is weak evidence.
#
# ## Implementation maturity
#
# - **Implemented:** the code path exists.
# - **Verified:** a defined software or mathematical property has executable evidence.
# - **Validated:** independent physical/reference evidence supports the declared case.
# - **Qualified:** the declared capability population and acceptance gates are complete.
#
# Documentation never changes numerical maturity.
#
# ## Change discipline
#
# A numerical change should close a coherent contract:
#
# \[
# \text{requirement}\rightarrow\text{equations}\rightarrow\text{source}\rightarrow\text{tests}\rightarrow\text{evidence}.
# \]
#
# Tolerances are part of the numerical contract. They must not be relaxed merely to obtain a green result.
#
# ## Failure handling
#
# Classify failures before changing code: modelling, mesh/geometry, discretisation, boundary condition, linear solve, nonlinear convergence, post-processing, or infrastructure.
#
# Do not silently replace the requested algorithm with another algorithm.
#
# ## Reproducibility
#
# Numerical claims should retain software revision, setup, mesh, solver configuration and metric definition. Record hardware or parallel configuration when it can affect the result.
#
# Theory is documented under docs/theory/ and V&V under docs/vv/.
