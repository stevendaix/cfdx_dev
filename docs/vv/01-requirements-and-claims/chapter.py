# %% [markdown]
# # 01 — Requirements and Claims
#
# V&V begins by making claims testable.
#
# ## Requirement record
#
# \[
# R=(ID,statement,scope,metric,criterion,evidence).
# \]
#
# The scope must identify the equation, physics, numerical method, mesh family and parameter population for which the claim is made.
#
# ## Traceability
#
# \[
# requirement\rightarrow mathematical\ property\rightarrow implementation\rightarrow test\rightarrow artifact\rightarrow decision.
# \]
#
# Every qualification row should therefore point to a stable requirement, a reproducible campaign and retained results.
#
# ## Good and weak claims
#
# “Gradient is second order” is incomplete unless the mesh family, norm, refinement sequence, boundary treatment and asymptotic criterion are specified. “Test passes” is incomplete unless the property and acceptance threshold are known.
#
# ## Change control
#
# A change to equations, discretisation, solver, tolerance or boundary treatment can invalidate existing evidence. Reassess affected claims rather than assuming regression tests are sufficient.
