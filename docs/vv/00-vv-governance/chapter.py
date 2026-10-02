# %% [markdown]
# # 00 — V&V Governance
#
# Verification and validation are evidence disciplines, not labels applied to a green test.
#
# ## Claim model
#
# A quantitative claim is represented by
#
# \[
# C=(R,O,M,A,P,E),
# \]
#
# where R=requirement, O=oracle/reference, M=metric, A=acceptance criterion, P=declared population and E=retained evidence.
#
# ## Maturity
#
# **Implemented** means the capability exists in source. **Verified** means a defined mathematical/software property has evidence. **Validated** means comparison with an independent physical/reference basis has been performed. **Qualified** means the declared population and all gates are complete.
#
# These states are not interchangeable.
#
# ## Evidence hierarchy
#
# Cheap invariants should precede isolated operators, MMS, convergence campaigns, coupled benchmarks and qualification. A coupled pass cannot compensate for a failed lower-level invariant.
#
# ## Acceptance criteria
#
# Criteria are defined before interpreting results:
#
# \[
# |Q-Q_{ref}|\le\epsilon_Q,\qquad
# p_{obs}\ge p_{min},\qquad
# |R_C|/Q_{ref}\le\epsilon_C.
# \]
#
# Criteria must have a documented numerical or engineering basis. They must not be widened after observing the result.
#
# ## Evidence package
#
# Retain revision, case/setup, mesh, physics, numerical schemes, solver criteria, raw outputs, post-processing definition, environment and acceptance decision.
#
# ## Status rule
#
# Missing evidence keeps the claim open. A green CI result is evidence for the jobs that ran; it is not automatically validation or qualification.
