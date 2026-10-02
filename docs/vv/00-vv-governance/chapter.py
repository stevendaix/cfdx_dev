# %% [markdown]
# # 00 — V&V Governance
#
# A V&V claim is not simply “test passes”. Define
#
# \[
# C=(R,O,M,A,P,E),
# \]
#
# where R is the requirement, O the oracle/reference, M the metric, A the acceptance criterion, P the declared population and E the retained evidence.
#
# ## Evidence levels
#
# | Level | Meaning |
# |---|---|
# | Code verification | implementation checked against a mathematical/software oracle |
# | Solution verification | numerical error and convergence controlled |
# | Validation | independent physical/reference comparison |
# | Qualification | complete declared population satisfies its gates |
#
# ## Acceptance criteria
#
# Criteria must be defined before interpreting the result:
#
# \[
# |Q-Q_{ref}|\le\epsilon_Q,\qquad
# p_{obs}\ge p_{min},\qquad
# \epsilon_C\le\epsilon_{C,max}.
# \]
#
# Thresholds must come from the requirement, numerical analysis or a documented engineering basis; they must not be selected after seeing the result.
#
# ## Evidence package
#
# Retain the inputs required to reproduce the claim, software revision, execution configuration, raw result, post-processing method and derived metric.
#
# Incomplete evidence remains incomplete; adjacent passing tests do not promote it automatically.
