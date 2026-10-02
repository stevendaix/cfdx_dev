# %% [markdown]
# # 09 — Linear Solver Verification
#
# true residual, exact systems, null spaces and preconditioner evidence.
#
# ## Development / evidence contract
#
# The chapter must connect the mathematical requirement to an implementation point and to executable evidence.
#
# \[
# \text{requirement}\rightarrow\text{implementation}\rightarrow\text{test}\rightarrow\text{evidence}.
# \]
#
# A passing test is evidence only for the property and population it explicitly defines. It does not automatically establish validation or qualification.
#
# ## Integrity rules
#
# Do not silently substitute algorithms, relax numerical criteria to obtain a pass, or present a reference value as a CFDX result.
