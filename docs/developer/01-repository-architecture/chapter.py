# %% [markdown]
# # 01 — Repository Architecture
#
# The implementation should be read as several contracts rather than one undifferentiated tree.
#
# \[
# \boxed{\text{core}\rightarrow\text{numerics}\rightarrow\text{physics}\rightarrow\text{API}}
# \]
#
# alongside
#
# \[
# \boxed{\text{tests}\rightarrow\text{verification/validation evidence}}.
# \]
#
# ## Data boundaries
#
# CFDX distinguishes:
#
# - case setup: <case>.cfdx.h5;
# - solver state/checkpoint: <case>.dat.h5;
# - visualisation/output snapshots: <case>_<time>.vtu.
#
# A checkpoint is state, not a configuration language. This distinction is required for reproducible restart.
#
# ## Documentation boundary
#
# Theory explains equations and numerical choices. Developer documentation explains implementation contracts and extension points. V&V documentation explains evidence and acceptance.
#
# ## Extension rule
#
# A new numerical or physical capability should identify its owner layer, public contract, tests and documentation before merging.
