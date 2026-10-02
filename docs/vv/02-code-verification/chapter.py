# %% [markdown]
# # 02 — Code Verification
#
# Code verification asks whether the implementation reproduces a specified algorithm or invariant.
#
# ## Conservation example
#
# For an internal face shared by owner P and neighbour N,
#
# \[
# F_{P,f}+F_{N,f}=0.
# \]
#
# The test should evaluate the same physical face flux through both orientations and assert antisymmetry within a declared floating-point tolerance.
#
# ## Verification hierarchy
#
# 1. algebraic invariants;
# 2. constant and linear manufactured fields;
# 3. isolated operators;
# 4. MMS;
# 5. coupled manufactured problems.
#
# A coupled benchmark does not replace lower-level operator evidence because compensating errors can hide each other.
#
# ## Evidence
#
# Record test identifier, source revision, numerical configuration, measured quantity and acceptance criterion.
