# %% [markdown]
# # 02 — Code Verification
#
# Code verification asks whether the implementation reproduces the specified mathematical algorithm or invariant.
#
# ## Verification ladder
#
# 1. geometry/topology invariants;
# 2. algebraic invariants;
# 3. constant and linear manufactured fields;
# 4. isolated operators;
# 5. MMS;
# 6. coupled manufactured problems.
#
# ## Example: face flux antisymmetry
#
# For an internal face shared by P and N,
#
# \[
# F_{P,f}+F_{N,f}=0.
# \]
#
# This verifies a discrete conservation contract independently of the global flow solution.
#
# ## Oracle types
#
# Use exact algebraic results, analytical solutions, manufactured solutions, independently implemented reference calculations, or invariant identities. A second copy of the same implementation is not an independent oracle.
#
# ## Reporting
#
# Record test ID, software revision, configuration, oracle, metric, tolerance and raw diagnostic. Distinguish numerical round-off from algorithmic failure.
