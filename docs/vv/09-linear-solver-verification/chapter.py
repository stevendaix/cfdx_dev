# %% [markdown]
# # 09 — Linear Solver Verification
#
# For
#
# \[
# Ax=b,
# \]
#
# verify the solver independently of the CFD application.
#
# ## Metrics
#
# True residual:
#
# \[
# r=b-Ax.
# \]
#
# Relative residual:
#
# \[
# \eta=\frac{\|r\|}{\|b\|}.
# \]
#
# Backward error may be more meaningful for ill-conditioned systems.
#
# ## Test ladder
#
# Exact diagonal systems → small dense systems → sparse systems → nonsymmetric systems → ill-conditioned systems → null-space systems → CFD matrices.
#
# ## Krylov methods
#
# Verify CG assumptions separately from GMRES/FGMRES/BiCGStab. Variable preconditioning must not be tested with an incompatible solver contract.
#
# ## AMG/MGR/Schur
#
# Test the preconditioner against exact small-system oracles before using large CFD cases. Energy contraction is distinct from residual reduction and should be used where the method's contract requires it.
