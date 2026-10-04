# %% [markdown]
# # 10 — Pressure–Velocity Verification
#
# In incompressible flow the algebraic system has a saddle-point structure:
#
# \[
# \begin{bmatrix}A&G\\D&0\end{bmatrix}
# \begin{bmatrix}u\\p\end{bmatrix}
# =
# \begin{bmatrix}b_u\\b_p\end{bmatrix}.
# \]
#
# ## Required checks
#
# - pressure gauge/null space;
# - compatibility of mass source and boundary fluxes;
# - divergence reduction;
# - pressure-correction consistency;
# - face-flux consistency;
# - convergence independent of residual-only claims.
#
# ## Algorithms
#
# Verify SIMPLE/SIMPLEC, PISO, PIMPLE, fractional-step projection and coupled methods with controlled problems before comparing full CFD benchmarks.
#
# ## Schur oracle
#
# For small block systems compare the implemented pressure operator with an independently constructed Schur complement:
#
# \[
# S=D A^{-1}G.
# \]
#
# Approximate Schur operators must document the approximation and its qualification scope.
