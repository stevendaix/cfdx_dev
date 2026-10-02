# %% [markdown]
# 09 — AMG, MGR and Schur Complement Methods
#
# ## 9.1 Saddle-point structure
#
# A coupled incompressible system can be written
# \[
# A=
# \begin{bmatrix}
# A_{uu}&A_{up}\\
# A_{pu}&A_{pp}
# \end{bmatrix}.
# \]
# Eliminating the velocity block gives
# \[
# S=A_{pp}-A_{pu}A_{uu}^{-1}A_{up}.
# \]
# This exact Schur complement is generally expensive because applying it requires solving with \(A_{uu}\).
#
# ## 9.2 Approximate Schur
#
# Replace \(A_{uu}^{-1}\) by \(M_u^{-1}\):
# \[
# \widetilde S=A_{pp}-A_{pu}M_u^{-1}A_{up}.
# \]
# A useful verification is to compare \(\widetilde Sx\) against an exact/reference Schur action on controlled matrices.
#
# ## 9.3 Block factorisation
#
# When \(A_{uu}\) is invertible:
# \[
# A=
# \begin{bmatrix}I&0\\A_{pu}A_{uu}^{-1}&I\end{bmatrix}
# \begin{bmatrix}A_{uu}&0\\0&S\end{bmatrix}
# \begin{bmatrix}I&A_{uu}^{-1}A_{up}\\0&I\end{bmatrix}.
# \]
# This identity explains why block preconditioners can solve velocity and pressure subproblems separately.
#
# ## 9.4 AMG principle
#
# Algebraic multigrid constructs a hierarchy
# \[
# A_h,\;A_H,\;A_{H_2},\ldots
# \]
# and uses smoothing to reduce high-frequency error and coarse correction to reduce low-frequency error.
#
# ## 9.5 Galerkin coarse operator
#
# \[
# A_H=R A_h P.
# \]
# With \(R=P^T\) in a symmetric setting, the coarse operator inherits a variational relationship to the fine operator. Any non-Galerkin approximation is a different algorithm and needs its own verification.
#
# ## 9.6 V-cycle
#
# A conceptual V-cycle is
# \[
# x\leftarrow S_{\rm pre}(A,b,x)
# \rightarrow r=b-Ax
# \rightarrow r_H=Rr
# \rightarrow A_H e_H=r_H
# \rightarrow x\leftarrow x+Pe_H
# \rightarrow S_{\rm post}.
# \]
# Coarse solves, smoothing counts and transfer operators are part of the preconditioner definition.
#
# ## 9.7 Energy norm
#
# For SPD \(A\):
# \[
# \|e\|_A=\sqrt{e^TAe}.
# \]
# A V-cycle contraction diagnostic is
# \[
# q_E=\frac{\|e_{out}\|_A}{\|e_{in}\|_A}.
# \]
# Energy contraction is more structural than simply observing a residual decrease on one right-hand side.
#
# ## 9.8 Coarsening and interpolation
#
# Coarsening selects coarse variables \(C\) and fine variables \(F\). Interpolation approximates
# \[
# x_F\approx P_{FC}x_C.
# \]
# Ruge–Stüben-style strength of connection, classical interpolation and smoothed aggregation make different mathematical assumptions and must not be conflated.
#
# ## 9.9 MGR
#
# Multigrid reduction selects a subset of variables to retain while reducing other blocks. In block CFD systems the variable ordering and block elimination strategy are essential:
# \[
# x=
# \begin{bmatrix}x_C\\x_F\end{bmatrix}.
# \]
# The resulting hierarchy depends on the chosen coarse variables and relaxation.
#
# ## 9.10 Matrix update semantics
#
# If only values change:
# \[
# A^{new}=A^{old}+\Delta A
# \]
# while sparsity is unchanged, a hierarchy may need numerical refresh. If connectivity changes, \(P,R\), graph strength and coarse sparsity may all change, requiring structural rebuild. Silent reuse of stale hierarchy data is a correctness risk.
#
# ## 9.11 Verification
#
# The minimum hierarchy ladder is:
# exact Schur oracle → approximate Schur action → Galerkin identity → transfer dimensions → V-cycle contraction → matrix-update semantics → coupled solver evidence.
#
# ## 9.12 CFDX implementation
#
# The current source tree contains [amg_preconditioner.h](../../../src/cfdx/core/linalg/amg_preconditioner.h), [coupled_amg_schur.h](../../../src/cfdx/core/linalg/coupled_amg_schur.h), [block_schur.h](../../../src/cfdx/core/linalg/block_schur.h) and HYPRE AMG integration. Tests include [test_exact_schur.cpp](../../../tests/unit/test_exact_schur.cpp), [test_schur_preconditioner.cpp](../../../tests/unit/test_schur_preconditioner.cpp), [test_mgr_preconditioner.cpp](../../../tests/unit/test_mgr_preconditioner.cpp), [test_amg_preconditioner_qualification.cpp](../../../tests/validation/test_amg_preconditioner_qualification.cpp) and [test_coupled_block_schur_amg.cpp](../../../tests/validation/test_coupled_block_schur_amg.cpp).
#
# ## 9.13 Executable Schur oracle
# %%
import numpy as np
Auu=np.diag([2.,4.])
Aup=np.array([[1.],[2.]])
Apu=Aup.T
App=np.array([[0.]])
S=App-Apu@np.linalg.solve(Auu,Aup)
assert np.isfinite(S).all()
