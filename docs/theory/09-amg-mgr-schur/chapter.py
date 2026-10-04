# %% [markdown]
# 09 — AMG, MGR and Schur Complement Methods
#
# ## 9.1 Saddle-point structure
#
# \[
# A=
# \begin{bmatrix}
# A_{uu}&A_{up}\\
# A_{pu}&A_{pp}
# \end{bmatrix},
# \qquad
# S=A_{pp}-A_{pu}A_{uu}^{-1}A_{up}.
# \]
# The exact Schur action is expensive because applying it requires an \(A_{uu}\) solve.
#
# ## 9.2 Block factorisation
#
# \[
# A=
# \begin{bmatrix}I&0\\A_{pu}A_{uu}^{-1}&I\end{bmatrix}
# \begin{bmatrix}A_{uu}&0\\0&S\end{bmatrix}
# \begin{bmatrix}I&A_{uu}^{-1}A_{up}\\0&I\end{bmatrix}.
# \]
# This identity is the algebraic foundation of block-Schur preconditioning.
#
# ## 9.3 Approximate Schur and LSC/BFBt concepts
#
# With \(M_u^{-1}\approx A_{uu}^{-1}\),
# \[
# \widetilde S=A_{pp}-A_{pu}M_u^{-1}A_{up}.
# \]
# Pressure-convection-diffusion and least-squares-commutator-type approximations exploit relationships between velocity and pressure operators. A BFBt/LSC implementation must be qualified from the actual discrete \(D,G,A_u\) blocks; naming the preconditioner is not sufficient.
#
# ## 9.4 AMG two-level principle
#
# Let \(P\) prolong coarse vectors and \(R\) restrict fine residuals:
# \[
# A_H=R A_hP.
# \]
# A two-level correction is
# \[
# x\leftarrow x+P A_H^{-1}R(b-A_hx).
# \]
# A smoother reduces error components poorly represented on the coarse grid; the coarse correction removes low-frequency/algebraically smooth components.
#
# ## 9.5 V-cycle
#
# \[
# x\rightarrow S_{pre}
# \rightarrow r=b-Ax
# \rightarrow r_H=Rr
# \rightarrow A_He_H=r_H
# \rightarrow x+Pe_H
# \rightarrow S_{post}.
# \]
# Coarse solve, smoother, transfer, damping and stopping rules together define the V-cycle.
#
# ## 9.6 Energy norm and contraction
#
# For SPD \(A\),
# \[
# \|e\|_A=\sqrt{e^TAe}.
# \]
# Define
# \[
# q_E=\frac{\|e_{out}\|_A}{\|e_{in}\|_A}.
# \]
# A hierarchy qualification should use an energy-norm contraction or another operator-level metric, not only a single right-hand-side residual ratio.
#
# ## 9.7 Ruge–Stüben and strength of connection
#
# AMG coarsening defines a strength relation, for example from a threshold applied to matrix couplings. Classical interpolation then builds
# \[
# x_F\approx P_{FC}x_C.
# \]
# A bug in strength/interpolation can produce a hierarchy that looks structurally valid but has poor or non-contracting error reduction; hence transfer and contraction tests belong in the qualification.
#
# ## 9.8 Smoothed aggregation
#
# Aggregation first forms groups of fine unknowns, constructs tentative prolongation \(P_0\), then smooths it:
# \[
# P=(I-\omega D^{-1}A)P_0
# \]
# in a representative Jacobi smoothing form. The exact smoother and near-nullspace vectors are part of the method.
#
# ## 9.9 MGR
#
# Split variables into retained \(C\) and eliminated \(F\):
# \[
# x=\begin{bmatrix}x_C\\x_F\end{bmatrix}.
# \]
# MGR recursively reduces the system while retaining selected physical variables. For incompressible flow this enables, for example, velocity-block reduction followed by a pressure/Schur hierarchy. Variable ordering and block maps are therefore first-class correctness data.
#
# ## 9.10 Matrix-value versus graph updates
#
# If only coefficients change,
# \[
# A^{new}=A^{old}+\Delta A
# \]
# with unchanged sparsity, numerical hierarchy data may need refresh. If connectivity changes, the graph and coarse sparsity may change and a structural rebuild can be required. Reusing a stale hierarchy after a topology change is not a valid optimisation.
#
# ## 9.11 Exact-oracle qualification
#
# A robust sequence is
# \[
# \boxed{
# \text{exact Schur}
# \rightarrow
# \text{approximate Schur}
# \rightarrow
# \text{Galerkin identity}
# \rightarrow
# \text{transfer correctness}
# \rightarrow
# \text{energy contraction}
# \rightarrow
# \text{matrix updates}
# \rightarrow
# \text{coupled CFD}
# }.
# \]
# Each stage isolates a different failure mode.
#
# ## 9.12 CFDX implementation
#
# Relevant sources include [amg_preconditioner.h](../../../src/cfdx/core/linalg/amg_preconditioner.h), [coupled_amg_schur.h](../../../src/cfdx/core/linalg/coupled_amg_schur.h), [block_schur.h](../../../src/cfdx/core/linalg/block_schur.h), [mgr_preconditioner.h](../../../src/cfdx/core/linalg/mgr_preconditioner.h), [linear_solver_dispatch.h](../../../src/cfdx/core/linalg/linear_solver_dispatch.h).
#
# Qualification tests include [test_exact_schur.cpp](../../../tests/unit/test_exact_schur.cpp), [test_schur_preconditioner.cpp](../../../tests/unit/test_schur_preconditioner.cpp), [test_mgr_preconditioner.cpp](../../../tests/unit/test_mgr_preconditioner.cpp), [test_amg_preconditioner_qualification.cpp](../../../tests/validation/test_amg_preconditioner_qualification.cpp) and [test_coupled_block_schur_amg.cpp](../../../tests/validation/test_coupled_block_schur_amg.cpp).
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
