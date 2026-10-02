# %% [markdown]
# 07 — Pressure–Velocity Coupling
#
# ## 7.1 Incompressible constraint
#
# For constant-density incompressible flow:
# \[
# \nabla\cdot\mathbf u=0,
# \]
# while momentum is
# \[
# \rho\left(
# \frac{\partial\mathbf u}{\partial t}
# +\nabla\cdot(\mathbf u\otimes\mathbf u)
# \right)
# =-\nabla p+\nabla\cdot(2\mu\mathbf S)+\mathbf f,
# \quad
# \mathbf S=\frac12(\nabla\mathbf u+\nabla\mathbf u^T).
# \]
# Pressure is a constraint variable: it adjusts so that the discrete velocity satisfies the discrete continuity equation.
#
# ## 7.2 Discrete saddle-point system
#
# Linearisation gives
# \[
# \begin{bmatrix}
# A_u&G\\
# D&C
# \end{bmatrix}
# \begin{bmatrix}\mathbf u\\p\end{bmatrix}
# =
# \begin{bmatrix}b_u\\b_p\end{bmatrix}.
# \]
# \(G\) is a pressure-gradient block and \(D\) a velocity-divergence block. In an ideal conservative pairing their discrete relationship is constrained by the chosen inner products and boundary treatment.
#
# Eliminating velocity:
# \[
# S=C-DA_u^{-1}G,
# \]
# \[
# Sp=b_p-DA_u^{-1}b_u.
# \]
# This is the mathematical basis for pressure-Schur preconditioning.
#
# ## 7.3 SIMPLE
#
# Split the momentum equation:
# \[
# A_u u=b_u-Gp.
# \]
# Using an approximate inverse \(H\approx A_u^{-1}\):
# \[
# u\approx H(b_u-Gp).
# \]
# A pressure correction \(p'=p^{new}-p\) is obtained from continuity:
# \[
# D H G\,p'=D H(b_u-Gp)-b_p.
# \]
# The actual CFDX implementation may use diagonal or another approximation; the approximation defines the algorithm.
#
# ## 7.4 SIMPLEC
#
# SIMPLEC changes the approximation used for the velocity correction so that less of the momentum-correction coupling is discarded. It is not simply a renamed SIMPLE implementation. The exact coefficient correction must be documented and verified against the implemented matrix row.
#
# ## 7.5 PISO
#
# PISO performs multiple pressure-correction stages within a time step:
# \[
# u^{*}=A_u^{-1}(b_u-Gp^n),
# \]
# followed by pressure correction and velocity correction, then repeated correction using the updated flux/velocity. The number and ordering of corrections are part of the algorithm contract.
#
# ## 7.6 PIMPLE
#
# PIMPLE introduces outer nonlinear iterations around PISO-like inner corrections:
# \[
# \text{outer loop}\rightarrow
# \text{momentum}\rightarrow
# \text{pressure corrections}\rightarrow
# \text{convergence test}.
# \]
# Relaxation and stopping criteria must be separated from the underlying pressure-correction equations.
#
# ## 7.7 Fractional step
#
# A projection method first computes
# \[
# \frac{u^*-u^n}{\Delta t}=R(u^n),
# \]
# then solves
# \[
# \nabla^2p^{n+1}
# =\frac{\rho}{\Delta t}\nabla\cdot u^*,
# \]
# and projects:
# \[
# u^{n+1}=u^*-\frac{\Delta t}{\rho}\nabla p^{n+1}.
# \]
# Boundary conditions for \(p\) and \(u\) must be derived consistently; the projection is not automatically equivalent to a collocated FVM SIMPLE formulation.
#
# ## 7.8 Pressure null space
#
# With pure Neumann pressure boundaries,
# \[
# p'=p+C
# \]
# leaves \(\nabla p\) unchanged. The discrete pressure matrix is singular unless a gauge or explicit null-space constraint is imposed. A solver that converges only after arbitrary pressure shifting has not necessarily demonstrated correct gauge handling.
#
# ## 7.9 Rhie–Chow / face flux
#
# Collocated arrangements can admit pressure checkerboarding if pressure and velocity are interpolated independently. A Rhie–Chow-type momentum interpolation introduces a pressure-difference correction into the face mass flux. The exact formula, coefficient scaling and boundary behaviour must be taken from the CFDX implementation before claiming equivalence to a named method.
#
# ## 7.10 Verification
#
# Required operator checks include:
# \[
# D u=0,
# \qquad
# \|r_u\|\rightarrow0,
# \qquad
# \|r_p\|\rightarrow0,
# \]
# plus pressure-gauge invariance, manufactured solutions and canonical incompressible flows. Couette and Poiseuille provide particularly useful low-complexity checks; Ghia tests the coupled steady solver at a more demanding level.
#
# ## 7.11 CFDX implementation
#
# See [pressure_velocity.h](../../../src/cfdx/physics/pressure_velocity.h), [pressure_velocity_algorithms.h](../../../src/cfdx/physics/pressure_velocity_algorithms.h) and [steady_incompressible_solver.h](../../../src/cfdx/physics/steady_incompressible_solver.h).
#
# Verification includes [test_steady_incompressible_solver.cpp](../../../tests/validation/test_steady_incompressible_solver.cpp) and incompressible unit tests.
#
# ## 7.12 Executable Schur check
# %%
import numpy as np
Au=np.diag([2.,3.]); G=np.array([[1.],[2.]])
D=G.T; C=np.zeros((1,1))
S=C-D@np.linalg.solve(Au,G)
assert S.shape==(1,1)
assert np.isfinite(S).all()
