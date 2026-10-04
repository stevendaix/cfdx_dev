# %% [markdown]
# 07 — Pressure–Velocity Coupling
#
# ## 7.1 Governing constraint
#
# For constant-density incompressible flow:
# \[
# \nabla\cdot\mathbf u=0,
# \]
# \[
# \rho\left(
# \frac{\partial\mathbf u}{\partial t}
# +\nabla\cdot(\mathbf u\otimes\mathbf u)\right)
# =-\nabla p+\nabla\cdot(2\mu\mathbf S)+\mathbf f.
# \]
# The pressure is not obtained from an independent equation of state; it is the Lagrange multiplier enforcing discrete continuity.
#
# ## 7.2 Discrete block system and Schur complement
#
# \[
# \begin{bmatrix}
# A_u&G\\
# D&C
# \end{bmatrix}
# \begin{bmatrix}\mathbf u\\p\end{bmatrix}
# =
# \begin{bmatrix}b_u\\b_p\end{bmatrix}.
# \]
# Eliminating velocity gives
# \[
# S=C-DA_u^{-1}G,
# \qquad
# Sp=b_p-DA_u^{-1}b_u.
# \]
# This connects segregated pressure correction to the exact coupled algebraic problem.
#
# ## 7.3 SIMPLE derivation
#
# Write
# \[
# A_u\mathbf u=b_u-Gp.
# \]
# Let \(H\) approximate \(A_u^{-1}\). Then
# \[
# \mathbf u\approx H(b_u-Gp).
# \]
# Introduce corrections
# \[
# p^{new}=p+p',\qquad
# \mathbf u^{new}=\mathbf u+\mathbf u',
# \]
# with
# \[
# \mathbf u'\approx-HGp'.
# \]
# Continuity gives
# \[
# D\mathbf u'=-D\mathbf u,
# \]
# hence
# \[
# D H G\,p'=D\mathbf u.
# \]
# Sign changes in this equation occur if the code defines the pressure-gradient block with the opposite sign; the implementation must be audited against the assembled row rather than against a memorised textbook formula.
#
# ## 7.4 SIMPLEC
#
# SIMPLEC changes the approximation for the velocity correction so that additional off-diagonal momentum coupling is retained. Its identity is therefore encoded in the exact approximation \(H\), not merely in an enum value.
#
# ## 7.5 PISO
#
# Starting from a predictor
# \[
# A_u\mathbf u^*=b_u-Gp^n,
# \]
# PISO applies successive pressure corrections within the same time level. Each correction updates the face flux/velocity and generates a new continuity defect. The number of correction stages is an algorithmic parameter.
#
# ## 7.6 PIMPLE
#
# PIMPLE combines outer nonlinear iterations with pressure-correction stages:
# \[
# \boxed{\text{outer iteration}
# \rightarrow\text{momentum}
# \rightarrow\text{pressure corrections}
# \rightarrow\text{relaxation}
# \rightarrow\text{convergence}}.
# \]
# Outer iteration and inner pressure-correction convergence must be reported separately.
#
# ## 7.7 Fractional-step projection
#
# Predictor:
# \[
# \frac{\mathbf u^*-\mathbf u^n}{\Delta t}=R(\mathbf u^n).
# \]
# Pressure equation:
# \[
# \nabla^2p^{n+1}
# =\frac{\rho}{\Delta t}\nabla\cdot\mathbf u^*.
# \]
# Projection:
# \[
# \mathbf u^{n+1}
# =\mathbf u^*-\frac{\Delta t}{\rho}\nabla p^{n+1}.
# \]
# Taking the divergence yields
# \[
# \nabla\cdot\mathbf u^{n+1}=0
# \]
# if the pressure equation and discrete gradient/divergence pair are exactly consistent.
#
# ## 7.8 Pressure gauge and compatibility
#
# Under pure Neumann pressure conditions,
# \[
# p\rightarrow p+C
# \]
# is a null-space transformation. The pressure Poisson equation also requires a compatibility condition:
# \[
# \int_\Omega b_p\,dV
# =\int_{\partial\Omega}\frac{\partial p}{\partial n}\,dA
# \]
# in the continuous case. A discrete incompatibility can prevent convergence even when the operator itself is correct.
#
# ## 7.9 Collocated face mass flux
#
# A collocated scheme needs a pressure-velocity interpolation that avoids an odd-even pressure mode. A Rhie–Chow-type flux has the generic structure
# \[
# \dot m_f=
# \dot m_f^{interp}
# -D_f\left[
# (p_N-p_P)-(\nabla p)_f\cdot\mathbf d_{PN}
# \right],
# \]
# where \(D_f\) is derived from the momentum diagonal/coefficient. The exact CFDX formula, including signs and boundary handling, must be taken from the implementation.
#
# ## 7.10 Coupling qualification
#
# A complete test is not just \(\|r\|\rightarrow0\). It should include
# \[
# \|D\mathbf u\|,\quad
# \|r_u\|,\quad
# \|r_p\|,\quad
# \Delta p\text{ gauge invariance},
# \]
# plus conservation and canonical flow quantities. The coupling method must be tested on more than one initial condition because a single converged state can hide a checkerboard or null-space defect.
#
# ## 7.11 CFDX implementation
#
# See [pressure_velocity.h](../../../src/cfdx/physics/pressure_velocity.h), [pressure_velocity_algorithms.h](../../../src/cfdx/physics/pressure_velocity_algorithms.h), [steady_incompressible_solver.h](../../../src/cfdx/physics/steady_incompressible_solver.h) and [numerical_method_registry.h](../../../src/cfdx/core/numerics/numerical_method_registry.h).
#
# ## 7.12 Executable Schur check
# %%
import numpy as np
Au=np.diag([2.,3.]); G=np.array([[1.],[2.]])
D=G.T; C=np.zeros((1,1))
S=C-D@np.linalg.solve(Au,G)
assert S.shape==(1,1)
assert np.isfinite(S).all()
