# %% [markdown]
# 17 — How CFDX Solves the Fluid-Mechanics Equations
#
# This chapter is the **numerical-resolution course**: it starts from Navier–Stokes and follows every transformation until a converged CFD result exists.
#
# ## 17.1 The complete chain
#
# \[
# physical\ problem
# \rightarrow
# mathematical\ model
# \rightarrow
# PDE
# \rightarrow
# integral\ balance
# \rightarrow
# mesh/control\ volumes
# \rightarrow
# discrete\ operators
# \rightarrow
# algebraic\ system
# \rightarrow
# nonlinear\ iterations
# \rightarrow
# linear\ solves
# \rightarrow
# convergence\ diagnostics
# \rightarrow
# V\&V.
# \]
#
# Each arrow is a mathematical or software contract.
#
# ## 17.2 Step 1 — Define the physical problem
#
# Specify geometry, fluid/material properties, operating conditions, reference scales, initial state, boundary conditions and required outputs.
#
# The first question is not "which solver?". It is:
# \[
# \boxed{\text{What mathematical problem represents the physical question?}}
# \]
#
# ## 17.3 Step 2 — Choose the governing equations
#
# For incompressible Newtonian flow:
# \[
# \nabla\cdot\mathbf u=0,
# \]
# \[
# \rho\left(\partial_t\mathbf u+\mathbf u\cdot\nabla\mathbf u\right)
# =-\nabla p+\mu\nabla^2\mathbf u+\rho\mathbf f.
# \]
# Additional equations are introduced only when additional physics requires them: energy, species, turbulence, radiation, solid conduction, etc.
#
# ## 17.4 Step 3 — Build the mesh
#
# The domain is partitioned into control volumes (V_P). Each face has
# \[
# \mathbf S_f=A_f\mathbf n_f.
# \]
# Closed-cell geometry requires
# \[
# \sum_f\mathbf S_f=0.
# \]
# Mesh quality affects both truncation error and algebraic conditioning.
#
# ## 17.5 Step 4 — Integrate the PDE
#
# A generic transport equation
# \[
# \partial_t(\rho\phi)+\nabla\cdot(\rho\mathbf u\phi)
# =\nabla\cdot(\Gamma\nabla\phi)+S
# \]
# becomes
# \[
# \frac d{dt}\int_{V_P}\rho\phi\,dV
# +\sum_f\int_{A_f}\rho\mathbf u\phi\cdot\mathbf n_f\,dA
# =\sum_f\int_{A_f}\Gamma\nabla\phi\cdot\mathbf n_f\,dA
# +\int_{V_P}S\,dV.
# \]
# This is the point at which finite volume makes conservation explicit.
#
# ## 17.6 Step 5 — Reconstruct values and gradients
#
# Cell-centred unknowns need face values. Typical reconstruction is
# \[
# \phi_f=(1-w)\phi_P+w\phi_N.
# \]
# Gradients may use Green–Gauss:
# \[
# (\nabla\phi)_P\approx\frac1{V_P}\sum_f\phi_f\mathbf S_f,
# \]
# or least squares:
# \[
# A\mathbf g=\mathbf b.
# \]
# Higher-order schemes require additional reconstruction and often limiting.
#
# ## 17.7 Step 6 — Assemble fluxes
#
# Mass flux:
# \[
# \dot m_f=\rho_f\mathbf u_f\cdot\mathbf S_f.
# \]
# Convective scalar flux:
# \[
# F_{c,f}=\dot m_f\phi_f.
# \]
# Diffusive flux:
# \[
# F_{d,f}=-\Gamma_f(\nabla\phi)_f\cdot\mathbf S_f.
# \]
# Pressure force:
# \[
# \mathbf F_{p,f}=-p_f\mathbf S_f.
# \]
#
# ## 17.8 Step 7 — Obtain the algebraic system
#
# After discretisation:
# \[
# A(U)U=b(U).
# \]
# For a linearised scalar equation:
# \[
# a_P\phi_P=\sum_Na_{PN}\phi_N+b_P.
# \]
# Convection, diffusion, source terms, boundary conditions and temporal terms all contribute to the coefficients.
#
# ## 17.9 Step 8 — Treat nonlinear convection
#
# Navier–Stokes is nonlinear because of
# \[
# \mathbf u\cdot\nabla\mathbf u.
# \]
# A nonlinear iteration produces
# \[
# U^{k+1}=\mathcal S(U^k).
# \]
# Picard/segregated methods linearise coefficients using the previous iterate; Newton methods retain Jacobian coupling:
# \[
# J(U^k)\delta U=-R(U^k).
# \]
#
# ## 17.10 Step 9 — Pressure–velocity coupling
#
# In incompressible flow the system is a saddle-point problem:
# \[
# \begin{bmatrix}A_u&G\\D&C\end{bmatrix}
# \begin{bmatrix}u\\p\end{bmatrix}
# =
# \begin{bmatrix}b_u\\b_p\end{bmatrix}.
# \]
# Segregated SIMPLE/SIMPLEC/PISO/PIMPLE methods approximate this coupling. Coupled methods solve the block system more directly. The Schur complement is
# \[
# S=C-DA_u^{-1}G.
# \]
#
# ## 17.11 Step 10 — Linear solution
#
# Each linearised equation gives
# \[
# Ax=b.
# \]
# The true algebraic residual is
# \[
# r=b-Ax.
# \]
# Krylov methods such as CG and GMRES are combined with preconditioners. AMG/MGR/Schur methods address the multiscale and block structure of CFD systems.
#
# ## 17.12 Step 11 — Relaxation and stability
#
# A segregated update may use
# \[
# U^{k+1}\leftarrow(1-\omega)U^k+\omega U^{k+1}_{raw}.
# \]
# Under-relaxation can stabilise an iteration but does not make an incorrect discretisation correct. Convergence must be checked on the actual equations and physical balances.
#
# ## 17.13 Step 12 — Time advancement
#
# For
# \[
# M\frac{dU}{dt}=R(U),
# \]
# implicit Euler gives
# \[
# M\frac{U^{n+1}-U^n}{\Delta t}=R(U^{n+1}),
# \]
# while BDF2 gives
# \[
# M\frac{3U^{n+1}-4U^n+U^{n-1}}{2\Delta t}=R(U^{n+1}).
# \]
# Physical time and pseudo-time must never be conflated.
#
# ## 17.14 Step 13 — Convergence
#
# A robust convergence assessment uses several independent indicators:
#
# \[
# \|r\| \rightarrow 0,
# \qquad
# |Q^{k+1}-Q^k|\rightarrow0,
# qquad
# R_{conservation}\rightarrow0.
# \]
# Residual convergence alone is insufficient. A solver can reduce an algebraic residual while a quantity of interest, mass balance or boundary flux remains wrong.
#
# ## 17.15 Step 14 — Verification
#
# Verification asks whether the implemented numerical method solves the specified mathematical problem. Use:
#
# - exact polynomial fields;
# - manufactured solutions;
# - analytical laminar flows;
# - conservation identities;
# - mesh/time refinement;
# - linear-solver reference problems.
#
# For a refinement study:
# \[
# E_h\approx Ch^p,
# \qquad
# p_{obs}=\frac{\ln(E_h/E_{h/r})}{\ln r}.
# \]
#
# ## 17.16 Step 15 — Validation
#
# Validation compares the mathematical model against independent physical evidence. A comparison with another CFD implementation can reveal implementation differences, but it is not automatically experimental validation.
#
# ## 17.17 Step 16 — Post-processing
#
# Quantities of interest should be derived from the conservative fields:
# \[
# Q=\int_\Gamma\mathbf F\cdot\mathbf n\,dA,
# \]
# or appropriate volume/point functionals. The definition of Q must be fixed before comparing cases.
#
# ## 17.18 Step 17 — Restart and provenance
#
# A restart must reproduce the mathematical state required by the selected algorithms. Provenance should identify mesh, physics, numerics, solver settings and software revision.
#
# ## 17.19 CFDX software chain
#
# Case/data handling, mesh/geometry, numerical operators, physics assembly, linear algebra, output and diagnostics are separated in the source tree. The exact paths should be linked from the chapter sections after checking the current repository tree.
#
# ## 17.20 Relationship to the other chapters
#
# \[
# \boxed{
# 00\rightarrow01\rightarrow02\rightarrow03\rightarrow04\rightarrow05
# \rightarrow06\rightarrow07\rightarrow08\rightarrow09
# }
# \]
# is the main derivation path from fluid mechanics to a working incompressible CFD solver.
#
# Chapters 10–13 add closure and multiphysics; chapters 14–15 establish numerical evidence; chapters 16–18 connect the scientific model to data, execution and source traceability.
#
# %%
from __future__ import annotations
import numpy as np
A=np.array([[4.0,1.0],[1.0,3.0]])
b=np.array([1.0,2.0])
x=np.linalg.solve(A,b)
assert np.allclose(b-A@x,0.0)
