# CFDX Theory — Complete Scientific Reference

> **Purpose.** This document is the mathematical spine of CFDX. It is written as a course and as a traceability document: every important numerical statement is connected to a continuous equation, its finite-volume discretisation, its numerical properties, the corresponding CFDX source path, and the relevant verification/validation evidence.
>
> **Evidence rule.** *Implemented* means code exists. *Verified* means a defined mathematical property has executable evidence. *Validated* means an independent physical/reference comparison has been completed. *Qualified* means the declared population and acceptance gates support the intended engineering scope. Documentation alone never changes these states.

## 0. Foundations

### 0.1 Continuum mechanics

CFDX represents macroscopic fields
\[
\rho=\rho(\mathbf x,t),\qquad
\mathbf u=\mathbf u(\mathbf x,t),\qquad
p=p(\mathbf x,t),\qquad
T=T(\mathbf x,t).
\]
The material derivative is
\[
\frac{D\phi}{Dt}=\frac{\partial\phi}{\partial t}+\mathbf u\cdot\nabla\phi.
\]
For a material volume, conservation means that the rate of change of an extensive quantity equals fluxes plus sources.

**Repository:** [field](../../src/cfdx/core/field/field.h), [mesh](../../src/cfdx/core/mesh/mesh.h), [physics](../../src/cfdx/physics/).

### 0.2 Kinematics

The velocity gradient is decomposed as
\[
\nabla\mathbf u=\mathbf D+\mathbf W,
\quad
\mathbf D=\frac12(\nabla\mathbf u+\nabla\mathbf u^T),
\quad
\mathbf W=\frac12(\nabla\mathbf u-\nabla\mathbf u^T).
\]
\(\mathbf D\) is the rate-of-deformation tensor; \(\mathbf W\) is the spin tensor. The divergence \(\nabla\cdot\mathbf u=\mathrm{tr}(\mathbf D)\) measures volumetric expansion.

### 0.3 Conservation laws

For a conserved density q and flux \(\mathbf F_q\),
\[
\frac{\partial q}{\partial t}+\nabla\cdot\mathbf F_q=s_q.
\]
Integrating over a fixed control volume V and applying Gauss' theorem gives
\[
\frac{d}{dt}\int_Vq\,dV+
\oint_{\partial V}\mathbf F_q\cdot\mathbf n\,dA
=\int_Vs_q\,dV.
\]
This integral equation is the starting point of the finite-volume method.

### 0.4 Constitutive laws

For a Newtonian fluid,
\[
\boldsymbol\sigma=-p\mathbf I+\boldsymbol\tau,
\]
\[
\boldsymbol\tau=2\mu\mathbf D+
\lambda(\nabla\cdot\mathbf u)\mathbf I.
\]
For incompressible constant-viscosity flow,
\[
\boldsymbol\tau=2\mu\mathbf D,
\qquad
\nabla\cdot\mathbf u=0.
\]
Fourier conduction is
\[
\mathbf q=-k\nabla T.
\]

### 0.5 Incompressible Navier–Stokes

\[
\nabla\cdot\mathbf u=0,
\]
\[
\rho\left(
\frac{\partial\mathbf u}{\partial t}
+\mathbf u\cdot\nabla\mathbf u
\right)
=-\nabla p+\mu\nabla^2\mathbf u+\rho\mathbf f.
\]
The pressure is a Lagrange multiplier enforcing continuity; it is not an independently transported scalar.

**Repository:** [incompressible.h](../../src/cfdx/physics/incompressible.h), [steady_incompressible_solver.h](../../src/cfdx/physics/steady_incompressible_solver.h).

### 0.6 Boundary and initial conditions

Dirichlet:
\[
\phi=\phi_b.
\]
Neumann:
\[
\nabla\phi\cdot\mathbf n=g_N.
\]
Robin:
\[
a\phi+b\nabla\phi\cdot\mathbf n=c.
\]
For a transient problem, \(\phi(\mathbf x,0)=\phi_0\) is also required. CFDX uses explicit boundary abstractions rather than silently inferring a physical condition.

**Repository:** [boundary_condition.h](../../src/cfdx/core/boundary/boundary_condition.h), [boundary_field.h](../../src/cfdx/core/boundary/boundary_field.h), [mathematical_condition.h](../../src/cfdx/core/boundary/mathematical_condition.h).

### 0.7 Dimensional analysis

\[
Re=\frac{\rho UL}{\mu},\qquad
Ma=\frac{U}{a},\qquad
Pe=\frac{UL}{\alpha},\qquad
Pr=\frac{\nu}{\alpha}.
\]
Nondimensionalisation exposes dominant balances and provides a consistency check: every term in a governing equation must have identical dimensions.

---

# 1. Conservation laws

For a generic finite-volume cell P,
\[
\frac{d}{dt}\int_{V_P}q\,dV+
\sum_{f\in\partial P}\int_{A_f}\mathbf F_q\cdot\mathbf n_f\,dA
=\int_{V_P}s_q\,dV.
\]

## 1.1 Mass

\[
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0.
\]
For a cell,
\[
\frac{d}{dt}(\rho_PV_P)+
\sum_f\dot m_f=0,
\qquad
\dot m_f=\rho_f\mathbf u_f\cdot\mathbf S_f,
\quad
\mathbf S_f=\mathbf n_fA_f.
\]
An internal face has opposite oriented area vectors for its two cells, so a single physical face flux must cancel exactly between neighbours.

## 1.2 Momentum

\[
\frac{\partial(\rho\mathbf u)}{\partial t}
+\nabla\cdot(\rho\mathbf u\otimes\mathbf u)
=-\nabla p+\nabla\cdot\boldsymbol\tau+\rho\mathbf f.
\]
The integral form is a sum of convective momentum flux, pressure traction, viscous traction and body force.

## 1.3 Energy

With \(E=e+|\mathbf u|^2/2\),
\[
\frac{\partial(\rho E)}{\partial t}
+\nabla\cdot[(\rho E+p)\mathbf u]
=\nabla\cdot(\boldsymbol\tau\mathbf u-\mathbf q)
+\rho\mathbf f\cdot\mathbf u+\dot q_v.
\]
For pure steady conduction with constant k,
\[
\nabla\cdot(k\nabla T)+\dot q_v=0.
\]

## 1.4 Species

For species i,
\[
\frac{\partial(\rho Y_i)}{\partial t}
+\nabla\cdot(\rho\mathbf uY_i)
=-\nabla\cdot\mathbf J_i+\dot\omega_i.
\]
A consistent mixture model must also satisfy
\[
\sum_iY_i=1,
\qquad
\sum_i\dot\omega_i=0
\]
when reactions conserve total mass.

## 1.5 Global conservation

Summing all cell balances gives
\[
\frac{d}{dt}\sum_Pq_PV_P+
\sum_{f\in\partial\Omega}\Phi_f
=\sum_Ps_PV_P.
\]
All internal faces cancel pairwise. Therefore a global conservation diagnostic can detect assembly errors independently of local solver residuals.

**Repository:** [conservation.h](../../src/cfdx/core/numerics/conservation.h), [flux.h](../../src/cfdx/core/numerics/flux.h), [integrate.h](../../src/cfdx/core/numerics/integrate.h), [finite_volume_transport.h](../../src/cfdx/physics/finite_volume_transport.h).

---

# 2. Finite-volume method

The FVM approximates the integral conservation equation directly. For a scalar transport equation
\[
\frac{\partial(\rho\phi)}{\partial t}
+\nabla\cdot(\rho\mathbf u\phi)
=\nabla\cdot(\Gamma\nabla\phi)+S,
\]
integration gives
\[
\frac{\rho_PV_P\phi_P^{n+1}-\rho_PV_P\phi_P^n}{\Delta t}
+\sum_f F_f\phi_f
=\sum_f\Gamma_f(\nabla\phi)_f\cdot\mathbf S_f+S_PV_P.
\]

For steady linearised transport this becomes
\[
a_P\phi_P=\sum_Na_{PN}\phi_N+b_P.
\]
The coefficients depend on convection, diffusion, source linearisation and boundary conditions.

## 2.1 Geometry and orientation

\[
\mathbf S_f=\mathbf n_fA_f,
\qquad
\sum_f\mathbf S_f=\mathbf0
\]
for a closed cell. The second identity is a geometric closure condition and should be exact up to the representation precision.

## 2.2 Diffusion

The face-normal diffusive flux is
\[
\Phi_d=-\Gamma_f(\nabla\phi)_f\cdot\mathbf S_f.
\]
For an orthogonal two-point approximation,
\[
(\nabla\phi)_f\cdot\mathbf n_f
\approx
\frac{\phi_N-\phi_P}{d_{PN}}.
\]
On non-orthogonal meshes, the centre-to-centre vector is not parallel to \(\mathbf n_f\); a correction must therefore be represented explicitly rather than hidden in a scalar distance.

## 2.3 Matrix properties

A robust discretisation must make the intended properties explicit: conservation, consistency, boundedness/monotonicity where applicable, diagonal dominance or another suitable stability property, and correct symmetry when the continuous operator is self-adjoint.

**Repository:** [laplacian.h](../../src/cfdx/core/numerics/laplacian.h), [divergence.h](../../src/cfdx/core/numerics/divergence.h), [convection.h](../../src/cfdx/core/numerics/convection.h), [source_term.h](../../src/cfdx/core/numerics/source_term.h).

---

# 3. Meshes

A polyhedral mesh consists of points, faces and cells. A cell is represented by an oriented collection of faces. For a valid closed polyhedron,
\[
\sum_f\mathbf S_f=0.
\]
Its volume can be evaluated by decomposing the polyhedron into signed pyramids/tetrahedra. For a planar face f with area vector \(\mathbf S_f\) and a reference point \(\mathbf x_0\), the signed volume contribution has the form
\[
V_f=\frac13A_f\,[(\mathbf x_f-\mathbf x_0)\cdot\mathbf n_f].
\]
The implementation must use a consistent geometric convention.

Quality is multidimensional. Typical indicators include non-orthogonality
\[
\theta_f=\cos^{-1}
\left(
\frac{\mathbf d_{PN}\cdot\mathbf n_f}
{|\mathbf d_{PN}|}
\right),
\]
skewness measured from the face-centre/centre-line mismatch, aspect ratio, and cell-volume positivity.

Non-planar polygonal faces require an area-vector and centroid convention that remains consistent with the cell closure equations.

**Repository:** [mesh.h](../../src/cfdx/core/mesh/mesh.h), [face.h](../../src/cfdx/core/mesh/face.h), [cell.h](../../src/cfdx/core/mesh/cell.h), [face_geometry.h](../../src/cfdx/core/geometry/face_geometry.h), [cell_geometry.h](../../src/cfdx/core/geometry/cell_geometry.h), [mesh_quality.h](../../src/cfdx/core/geometry/mesh_quality.h), [mesh_validator.h](../../src/cfdx/core/geometry/mesh_validator.h).

---

# 4. Gradients and reconstruction

This is a central N2 chapter. The gradient is needed by diffusion, viscous stresses, turbulence production, reconstruction and many limiters.

## 4.1 Taylor foundation

For a smooth scalar field,
\[
\phi(\mathbf x_P+\mathbf d)
=
\phi_P+\nabla\phi_P\cdot\mathbf d
+\frac12\mathbf d^T\mathbf H_P\mathbf d+O(h^3).
\]
Thus
\[
\Delta\phi=\nabla\phi_P\cdot\mathbf d+O(h^2).
\]
A gradient method is second-order in smooth regions only if the geometry/stencil and all reconstruction operations preserve the required cancellation of the Taylor remainder.

## 4.2 Green–Gauss

Gauss' theorem gives
\[
\int_{V_P}\nabla\phi\,dV
=
\oint_{\partial V_P}\phi\mathbf n\,dA,
\]
hence
\[
(\nabla\phi)_P
\approx
\frac1{V_P}\sum_f\phi_f\mathbf S_f.
\]
With linear cell interpolation,
\[
\phi_f=(1-w_f)\phi_P+w_f\phi_N,
\qquad
w_f=
\frac{|\mathbf x_f-\mathbf x_P|}
{|\mathbf x_N-\mathbf x_P|}.
\]
The method is conservative in the face-integral sense, but its accuracy depends on face interpolation and mesh geometry. OpenFOAM documents Gauss gradients as a standard second-order option and exposes the interpolation choice explicitly. citeturn0search0turn0search2

## 4.3 Least squares

For neighbours N define
\[
\mathbf d_N=\mathbf x_N-\mathbf x_P,
\qquad
\Delta\phi_N=\phi_N-\phi_P.
\]
Minimise
\[
J(\mathbf g)=
\sum_Nw_N
(\Delta\phi_N-\mathbf g\cdot\mathbf d_N)^2.
\]
The normal equations are
\[
\mathbf A\mathbf g=\mathbf b,
\]
\[
\mathbf A=\sum_Nw_N\mathbf d_N\mathbf d_N^T,
\qquad
\mathbf b=\sum_Nw_N\mathbf d_N\Delta\phi_N.
\]
In 3-D, a unique gradient requires
\[
\mathrm{rank}(\mathbf A)=3.
\]
A robust implementation must detect rank deficiency rather than silently inverting a singular matrix.

OpenFOAM exposes both least-squares gradients and arbitrary weighted least-squares stencils; SU2 exposes Green–Gauss and weighted least-squares reconstruction choices. citeturn0search4turn0search7turn0search12

## 4.4 Weighting and conditioning

Common choices include
\[
w_N=|\mathbf d_N|^{-p},
\qquad p>0.
\]
The conditioning indicator
\[
\kappa(\mathbf A)=
\frac{\lambda_{\max}(\mathbf A)}
{\lambda_{\min}(\mathbf A)}
\]
measures sensitivity to perturbations. For a symmetric positive-definite system,
\[
\frac{\|\delta\mathbf g\|}
{\|\mathbf g\|}
\lesssim
\kappa(\mathbf A)
\frac{\|\delta\mathbf b\|}
{\|\mathbf b\|}.
\]
Conditioning must therefore be diagnosed before claiming robustness on arbitrary polyhedra.

## 4.5 Boundary reconstruction

A boundary face does not provide a second cell value. The reconstruction policy must therefore distinguish:
1. internal-cell neighbours;
2. physical boundary constraints;
3. ghost/equivalent values, if used;
4. one-sided or constrained reconstruction.

The boundary treatment must be compatible with the physical BC and must not introduce an undocumented order reduction.

## 4.6 Face reconstruction

The face value and the cell gradient are separate abstractions:
\[
\phi_f=\mathcal R_f(\phi_P,\phi_N,\nabla\phi_P,\nabla\phi_N,\mathbf x_f,\ldots).
\]
A gradient implementation must not accidentally determine the face interpolation policy. This separation is necessary for future MUSCL/TVD schemes.

## 4.7 Limiting

A generic limited reconstruction is
\[
\phi_f=
\phi_P+
\alpha_P\nabla\phi_P\cdot
(\mathbf x_f-\mathbf x_P),
\qquad
0\le\alpha_P\le1.
\]
A limiter is designed to prevent extrapolated values from violating a prescribed local bound. It trades formal accuracy against nonlinear boundedness. OpenFOAM explicitly provides cell- and face-limited versions of base gradient schemes. citeturn0search0turn0search5

## 4.8 Verification

For a manufactured polynomial
\[
\phi(x,y,z)=c_0+c_xx+c_yy+c_zz,
\]
the exact gradient is constant:
\[
\nabla\phi=(c_x,c_y,c_z)^T.
\]
This is the minimum exactness test. For a smooth non-linear field, use
\[
E_h=
\left(
\frac{\sum_PV_P|\nabla\phi_h-\nabla\phi|^2}
{\sum_PV_P|\nabla\phi|^2}
\right)^{1/2},
\]
and
\[
p_{obs}=
\frac{\log(E_h/E_{h/2})}{\log 2}.
\]
The mesh sequence must be controlled; a single mesh cannot establish order.

**Repository:** [gradient.h](../../src/cfdx/core/numerics/gradient.h), [gradient_stencil.h](../../src/cfdx/core/numerics/gradient_stencil.h), [least_squares_gradient.h](../../src/cfdx/core/fvm/least_squares_gradient.h), [interpolation.h](../../src/cfdx/core/numerics/interpolation.h); verification: [test_gradient_verification.cpp](../../tests/validation/test_gradient_verification.cpp), [test_polyhedral_gradient_campaign.cpp](../../tests/validation/test_polyhedral_gradient_campaign.cpp).

**Current evidence status:** N2 remains **PARTIAL / not qualified** until the required polyhedral accuracy, conditioning, boundary reconstruction and face-reconstruction evidence is complete.

---

# 5. Fluxes

Fluxes are the conservative interface between adjacent control volumes.

Mass:
\[
F_m=\rho_f\mathbf u_f\cdot\mathbf S_f.
\]
Convective scalar flux:
\[
F_c=F_m\phi_f.
\]
Diffusive flux:
\[
F_d=-\Gamma_f(\nabla\phi)_f\cdot\mathbf S_f.
\]
Pressure force on a face:
\[
\mathbf F_p=-p_f\mathbf S_f.
\]
The total face contribution must use one oriented flux for both adjacent cells:
\[
F_{P,f}=-F_{N,f}.
\]
This is the discrete statement of internal conservation.

For convection, first-order upwind can be written
\[
\phi_f=
\begin{cases}
\phi_P,&F_m>0,\\
\phi_N,&F_m<0.
\end{cases}
\]
Central interpolation is
\[
\phi_f=(1-w)\phi_P+w\phi_N.
\]
Higher-order reconstruction adds gradient information but requires boundedness controls where the solution contains strong extrema.

**Repository:** [flux.h](../../src/cfdx/core/numerics/flux.h), [convection.h](../../src/cfdx/core/numerics/convection.h), [interpolation.h](../../src/cfdx/core/numerics/interpolation.h), [finite_volume_transport.h](../../src/cfdx/physics/finite_volume_transport.h); tests: [test_flux.cpp](../../tests/unit/test_flux.cpp), [test_conservation_boundedness.cpp](../../tests/unit/test_conservation_boundedness.cpp), [test_convection_scheme_verification.cpp](../../tests/validation/test_convection_scheme_verification.cpp).

---

# 6. Time integration

After spatial discretisation,
\[
\mathbf M\frac{d\mathbf U}{dt}=\mathbf R(\mathbf U,t).
\]
Explicit Euler:
\[
\mathbf U^{n+1}
=
\mathbf U^n+
\Delta t\,\mathbf M^{-1}\mathbf R(\mathbf U^n).
\]
Implicit Euler:
\[
\mathbf M\frac{\mathbf U^{n+1}-\mathbf U^n}{\Delta t}
=\mathbf R(\mathbf U^{n+1}).
\]
Crank–Nicolson:
\[
\mathbf M\frac{\mathbf U^{n+1}-\mathbf U^n}{\Delta t}
=
\frac12[\mathbf R(\mathbf U^{n+1})+\mathbf R(\mathbf U^n)].
\]
BDF2:
\[
\mathbf M
\frac{3\mathbf U^{n+1}-4\mathbf U^n+\mathbf U^{n-1}}
{2\Delta t}
=
\mathbf R(\mathbf U^{n+1}).
\]
For variable timesteps, the coefficients must be derived from the actual time nodes rather than reusing constant-step BDF2 coefficients.

A representative explicit convective restriction is
\[
CFL=\frac{|u|\Delta t}{\Delta x}\lesssim C_{max},
\]
with the exact limit determined by the spatial operator and dimensionality.

**Repository:** [temporal.h](../../src/cfdx/physics/temporal.h), [low_storage_time_integration.h](../../src/cfdx/physics/low_storage_time_integration.h), [adaptive_cfl.h](../../src/cfdx/physics/adaptive_cfl.h), [core temporal](../../src/cfdx/core/numerics/temporal.h); tests: [test_temporal_physics.cpp](../../tests/unit/test_temporal.cpp), [test_temporal_order_matrix.cpp](../../tests/validation/test_temporal_order_matrix.cpp).

---

# 7. Pressure–velocity coupling

The incompressible discrete system can be written
\[
\begin{bmatrix}
A_u&G\\
D&0
\end{bmatrix}
\begin{bmatrix}
\mathbf u\\p
\end{bmatrix}
=
\begin{bmatrix}
b_u\\b_p
\end{bmatrix}.
\]
Eliminating velocity gives a Schur complement:
\[
S p=
b_p-D A_u^{-1}b_u,
\qquad
S=-DA_u^{-1}G
\]
up to the sign convention used by the code.

SIMPLE constructs an approximate pressure correction from an approximate inverse of the momentum operator. SIMPLEC changes that approximation to improve coupling consistency. PISO performs multiple pressure corrections within a timestep/outer iteration. PIMPLE combines outer nonlinear iteration with PISO-like corrections. Fractional-step methods split momentum and projection. A coupled method solves the block system more directly.

The pressure Poisson equation is singular for pure Neumann pressure conditions:
\[
p\rightarrow p+C.
\]
A gauge or null-space treatment is therefore required.

Rhie–Chow-type pressure-velocity interpolation is relevant on collocated meshes because it suppresses checkerboard pressure modes, but its exact implementation must be audited rather than assumed from the function name.

**Repository:** [pressure_velocity.h](../../src/cfdx/physics/pressure_velocity.h), [pressure_velocity_algorithms.h](../../src/cfdx/physics/pressure_velocity_algorithms.h), [steady_incompressible_solver.h](../../src/cfdx/physics/steady_incompressible_solver.h); tests: [test_incompressible.cpp](../../tests/unit/test_incompressible.cpp), [test_steady_incompressible_solver.cpp](../../tests/validation/test_steady_incompressible_solver.cpp), [test_schur_infrastructure.cpp](../../tests/unit/test_schur_infrastructure.cpp).

---

# 8. Linear algebra

Every implicit discretisation produces
\[
A x=b.
\]
The true residual is
\[
r=b-Ax.
\]
A relative residual can be defined as
\[
\eta=
\frac{\|r\|}
{\|b\|+\|A\|\|x\|}
\]
or with another explicitly documented normalisation. Residual definitions must never be mixed silently.

For SPD systems, CG is appropriate:
\[
x_{k+1}=x_k+\alpha_kp_k,
\]
with conjugate directions. GMRES constructs a Krylov approximation
\[
x_k\in x_0+\mathcal K_k(A,r_0),
\qquad
\mathcal K_k=
\mathrm{span}\{r_0,Ar_0,\ldots,A^{k-1}r_0\}.
\]
Preconditioning replaces the difficult system by an equivalent better-conditioned one, e.g.
\[
M^{-1}Ax=M^{-1}b.
\]
Conditioning is a property of the mathematical operator; convergence rate is also controlled by non-normality, spectrum and preconditioner quality.

**Repository:** [sparse_matrix.h](../../src/cfdx/core/linalg/sparse_matrix.h), [linear_system.h](../../src/cfdx/core/linalg/linear_system.h), [cg_solver.h](../../src/cfdx/core/linalg/cg_solver.h), [gmres_solver.h](../../src/cfdx/core/linalg/gmres_solver.h), [linear_solver_dispatch.h](../../src/cfdx/core/linalg/linear_solver_dispatch.h); tests: [test_cg_solver.cpp](../../tests/unit/test_cg_solver.cpp), [test_gmres_solver.cpp](../../tests/unit/test_gmres_solver.cpp), [test_krylov_preconditioning.cpp](../../tests/unit/test_krylov_preconditioning.cpp).

---

# 9. AMG, MGR and Schur

For a block system
\[
A=
\begin{bmatrix}
A_{uu}&A_{up}\\
A_{pu}&A_{pp}
\end{bmatrix},
\]
the exact pressure Schur complement is
\[
S=A_{pp}-A_{pu}A_{uu}^{-1}A_{up}.
\]
Because applying \(A_{uu}^{-1}\) exactly is expensive, CFDX can use an approximation \(\tilde S\).

AMG represents a hierarchy
\[
A_0\rightarrow A_1\rightarrow\cdots\rightarrow A_L
\]
with restriction R, prolongation P and a coarse operator
\[
A_c=R A_f P.
\]
A two-grid correction is conceptually
\[
x\leftarrow x+P A_c^{-1}R(b-Ax).
\]
A V-cycle combines smoothing on fine/coarse levels and coarse correction.

The meaningful convergence metric for an elliptic preconditioner should not be limited to a residual norm. An energy norm is
\[
\|e\|_A=\sqrt{e^TAe},
\]
and an energy contraction factor can be defined as
\[
\rho_E=
\frac{\|e_{out}\|_A}{\|e_{in}\|_A}.
\]
MGR generalises multilevel reduction to block-variable systems; the exact ordering, coarse variables and relaxation strategy are part of the algorithm and must be documented.

**Repository:** [exact_schur.h](../../src/cfdx/core/linalg/exact_schur.h), [block_schur.h](../../src/cfdx/core/linalg/block_schur.h), [schur_approximation.h](../../src/cfdx/core/linalg/schur_approximation.h), [amg_preconditioner.h](../../src/cfdx/core/linalg/amg_preconditioner.h), [hypre_amg.h](../../src/cfdx/core/linalg/hypre_amg.h), [mgr_preconditioner.h](../../src/cfdx/core/linalg/mgr_preconditioner.h), [coupled_amg_schur.h](../../src/cfdx/core/linalg/coupled_amg_schur.h); tests: [test_exact_schur.cpp](../../tests/unit/test_exact_schur.cpp), [test_schur_preconditioner.cpp](../../tests/unit/test_schur_preconditioner.cpp), [test_mgr_preconditioner.cpp](../../tests/unit/test_mgr_preconditioner.cpp), [test_coupled_block_schur_amg.cpp](../../tests/validation/test_coupled_block_schur_amg.cpp).

---

# 10. Turbulence

Reynolds decomposition:
\[
\mathbf u=\overline{\mathbf u}+\mathbf u'.
\]
Averaging the nonlinear convection term introduces the Reynolds stress:
\[
-\rho\overline{u_i'u_j'}.
\]
The closure problem is to model this tensor.

The eddy-viscosity hypothesis gives
\[
-\rho\overline{u_i'u_j'}
=
2\mu_t\mathbf S
-\frac23\rho k\mathbf I,
\]
with
\[
k=\frac12\overline{u_i'u_i'}.
\]
RANS models add transport equations for turbulence quantities. A generic k-equation has the structure
\[
\frac{Dk}{Dt}
=P_k-\beta^*k\omega
+\nabla\cdot[(\nu+\sigma_k\nu_t)\nabla k].
\]
The SST family blends near-wall k-omega behaviour with outer-flow k-epsilon behaviour. Spalart–Allmaras solves a transported working variable with wall-distance dependence.

LES filters the equations:
\[
\bar u_i=G*u_i,
\]
creating subgrid stresses
\[
\tau_{ij}^{sgs}=\overline{u_iu_j}-\bar u_i\bar u_j.
\]
Hybrid RANS/LES models introduce a modelled length scale and require careful wall-distance and grid-dependence analysis.

**Repository:** [turbulence.h](../../src/cfdx/physics/turbulence.h), [turbulence_models.h](../../src/cfdx/physics/turbulence_models.h), [turbulence_transport.h](../../src/cfdx/physics/turbulence_transport.h), [turbulence_solver.h](../../src/cfdx/physics/turbulence_solver.h), [sst_solver.h](../../src/cfdx/physics/sst_solver.h), [spalart_allmaras.h](../../src/cfdx/physics/spalart_allmaras.h), [wall_distance.h](../../src/cfdx/physics/wall_distance.h); tests: [test_turbulence_qualification_equations.cpp](../../tests/unit/test_turbulence_qualification_equations.cpp), [test_sst_blending.cpp](../../tests/unit/test_sst_blending.cpp).

---

# 11. Heat transfer

The thermal equation can be written
\[
\rho c_p\frac{DT}{Dt}
=
\nabla\cdot(k\nabla T)+S_T
\]
for a simple constant-property model.

Fourier's law:
\[
\mathbf q=-k\nabla T.
\]
The normal heat flux is
\[
q_n=\mathbf q\cdot\mathbf n=-k\nabla T\cdot\mathbf n.
\]
At a conjugate interface,
\[
T_f^{(1)}=T_f^{(2)},
\qquad
\mathbf q^{(1)}\cdot\mathbf n
=
-\mathbf q^{(2)}\cdot\mathbf n
\]
for perfect thermal contact.

The Biot and Fourier numbers are
\[
Bi=\frac{hL_c}{k_s},
\qquad
Fo=\frac{\alpha t}{L_c^2}.
\]
The Péclet number for thermal transport is
\[
Pe=\frac{UL}{\alpha}.
\]

**Repository:** [thermal.h](../../src/cfdx/physics/thermal.h), [energy_solver.h](../../src/cfdx/physics/energy_solver.h), [cht_solver.h](../../src/cfdx/physics/cht_solver.h), [conductivity.h](../../src/cfdx/transport/conductivity/conductivity.h); tests: [test_thermal_vv.cpp](../../tests/validation/test_thermal_vv.cpp), [test_cht_validation.cpp](../../tests/validation/test_cht_validation.cpp), [test_thermal_radiation_model_matrix.cpp](../../tests/validation/test_thermal_radiation_model_matrix.cpp).

---

# 12. Radiation

A blackbody emits
\[
E_b=\sigma T^4.
\]
For a diffuse-gray surface,
\[
E=\varepsilon\sigma T^4.
\]
For surface-to-surface exchange, radiosity is
\[
J=E+(1-\varepsilon)G,
\]
where G is irradiation. The net radiative heat flux is
\[
q''=J-G.
\]

For view factors,
\[
A_iF_{ij}=A_jF_{ji},
\]
and for a closed enclosure,
\[
\sum_jF_{ij}=1.
\]
These are exact geometric/physical identities and are natural verification targets.

For P1 radiation, a representative diffusion-type equation is
\[
-\nabla\cdot
\left(\frac{1}{3\beta}\nabla G\right)
+\beta G
=
4\beta\sigma T^4
\]
for the appropriate absorption/scattering convention. DOM instead solves directional transport equations for discrete ordinates.

**Repository:** [radiation.h](../../src/cfdx/physics/radiation.h), [radiation_s2s.h](../../src/cfdx/physics/radiation_s2s.h), [radiation_solver.h](../../src/cfdx/physics/radiation_solver.h), [radiation_models.h](../../src/cfdx/physics/radiation_models.h); tests: [test_radiation_vv.cpp](../../tests/validation/test_radiation_vv.cpp), [test_s2s_radiation_vv.cpp](../../tests/validation/test_s2s_radiation_vv.cpp), [test_phase12_radiation_vv.cpp](../../tests/validation/test_phase12_radiation_vv.cpp).

---

# 13. Multiphysics

A coupled problem is represented abstractly as
\[
F(y)=0,
\qquad
y=(u,p,T,Y_1,\ldots).
\]
Newton linearisation gives
\[
J(y^k)\delta y=-F(y^k),
\qquad
y^{k+1}=y^k+\delta y.
\]
The Jacobian has block structure:
\[
J=
\begin{bmatrix}
J_{uu}&J_{up}&J_{uT}&\cdots\\
J_{pu}&J_{pp}&J_{pT}&\cdots\\
J_{Tu}&J_{Tp}&J_{TT}&\cdots
\end{bmatrix}.
\]
Segregated algorithms approximate these couplings by solving subsets sequentially; monolithic methods solve the coupled block system.

At every interface, conservation must be checked separately from convergence:
\[
\mathcal R_\Gamma=
\sum_{f\in\Gamma}
(q_f^{(1)}+q_f^{(2)})
\]
must approach the declared discretisation/round-off tolerance for an internal conservative interface.

**Repository:** [finite_volume_transport.h](../../src/cfdx/physics/finite_volume_transport.h), [cht_solver.h](../../src/cfdx/physics/cht_solver.h), [radiation_solver.h](../../src/cfdx/physics/radiation_solver.h), [steady_incompressible_solver.h](../../src/cfdx/physics/steady_incompressible_solver.h); tests: [test_level_c_coupled_verification.cpp](../../tests/validation/test_level_c_coupled_verification.cpp), [test_cht_validation.cpp](../../tests/validation/test_cht_validation.cpp).

---

# 14. Numerical analysis

Consistency measures whether the discrete operator approaches the differential operator:
\[
L_hu=L u+\tau_h,
\qquad
\tau_h\rightarrow0
\quad(h\rightarrow0).
\]
For a smooth solution,
\[
u_h=u+Ch^p+o(h^p).
\]
Given two refinement levels,
\[
p_{obs}
=
\frac{\ln(E_h/E_{h/2})}{\ln2}.
\]
Stability concerns amplification of perturbations. For a linear evolution system,
\[
\frac{dU}{dt}=AU,
\]
stability requires bounded behaviour of the numerical amplification matrix over the relevant timestep range.

Boundedness means, for a scalar field with neighbour bounds,
\[
\min_N\phi_N\le\phi_f\le\max_N\phi_N
\]
when the scheme declares such a discrete maximum-principle property. Monotonicity is stronger and depends on the full coefficient structure.

The total error budget can be decomposed conceptually as
\[
E_{total}
\lesssim
E_{model}+E_{discretisation}+E_{iteration}
+E_{roundoff}+E_{input}.
\]
Reducing a linear residual does not remove modelling or discretisation error.

**Repository:** numerical operators under [core/numerics](../../src/cfdx/core/numerics/), linear algebra under [core/linalg](../../src/cfdx/core/linalg/); tests: [test_mesh_refinement_order.cpp](../../tests/validation/test_mesh_refinement_order.cpp), [test_mms_scalar_diffusion.cpp](../../tests/validation/test_mms_scalar_diffusion.cpp), [test_conservation_boundedness.cpp](../../tests/unit/test_conservation_boundedness.cpp).

---

# 15. Verification and validation

Verification asks whether the equations were solved correctly; validation asks whether the equations/model represent the physical system adequately.

For grid convergence, if
\[
E_h\approx Ch^p,
\]
then three systematically refined grids allow an observed order and, with appropriate assumptions, a Richardson/GCI-style uncertainty estimate. A generic GCI expression is
\[
GCI_{12}
=
F_s
\frac{|\phi_1-\phi_2|}
{|\phi_1|}
\frac{1}{r^p-1},
\]
with safety factor \(F_s\), refinement ratio r and observed/order-assumed p explicitly documented.

Manufactured solutions choose an analytic field u and derive the source term:
\[
S_{MMS}=L(u_{exact}).
\]
The numerical solution can then be compared directly with u_exact.

Iterative error must be separated from discretisation error. If
\[
\|r_k\|\rightarrow0
\]
but the grid error remains large, solver convergence has not established solution accuracy.

CFDX evidence should therefore contain:
- analytical identities;
- unit tests;
- MMS;
- refinement studies;
- benchmark comparisons;
- conservation diagnostics;
- retained numerical settings and revision identifiers.

**Repository:** [V&V tests](../../tests/validation/), [validation documents](../validation/), [qualification registry](../validation/CFDX_QUALIFICATION_REGISTRY.json).

Known engineering qualification cases include Couette, Poiseuille, Ghia, thermal/radiation matrices and selected longer campaigns. The documentation must report their actual current status rather than infer maturity from the existence of a test file.

---

# 16. Data model and file formats

The case representation separates setup from runtime state.

Conceptually:
\[
\mathcal C=
\{
\mathcal M,\mathcal B,\mathcal P,\mathcal N,\mathcal F,\mathcal V
\},
\]
where mesh, boundaries, physical models, numerics, fields and metadata form the authoritative setup.

HDF5 provides hierarchical datasets and attributes. A dataset can be viewed abstractly as
\[
D=(shape,dtype,data,attributes).
\]
The schema must define:
- coordinate units;
- topology;
- field locations;
- component ordering;
- boundary identifiers;
- numerical-method identifiers;
- provenance;
- version/compatibility semantics.

A checkpoint is state, not a replacement for the case definition. Restart correctness requires
\[
\mathcal S_{restart}(t_n)
\equiv
\mathcal S_{reference}(t_n)
\]
within a declared numerical comparison tolerance, including all required fields and history.

VTU output is a visualization representation and must not silently become the authoritative case state.

**Repository:** [case_hdf5_io.h](../../src/cfdx/io/hdf5/case_hdf5_io.h), [hdf5_reader.h](../../src/cfdx/io/hdf5/hdf5_reader.h), [hdf5_writer.h](../../src/cfdx/io/hdf5/hdf5_writer.h), [schema.h](../../src/cfdx/io/hdf5/schema.h), [dat_restart.h](../../src/cfdx/io/restart/dat_restart.h), [vtu_writer.h](../../src/cfdx/io/vtu/vtu_writer.h); tests: [test_case_hdf5_io.cpp](../../tests/unit/test_case_hdf5_io.cpp), [test_hdf5_roundtrip.cpp](../../tests/unit/test_hdf5_roundtrip.cpp), [test_solver_vtu_provenance.cpp](../../tests/validation/test_solver_vtu_provenance.cpp).

---

# 17. Computational chain

The CFDX execution chain is:

\[
\boxed{
\text{case}
\rightarrow
\text{mesh}
\rightarrow
\text{fields/BC}
\rightarrow
\text{numerics}
\rightarrow
\text{discretisation}
\rightarrow
\text{assembly}
\rightarrow
\text{linear/nonlinear solve}
\rightarrow
\text{diagnostics}
\rightarrow
\text{output/restart}
}
\]

A production run should be explainable at every transition.

### Case ingestion

The authoritative setup is parsed and validated before solver construction.

### Mesh

Topology is loaded, geometric quantities are computed, and invalid cells/faces are rejected.

### Fields and boundaries

Fields acquire dimensions, locations and boundary constraints. Boundary conditions contribute explicit algebraic terms.

### Numerical selection

The numerical-method registry resolves named algorithms to explicit implementations/contracts.

**Repository:** [numerical_method_registry.h](../../src/cfdx/core/numerics/numerical_method_registry.h), [numerical_method_contract.h](../../src/cfdx/core/numerics/numerical_method_contract.h).

### Assembly

The residual/operator is built from face fluxes, source terms and temporal contributions.

### Solver

The resulting algebraic systems are solved using the selected linear and pressure-velocity algorithms.

### Diagnostics

Diagnostics include residuals, conservation defects, boundedness, solver status, iteration counts and provenance.

### Output

Fields are written to state/checkpoint or visualization formats without changing their semantic role.

**Repository:** [simulation_controller.cpp](../../src/cfdx/application/simulation_controller.cpp), [production solver](../../apps/production_solver/main.cpp), [finite_volume_transport.h](../../src/cfdx/physics/finite_volume_transport.h).

---

# 18. Source-tree physics audit and traceability

The documentation must be auditable from equation to code.

The required chain for a numerical feature is:

\[
\text{equation}
\rightarrow
\text{discrete operator}
\rightarrow
\text{implementation}
\rightarrow
\text{unit test}
\rightarrow
\text{verification campaign}
\rightarrow
\text{validation/qualification evidence}.
\]

A source path is not evidence of correctness. A passing unit test is not automatically physical validation. A benchmark result is not automatically qualification for every mesh/physics population.

## 18.1 Core map

| Domain | Primary source paths |
|---|---|
| Mesh topology | [src/cfdx/core/mesh](../../src/cfdx/core/mesh) |
| Geometry | [src/cfdx/core/geometry](../../src/cfdx/core/geometry) |
| Fields | [src/cfdx/core/field](../../src/cfdx/core/field) |
| Boundary conditions | [src/cfdx/core/boundary](../../src/cfdx/core/boundary) |
| Interpolation/gradients | [src/cfdx/core/numerics](../../src/cfdx/core/numerics) |
| FVM transport | [src/cfdx/physics/finite_volume_transport.h](../../src/cfdx/physics/finite_volume_transport.h) |
| Incompressible flow | [src/cfdx/physics/incompressible.h](../../src/cfdx/physics/incompressible.h) |
| Pressure-velocity coupling | [src/cfdx/physics/pressure_velocity_algorithms.h](../../src/cfdx/physics/pressure_velocity_algorithms.h) |
| Turbulence | [src/cfdx/physics/turbulence_models.h](../../src/cfdx/physics/turbulence_models.h) |
| Thermal | [src/cfdx/physics/thermal.h](../../src/cfdx/physics/thermal.h) |
| Radiation | [src/cfdx/physics/radiation.h](../../src/cfdx/physics/radiation.h) |
| Linear algebra | [src/cfdx/core/linalg](../../src/cfdx/core/linalg) |
| HDF5 | [src/cfdx/io/hdf5](../../src/cfdx/io/hdf5) |
| Restart | [src/cfdx/io/restart/dat_restart.h](../../src/cfdx/io/restart/dat_restart.h) |
| VTU | [src/cfdx/io/vtu/vtu_writer.h](../../src/cfdx/io/vtu/vtu_writer.h) |
| Production execution | [apps/production_solver/main.cpp](../../apps/production_solver/main.cpp) |
| Unit tests | [tests/unit](../../tests/unit) |
| Numerical tests | [tests/numerical](../../tests/numerical) |
| Validation tests | [tests/validation](../../tests/validation) |

## 18.2 Audit rules

1. Every equation has a defined variable convention.
2. Every implementation claim points to an existing repository path.
3. Every verification claim points to a test/campaign.
4. Validation claims identify the independent reference.
5. Qualification claims identify the population and acceptance gate.
6. Partial functionality remains explicitly marked partial.
7. A fallback algorithm is never described as equivalent to the requested method.
8. Numerical tolerances are not relaxed merely to obtain a passing status.
9. Failed or incomplete benchmarks remain visible in the evidence record.
10. Documentation is updated from repository evidence, not from intended architecture.

## 18.3 N2 traceability example

\[
\nabla\phi
\rightarrow
\begin{cases}
\text{Green--Gauss},\\
\text{LS/WLS},\\
\text{vertex reconstruction}
\end{cases}
\rightarrow
\text{face reconstruction}
\rightarrow
\text{diffusion/viscous terms}
\rightarrow
\text{verification}.
\]

The relevant source files are [gradient.h](../../src/cfdx/core/numerics/gradient.h), [gradient_stencil.h](../../src/cfdx/core/numerics/gradient_stencil.h), [least_squares_gradient.h](../../src/cfdx/core/fvm/least_squares_gradient.h) and [interpolation.h](../../src/cfdx/core/numerics/interpolation.h). The corresponding campaign includes [test_gradient_verification.cpp](../../tests/validation/test_gradient_verification.cpp) and [test_polyhedral_gradient_campaign.cpp](../../tests/validation/test_polyhedral_gradient_campaign.cpp).

The scientific acceptance question is not merely whether these files compile. It is whether the selected algorithms reproduce the expected polynomial exactness, convergence order, robustness under distorted/polyhedral geometry, conditioning behaviour, boundary reconstruction and face-value reconstruction required by N2.

---

# Cross-chapter numerical vocabulary

**Consistency**
\[
L_hu-Lu\rightarrow0.
\]

**Stability**
\[
\|G^n\|\le C
\]
over the declared time/resolution range for an appropriate amplification operator G.

**Convergence**
\[
u_h\rightarrow u
\quad\text{as}\quad h\rightarrow0.
\]

**Observed order**
\[
p_{obs}=\frac{\ln(E_h/E_{rh})}{\ln r}.
\]

**Conservation defect**
\[
\epsilon_C=
\frac{|\text{inflow}-\text{outflow}+\text{source}-dQ/dt|}
{Q_{ref}}.
\]

**Residual**
\[
r=b-Ax.
\]

**Energy norm**
\[
\|e\|_A=\sqrt{e^TAe}.
\]

**Condition number**
\[
\kappa(A)=\|A\|\|A^{-1}\|.
\]

These quantities are complementary. None is a universal substitute for the others.

# References

The theory should cite stable editions/identifiers in the central bibliography. Core references include Versteeg & Malalasekera, *An Introduction to Computational Fluid Dynamics*, 2nd ed.; Ferziger, Perić & Street, *Computational Methods for Fluid Dynamics*; Moukalled, Mangani & Darwish, *The Finite Volume Method in Computational Fluid Dynamics*; Jasak, *Error Analysis and Estimation for the Finite Volume Method with Applications to Fluid Flows*; Barth & Jespersen on multidimensional limiter construction; and the official OpenFOAM gradient documentation for implementation-level comparison. OpenFOAM currently documents Gauss, least-squares and limited gradient families explicitly, including weighted least-squares stencil infrastructure in its API. citeturn0search0turn0search4turn0search8

**Important status statement:** this document is a theory and traceability reference. It does not declare N1–N17 complete. In particular, N2 remains **PARTIAL** until its repository-specific verification and qualification gates are actually satisfied.
