# 01 — Conservation Laws

## 1. Generic balance

For an extensive quantity with density \(q\), flux \(\mathbf F\) and source \(s\),
\[
\frac{d}{dt}\int_Vq\,dV+\oint_{\partial V}\mathbf F\cdot\mathbf n\,dA=\int_Vs\,dV.
\]
Finite volume retains this integral form, which is why internal face fluxes can cancel exactly.

## 2. Mass

\[
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0.
\]
Integrating over cell \(P\):
\[
\frac{d(\rho_PV_P)}{dt}+\sum_{f\in\partial P}F_{m,f}=0,
\qquad
F_{m,f}=\rho_f\mathbf u_f\cdot\mathbf S_f.
\]
For an internal face shared by \(P\) and \(N\),
\[
\mathbf S_{f,N}=-\mathbf S_{f,P},
\]
so one physical flux is counted with opposite signs.

## 3. Momentum

\[
\frac{d}{dt}\int_V\rho\mathbf u\,dV+
\oint\rho\mathbf u(\mathbf u\cdot\mathbf n)\,dA
=
-\oint p\mathbf n\,dA+
\oint\boldsymbol\tau\mathbf n\,dA+
\int_V\rho\mathbf f\,dV.
\]
The pressure contribution is a surface force; viscous stress is also a surface traction.

**CFDX paths:** `src/cfdx/physics/incompressible.h`, `src/cfdx/core/numerics/flux.h`.

## 4. Energy

\[
\frac{d}{dt}\int_V\rho E\,dV+
\oint(\rho E+p)\mathbf u\cdot\mathbf n\,dA
=
\oint(\boldsymbol\tau\mathbf u-\mathbf q)\cdot\mathbf n\,dA+
\int_VS_E\,dV.
\]
With Fourier's law \(\mathbf q=-k\nabla T\), the conductive flux leaving a cell is \(-k\nabla T\cdot\mathbf n\).

**CFDX path:** `src/cfdx/physics/thermal.h`.

## 5. Species

For mass fraction \(Y_k\),
\[
\frac{\partial(\rho Y_k)}{\partial t}
+\nabla\cdot(\rho\mathbf uY_k)
=
-\nabla\cdot\mathbf J_k+\dot\omega_k,
\]
with
\[
\sum_kY_k=1,\qquad \sum_k\mathbf J_k=0
\]
for a consistent mixture model.

## 6. Conservation after discretisation

A cell equation has the form
\[
a_P\phi_P=\sum_Na_{PN}\phi_N+b_P.
\]
The coefficients arise from face fluxes and sources, not from an arbitrary algebraic fit. A global conservation diagnostic is
\[
R_{global}=\sum_P
\left[
\frac{d}{dt}(q_PV_P)+\sum_fF_{q,f}-S_PV_P
\right].
\]
For a closed steady system, \(R_{global}\) should vanish to the expected arithmetic/discretisation tolerance.

## 7. Flux antisymmetry

For an internal face,
\[
F_{P,f}+F_{N,f}=0.
\]
This is stronger than checking that a final global sum is small: two errors can cancel globally. CFDX verification should therefore inspect both per-face antisymmetry and global closure.

## 8. Boundary terms

Boundary faces have no neighbouring cell. Their contribution must be generated from the declared physical boundary condition. The sign is always defined using the outward normal of the current control volume.

**CFDX paths:** `src/cfdx/core/boundary/`, `src/cfdx/physics/boundary_constraint_fvm.h`.

## 9. Conservation and boundedness

Conservation and boundedness are different properties. A conservative discretisation can still produce an unphysical overshoot. Positivity of a scalar may require coefficient and limiter conditions in addition to exact flux cancellation.

## 10. Verification

Required evidence includes:
1. constant-field flux cancellation;
2. internal-face antisymmetry;
3. closed-domain conservation;
4. source-term balance;
5. refinement of conservation error where appropriate;
6. independent reconstruction of reported residuals.

**Implementation:** `src/cfdx/core/numerics/flux.h`, `src/cfdx/core/numerics/convection.h`, `src/cfdx/core/numerics/laplacian.h`.
