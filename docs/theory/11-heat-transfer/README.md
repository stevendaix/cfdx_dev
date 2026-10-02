# 11 — Heat Transfer

## 1. Energy balance
For a constant-property fluid,
\[
\rho c_p\left(\frac{\partial T}{\partial t}+u\cdot\nabla T\right)
=\nabla\cdot(k\nabla T)+\dot q_v.
\]
For a stationary solid,
\[
\nabla\cdot(k\nabla T)+\dot q_v=0.
\]

## 2. Fourier conduction
\[
q_c=-k\nabla T,
\qquad
Q_f=q_c\cdot S_f.
\]
The sign is physical: heat flows down the temperature gradient.

**CFDX:** `src/cfdx/physics/thermal.h`, `src/cfdx/core/numerics/laplacian.h`.

## 3. Convection
A Newton boundary condition is
\[
q_n=h(T_s-T_\infty).
\]
The coefficient \(h\) is a model input/correlation, not an intrinsic constant of the wall.

## 4. Dimensionless groups
\[
\alpha=\frac{k}{\rho c_p},
\quad
Pe=\frac{UL}{\alpha},
\quad
Pr=\frac{\nu}{\alpha},
\quad
Nu=\frac{hL}{k}.
\]
The relative magnitude of advection and diffusion is controlled by \(Pe\).

## 5. Conjugate heat transfer
Fluid:
\[
\rho_fc_f\frac{DT_f}{Dt}=\nabla\cdot(k_f\nabla T_f)+S_f.
\]
Solid:
\[
\rho_sc_s\frac{\partial T_s}{\partial t}=\nabla\cdot(k_s\nabla T_s)+S_s.
\]
At a perfect interface,
\[
T_f=T_s,
\qquad
-k_f\nabla T_f\cdot n=-k_s\nabla T_s\cdot n.
\]
The second equation is conservation of interfacial heat flux.

## 6. Nonlinear properties
If \(k=k(T)\), then
\[
\nabla\cdot(k(T)\nabla T)
\]
contains both coefficient variation and temperature gradient. Picard and Newton treatments have different convergence properties and must be documented.

## 7. Verification
Use 1-D conduction with known solution, linear manufactured temperature, interface flux antisymmetry, energy conservation and mesh refinement. Regression evidence is indexed by `docs/validation/THERMAL_RADIATION_REGRESSION_MATRIX.md`.