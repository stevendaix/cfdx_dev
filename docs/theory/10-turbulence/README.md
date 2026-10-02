# 10 — Turbulence

## 1. Reynolds decomposition

Write
\[
u_i=\overline u_i+u_i',
\qquad
\overline{u_i'}=0.
\]
Averaging the Navier–Stokes equations produces the Reynolds-stress term
\[
-\rho\overline{u_i'u_j'},
\]
which creates the closure problem.

## 2. Eddy viscosity

The Boussinesq hypothesis writes
\[
-\rho\overline{u_i'u_j'}
=
2\mu_t S_{ij}
-\frac23\rho k\delta_{ij},
\]
where
\[
S_{ij}=\frac12
\left(
\frac{\partial\overline u_i}{\partial x_j}
+
\frac{\partial\overline u_j}{\partial x_i}
\right),
\quad
k=\frac12\overline{u_i'u_i'}.
\]

## 3. RANS closure

The turbulent viscosity \(\mu_t\) is supplied by a model. A generic transported quantity satisfies
\[
\frac{\partial(\rho\phi)}{\partial t}
+\nabla\cdot(\rho\mathbf u\phi)
=
\nabla\cdot(\Gamma_\phi\nabla\phi)+S_\phi.
\]

## 4. SST and k-omega

The two-equation family solves transport equations for \(k\) and a frequency-like variable \(\omega\). A representative production term is
\[
P_k=2\mu_t S:S.
\]
The exact source, cross-diffusion, limiter and blending functions must be taken from the CFDX implementation rather than inferred from the model name.

## 5. Spalart–Allmaras

SA uses a transported modified viscosity variable \(\tilde\nu\). Its model is sensitive to wall distance, destruction terms and production limiting. Therefore turbulence verification cannot be separated from the wall-distance contract.

## 6. LES

LES resolves large scales and models subgrid stresses:
\[
\tau_{ij}^{sgs}=\overline{u_i u_j}-\bar u_i\bar u_j.
\]
A Smagorinsky-type closure is
\[
\nu_t=(C_s\Delta)^2|S|,
\qquad
|S|=\sqrt{2S:S}.
\]

## 7. DES/DDES/IDDES

Hybrid RANS/LES models modify the length scale so that attached boundary layers retain RANS behaviour while separated regions can resolve larger eddies. Their validity depends strongly on mesh, wall treatment and shielding functions.

## 8. Wall treatment

Near a no-slip wall,
\[
u_t=0.
\]
The wall shear stress is
\[
\tau_w=\mu\left.\frac{\partial u_t}{\partial n}\right|_w,
\]
and
\[
u^+=\frac{u_t}{u_\tau},
\qquad
y^+=\frac{\rho u_\tau y}{\mu}.
\]
These quantities determine whether a wall-resolved or wall-function regime is being represented.

**CFDX paths:** `src/cfdx/physics/turbulence.h`, `src/cfdx/physics/turbulence_transport.h`, `src/cfdx/physics/wall_distance_turbulence.h`.

## 9. Verification

Verify source-term signs, positivity, dimensions, wall-distance dependence, limiting behaviour and manufactured transport solutions before physical validation. See `docs/validation/TURBULENCE_CLOSURE_EQUATIONS.md`.
