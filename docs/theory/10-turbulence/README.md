# 10 — Turbulence

## 1. Reynolds decomposition
\[
u_i=\bar u_i+u_i',\qquad\overline{u_i'}=0.
\]
Averaging creates
\[
-\rho\overline{u_i'u_j'},
\]
the Reynolds stress.

## 2. Eddy viscosity
The Boussinesq hypothesis is
\[
-\rho\overline{u_i'u_j'}
=2\mu_tS_{ij}-\frac23\rho k\delta_{ij},
\]
with
\[
S_{ij}=\frac12
\left(
\frac{\partial\bar u_i}{\partial x_j}
+
\frac{\partial\bar u_j}{\partial x_i}
\right),
\qquad
k=\frac12\overline{u_i'u_i'}.
\]
The model therefore turns unresolved stress into a closure for \(\mu_t\).

## 3. Generic transport equation
\[
\frac{\partial(\rho\phi)}{\partial t}
+\nabla\cdot(\rho u\phi)
=
\nabla\cdot(\Gamma_\phi\nabla\phi)+S_\phi.
\]
Production, destruction and cross-diffusion terms must be dimensionally and sign-consistently implemented.

## 4. k–epsilon / k–omega / SST
The two-equation family transports turbulence scales. A representative production is
\[
P_k=2\mu_tS:S.
\]
The exact constants, blending functions, compressibility corrections and limiters are model-specific and must be read from the CFDX implementation.

## 5. Spalart–Allmaras
SA transports a modified viscosity \(\tilde\nu\). Its wall-distance dependence makes
\[
d_w=\text{distance to wall}
\]
part of the numerical model. Wall-distance errors therefore propagate into turbulence production/destruction.

## 6. LES
Filtering creates subgrid stress
\[
\tau_{ij}=\overline{u_iu_j}-\bar u_i\bar u_j.
\]
Smagorinsky gives
\[
\nu_t=(C_s\Delta)^2|S|,
\qquad
|S|=\sqrt{2S:S}.
\]

## 7. Wall treatment
\[
u_t|_w=0,\qquad
\tau_w=\mu\left.\frac{\partial u_t}{\partial n}\right|_w,
\]
and
\[
u^+=u_t/u_\tau,\qquad
y^+=\rho u_\tau y/\mu.
\]
These nondimensional coordinates determine wall-resolution requirements.

**CFDX paths:** `src/cfdx/physics/turbulence.h`, `turbulence_transport.h`, `wall_distance_turbulence.h`.

## 8. Verification
Verify transport operators and source terms independently, then manufactured turbulence fields, positivity, wall-distance sensitivity and canonical validation cases. See `docs/validation/TURBULENCE_CLOSURE_EQUATIONS.md`.