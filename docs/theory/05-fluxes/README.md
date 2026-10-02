# 05 — Fluxes

Fluxes are the numerical interface between neighbouring control volumes and therefore the primary carriers of conservation.

## 1. Mass flux

\[
F_m=\rho_f\mathbf u_f\cdot\mathbf S_f.
\]
For incompressible constant density, mass conservation requires
\[
\sum_fF_m=0.
\]

## 2. Convective flux

For transported scalar \(\phi\),
\[
F_\phi=F_m\phi_f.
\]
The numerical question is therefore how \(\phi_f\) is reconstructed from cell values. This separates the conservative flux from the reconstruction scheme.

## 3. Diffusive flux

\[
F_d=\Gamma_f(\nabla\phi)_f\cdot\mathbf S_f.
\]
The sign convention is outward flux. Fourier heat flux is the special case \(\Gamma=k\) with \(\mathbf q=-k\nabla T\).

## 4. Pressure flux

The pressure force on a cell is
\[
\mathbf F_p=-\sum_f p_f\mathbf S_f.
\]
The minus sign follows from \(\boldsymbol\sigma_p=-p\mathbf I\). Pressure interpolation and pressure-gradient reconstruction must preserve the intended momentum balance.

## 5. Energy flux

The total-energy advective flux is
\[
F_E=(\rho E+p)\mathbf u\cdot\mathbf S,
\]
while conductive and viscous contributions are added consistently with the energy equation.

## 6. Face conservation

For every internal face,
\[
F_{P,f}=-F_{N,f}.
\]
This must be true for the numerical flux itself, not merely after summing cells.

## 7. Bounded convection

A high-order reconstruction may reduce truncation error but generate overshoots. Limiter functions commonly use
\[
r=\frac{\Delta^-\phi}{\Delta^+\phi},
\qquad
\phi_f=\phi_{up}+\frac12\Psi(r)(\phi_{down}-\phi_{up}),
\]
with a limiter \(\Psi\) designed to satisfy a specified boundedness/TV​D property.

**CFDX paths:** `src/cfdx/core/numerics/flux.h`, `src/cfdx/core/numerics/convection.h`.

## 8. Verification

Use constant fields, linear manufactured fields, internal-face antisymmetry, diffusion of a known linear field, conservation and refinement. Never infer flux correctness solely from solver convergence.

**Tests:** `tests/unit/test_interpolation.cpp` and numerical-method tests under `tests/`.
