# 02 — Finite-Volume Method

## 1. Principle

FVM integrates the PDE over each control volume. For a scalar transport equation
\[
\frac{\partial(\rho\phi)}{\partial t}+\nabla\cdot(\rho\mathbf u\phi)
=\nabla\cdot(\Gamma\nabla\phi)+S_\phi,
\]
integration gives
\[
\frac{d}{dt}\int_V\rho\phi\,dV+
\sum_fF_f\phi_f
=
\sum_f\Gamma_f(\nabla\phi)_f\cdot\mathbf S_f+
S_\phi V.
\]

## 2. Geometry

For face \(f\),
\[
\mathbf S_f=A_f\mathbf n_f.
\]
For a closed polyhedral cell,
\[
\sum_f\mathbf S_f=\mathbf0.
\]
The cell volume may be obtained from an oriented surface representation; its positivity and consistency with the face orientations are mandatory geometric invariants.

**CFDX paths:** `src/cfdx/core/mesh/mesh.h`, `src/cfdx/core/mesh/mesh.cpp`, `src/cfdx/core/mesh/solver_mesh.h`.

## 3. Face interpolation

Linear interpolation between cell centres is
\[
\phi_f=(1-w)\phi_P+w\phi_N,
\qquad
w=\frac{|\mathbf C_f-\mathbf C_P|}
{|\mathbf C_N-\mathbf C_P|}.
\]
For a midpoint, \(w=1/2\). On skewed/non-orthogonal meshes, the geometric face point and the line joining cell centres need not coincide, so the interpolation policy must be explicit.

**CFDX path:** `src/cfdx/core/numerics/interpolation.h`.

## 4. Diffusion

Orthogonal diffusion is approximated by
\[
(\nabla\phi)_f\cdot\mathbf S_f
\approx
\frac{\phi_N-\phi_P}{d_{PN}}A_f.
\]
For non-orthogonal geometry,
\[
\mathbf S_f=\mathbf d_f+\mathbf S_f^{\perp},
\]
and
\[
(\nabla\phi)_f\cdot\mathbf S_f
=
(\nabla\phi)_f\cdot\mathbf d_f+
(\nabla\phi)_f\cdot\mathbf S_f^{\perp}.
\]
The correction requires a gradient/reconstruction method consistent with the documented policy.

**CFDX path:** `src/cfdx/core/numerics/laplacian.h`.

## 5. Convection

Define
\[
F_f=\rho_f\mathbf u_f\cdot\mathbf S_f.
\]
The convective contribution is \(F_f\phi_f\). Upwind uses the sign of \(F_f\); central differencing uses a centred interpolation; higher-order schemes require boundedness analysis.

**CFDX path:** `src/cfdx/core/numerics/convection.h`.

## 6. Sources

A linearised source is commonly written
\[
S_\phi=S_C+S_P\phi_P.
\]
The discrete contribution changes the diagonal and RHS. A source linearisation is numerically safer when \(S_P\le0\), because it increases diagonal dominance rather than destabilising it.

## 7. Matrix structure

After discretisation,
\[
a_P\phi_P-\sum_Na_{PN}\phi_N=b_P.
\]
The signs and diagonal dominance depend on the PDE and discretisation. One must inspect the actual assembled matrix rather than infer properties from a generic stencil.

## 8. Boundedness

A sufficient local form for a scalar equation is often an M-matrix-like structure:
\[
a_P\ge\sum_Na_{PN}\ge0,
\]
with appropriate RHS/source conditions. These are not universal necessary conditions, but they are useful diagnostics for diffusion-dominated transport.

## 9. Polyhedral FVM

For arbitrary polyhedra, all faces participate in the same balance:
\[
\sum_fF_f=0
\]
for a closed stationary incompressible cell. No hexahedral assumption is allowed in the algebraic formulation.

## 10. Verification

Check geometry closure, constant-field gradient/diffusion behaviour, face antisymmetry, matrix assembly against a hand-computed cell, conservation, boundedness and mesh refinement.

**Implementation:** `src/cfdx/core/numerics/`, `src/cfdx/core/fvm/`, `tests/unit/test_interpolation.cpp`, `tests/unit/test_laplacian.cpp`.
