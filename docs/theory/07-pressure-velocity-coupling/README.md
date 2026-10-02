# 07 — Pressure–Velocity Coupling

For incompressible flow, pressure enforces the divergence constraint rather than being supplied by an independent thermodynamic equation.

## 1. Block system

Linearisation gives
\[
\begin{bmatrix}
A & G\\
D & 0
\end{bmatrix}
\begin{bmatrix}
u\\p
\end{bmatrix}
=
\begin{bmatrix}
b\\0
\end{bmatrix}.
\]
Here \(A\) is the momentum operator, \(G\) pressure gradient and \(D\) discrete divergence.

## 2. Schur complement

From
\[
Au+Gp=b,
\qquad
Du=0,
\]
we obtain
\[
u=A^{-1}(b-Gp),
\]
then
\[
(DA^{-1}G)p=DA^{-1}b.
\]
The Schur operator is
\[
S=DA^{-1}G.
\]
Its quality controls pressure correction and coupled convergence.

**CFDX paths:** `src/cfdx/physics/pressure_velocity.h`, `src/cfdx/core/linalg/exact_schur.h`, `src/cfdx/core/linalg/block_schur.h`.

## 3. SIMPLE

With approximate momentum inverse \(\tilde A^{-1}\),
\[
\tilde S=D\tilde A^{-1}G.
\]
A pressure correction is obtained from the approximate Schur system and relaxed before updating velocity and pressure. The approximation makes SIMPLE inexpensive but iterative.

## 4. SIMPLEC

SIMPLEC modifies the momentum-correction approximation to reduce the influence of neighbour velocity corrections. The exact algebra must follow the CFDX implementation rather than a generic label.

**CFDX path:** `src/cfdx/core/linalg/simplerc_schur.h`.

## 5. PISO

PISO performs multiple pressure-correction stages within a time step so that the corrected velocity satisfies continuity more tightly without a complete outer nonlinear iteration.

## 6. PIMPLE

PIMPLE combines pressure-correction iterations with outer under-relaxed/nonlinear loops. It is particularly useful for transient problems where multiple pressure corrections and outer iterations are required.

## 7. Fractional step

A projection method first computes an intermediate velocity,
\[
\frac{u^*-u^n}{\Delta t}=N(u^n)-\nabla p^*,
\]
then corrects it:
\[
u^{n+1}=u^*-\Delta t\nabla\delta p,
\]
with
\[
\nabla^2\delta p=\frac{1}{\Delta t}\nabla\cdot u^*.
\]

## 8. Coupled solve

A monolithic approach solves the saddle-point system directly or through block preconditioning. This avoids some splitting errors but requires robust block linear algebra.

## 9. Rhie–Chow and checkerboarding

On collocated grids, pressure and velocity can decouple through a checkerboard mode. Face mass flux reconstruction must therefore couple pressure and momentum consistently. Any Rhie–Chow-like correction must be documented algebraically and verified on a pressure-mode test.

## 10. Pressure null space

Pressure is defined only up to an additive constant for closed incompressible domains:
\[
p'=p+C.
\]
A gauge constraint such as \(\int_\Omega p\,d\Omega=0\) or one reference value removes the null space.

## 11. Verification

Use discrete divergence, pressure null-space, manufactured pressure-gradient, lid-driven cavity, Poiseuille and coupling-iteration diagnostics.

**Implementation:** `src/cfdx/physics/pressure_velocity.h`, `src/cfdx/physics/pressure_velocity_algorithms.h`, `src/cfdx/core/linalg/`.
