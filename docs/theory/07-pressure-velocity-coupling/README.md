# 07 — Pressure–Velocity Coupling

For incompressible flow, pressure is a Lagrange multiplier enforcing
\[
D u=0,
\]
not an independent thermodynamic variable.

## 1. Block equations
After discretisation,
\[
\begin{bmatrix}A&G\\D&0\end{bmatrix}
\begin{bmatrix}u\\p\end{bmatrix}
=
\begin{bmatrix}b\\0\end{bmatrix}.
\]
Here \(A\) contains momentum, \(G\) pressure gradient and \(D\) divergence.

## 2. Schur complement
From \(Au+Gp=b\),
\[
u=A^{-1}(b-Gp),
\]
and therefore
\[
DA^{-1}Gp=DA^{-1}b.
\]
The pressure operator is
\[
S=DA^{-1}G.
\]
An approximate pressure solver replaces \(A^{-1}\) by \(\tilde A^{-1}\).

**CFDX paths:** `src/cfdx/physics/pressure_velocity.h`, `src/cfdx/core/linalg/exact_schur.h`, `src/cfdx/core/linalg/block_schur.h`.

## 3. SIMPLE
Linearising momentum gives
\[
A\delta u=-G\delta p+r_u,
\qquad
D\delta u=-r_p.
\]
With \(A^{-1}\approx\tilde A^{-1}\),
\[
(D\tilde A^{-1}G)\delta p=
r_p+D\tilde A^{-1}r_u.
\]
Relaxation is then applied to pressure and/or velocity. SIMPLE convergence depends strongly on the Schur approximation.

## 4. SIMPLEC
SIMPLEC modifies the approximation to the velocity correction so that neighbour corrections are represented differently. The exact coefficient formula is implementation-defined and must be traced to `src/cfdx/core/linalg/simplerc_schur.h`.

## 5. PISO
PISO performs successive pressure corrections inside a time step. Each correction attempts to reduce the continuity defect without requiring a complete outer momentum solve.

## 6. PIMPLE
PIMPLE combines inner pressure corrections with outer nonlinear iterations. Its usefulness comes from separating fast pressure correction from slower nonlinear coupling.

## 7. Fractional step
A projection method computes
\[
\frac{u^*-u^n}{\Delta t}=N(u^n),
\]
then
\[
u^{n+1}=u^*-\Delta tG\delta p.
\]
Applying \(D\) and enforcing \(Du^{n+1}=0\) gives
\[
(DG)\delta p=\frac{1}{\Delta t}Du^*.
\]
Boundary and pressure-gauge treatment are part of this equation.

## 8. Collocated pressure/velocity
A collocated arrangement can admit checkerboard pressure modes. Face mass fluxes therefore need a pressure-velocity coupling correction consistent with the momentum equation. Any Rhie–Chow-like implementation must be tested against a manufactured checkerboard mode.

## 9. Pressure null space
For closed incompressible flow,
\[
p\mapsto p+C
\]
leaves the velocity unchanged. A gauge such as \(p(x_0)=0\) or \(\int p\,dV=0\) removes the null space.

## 10. Verification
Check discrete continuity, pressure gauge, Schur oracle, manufactured pressure gradient, Poiseuille and Ghia cavity. Solver residual reduction alone is not a coupling verification criterion.