# 13 — Multiphysics

## 1. Residual formulation
All coupled variables are collected into
\[
y=[u,p,T,\phi_1,\ldots]^T,
\qquad
F(y)=0.
\]
Newton linearisation gives
\[
J(y^k)\delta y=-F(y^k),
\qquad
y^{k+1}=y^k+\delta y.
\]

## 2. Block Jacobian
\[
J=
\begin{bmatrix}
J_{uu}&J_{up}&J_{uT}\\
J_{pu}&J_{pp}&J_{pT}\\
J_{Tu}&J_{Tp}&J_{TT}
\end{bmatrix}.
\]
Diagonal blocks represent self-physics; off-diagonal blocks represent coupling. Omitting a required off-diagonal derivative can alter nonlinear convergence or even the solution path.

## 3. Segregated solution
A block Gauss–Seidel iteration has the form
\[
u^{k+1}=F_u^{-1}(p^k,T^k),
\]
\[
p^{k+1}=F_p^{-1}(u^{k+1},T^k),
\]
\[
T^{k+1}=F_T^{-1}(u^{k+1},p^{k+1}).
\]
Under-relaxation is
\[
y^{k+1}=(1-\omega)y^k+\omega y^{k+1}_{raw}.
\]

## 4. Monolithic solution
\[
J\delta y=-F
\]
solves all variables simultaneously. Block preconditioning is then essential for scalability.

**CFDX paths:** `src/cfdx/physics/`, `src/cfdx/core/linalg/`.

## 5. Interface conservation
A conservative interface satisfies
\[
q_A+q_B=0.
\]
For thermal perfect contact,
\[
T_A=T_B.
\]
Equivalent pairwise balance checks should be applied to momentum, heat and radiation exchange.

## 6. Convergence
Monitor
\[
\eta_F=\frac{\|F(y^k)\|}{\|F(y^0)\|},
\qquad
\eta_y=\frac{\|\delta y\|}
{\max(\|y^k\|,y_{scale})}.
\]
A linear solver convergence flag cannot substitute for the nonlinear physical residual.

## 7. Verification ladder
First verify each physics alone, then each pairwise coupling, then the full coupled system. Conservation checks must be retained at every level.