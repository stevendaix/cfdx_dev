# 06 — Time Integration

## 1. Semi-discrete equation

After spatial discretisation,
\[
M\frac{d\mathbf u}{dt}=R(\mathbf u,t),
\]
or, after moving terms,
\[
\frac{d\mathbf u}{dt}=F(\mathbf u,t).
\]

**CFDX path:** `src/cfdx/core/numerics/temporal.h`.

## 2. Explicit Euler

\[
\mathbf u^{n+1}=\mathbf u^n+\Delta tF(\mathbf u^n,t^n).
\]
The local truncation error is \(O(\Delta t^2)\), giving first-order global accuracy. Stability depends on the eigenvalues of the spatial operator.

## 3. Implicit Euler

\[
\mathbf u^{n+1}=\mathbf u^n+\Delta tF(\mathbf u^{n+1},t^{n+1}).
\]
For a linear system \(F=-A\mathbf u+b\),
\[
(I+\Delta tA)\mathbf u^{n+1}=\mathbf u^n+\Delta t b.
\]
It is first-order accurate but often substantially more stable for stiff diffusion.

## 4. Crank–Nicolson

\[
\mathbf u^{n+1}=\mathbf u^n+
\frac{\Delta t}{2}
[F(\mathbf u^n,t^n)+F(\mathbf u^{n+1},t^{n+1})],
\]
with second-order temporal accuracy for sufficiently smooth solutions.

## 5. BDF2

For constant \(\Delta t\),
\[
\frac{3u^{n+1}-4u^n+u^{n-1}}{2\Delta t}=F(u^{n+1},t^{n+1}).
\]
It is second-order and implicit.

## 6. CFL

For advection, a representative Courant number is
\[
Co=\frac{U\Delta t}{\Delta x}.
\]
On a finite-volume mesh a cell-based definition is
\[
Co_P=\frac{\Delta t}{V_P}\sum_f\max(F_{m,f},0)/\rho_f,
\]
up to the exact convention used by the implementation. The convention must be documented rather than assumed.

## 7. Pseudo-transient continuation

A steady nonlinear problem
\[
R(u)=0
\]
may be solved through
\[
M\frac{u^{n+1}-u^n}{\Delta t}+R(u^{n+1})=0.
\]
The pseudo-time step is a nonlinear stabilisation parameter, not physical time.

## 8. Verification

Use manufactured ODE/PDE solutions, temporal refinement at fixed spatial resolution, order estimation
\[
p=\frac{\log(E_{\Delta t}/E_{\Delta t/2})}{\log2},
\]
and sensitivity to solver tolerances.

**Implementation:** `src/cfdx/core/numerics/temporal.h`, with time history/checkpoint handling in `src/cfdx/application/`.
