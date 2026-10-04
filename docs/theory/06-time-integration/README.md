# 06 — Time Integration

After spatial discretisation the semi-discrete problem is
\[
M\dot u=R(u,t).
\]
For a linear diffusion/advection operator, \(\dot u=Lu+f\). The temporal scheme must be analysed separately from spatial accuracy.

## 1. Explicit Euler
\[
u^{n+1}=u^n+\Delta t,F(u^n,t^n).
\]
Taylor expansion gives
\[
u(t+\Delta t)=u+\Delta t\dot u+\frac{\Delta t^2}{2}\ddot u+O(\Delta t^3),
\]
so the global order is one. For \(\dot u=\lambda u\), the amplification factor is
\[
G=1+\lambda\Delta t,
\]
and stability requires \(|G|\le1\).

## 2. Implicit Euler
\[
u^{n+1}=u^n+\Delta tF(u^{n+1},t^{n+1}).
\]
For \(F=-Au+b\),
\[
(I+\Delta t A)u^{n+1}=u^n+\Delta tb.
\]
It is first-order but robust for stiff diffusion.

## 3. Crank–Nicolson
\[
u^{n+1}=u^n+\frac{\Delta t}{2}
[F(u^n,t^n)+F(u^{n+1},t^{n+1})].
\]
The method is second-order for smooth solutions but may oscillate for strongly stiff/non-smooth problems.

## 4. BDF2
\[
\frac{3u^{n+1}-4u^n+u^{n-1}}{2\Delta t}=F(u^{n+1},t^{n+1}).
\]
Its global temporal order is two for constant step size. Variable-step coefficients must be derived from the actual step ratio rather than copied from the constant-step formula.

## 5. CFL and diffusion restrictions
For advection,
\[
Co=\frac{U\Delta t}{\Delta x}.
\]
For diffusion,
\[
Fo=\frac{\alpha\Delta t}{\Delta x^2}.
\]
The precise finite-volume definition must use the cell volume and face fluxes:
\[
Co_P\sim\frac{\Delta t}{V_P}\sum_f\frac{|F_f|}{\rho_f}.
\]
Stability limits are properties of the complete discretisation, not universal constants.

## 6. Pseudo-transient continuation
A steady residual \(R(u)=0\) can be approached through
\[
M\frac{u^{n+1}-u^n}{\Delta t}+R(u^{n+1})=0.
\]
Here \(\Delta t\) is numerical, not physical.

## 7. Restart
A restart must preserve the discrete state at \(t_n\), including all history variables required by multi-step schemes. BDF2 therefore requires enough previous states to reconstruct its derivative.

**CFDX path:** `src/cfdx/core/numerics/temporal.h`; application state/checkpoint: `src/cfdx/application/`, `python/cfdx/checkpoint.py`.

## 8. Verification
Use a manufactured ODE, diffusion eigenmode, temporal refinement and restart equivalence. For fixed spatial error,
\[
p_t=\frac{\ln(E_{\Delta t}/E_{\Delta t/r})}{\ln r}.
\]
Do not infer temporal order from a mixed space/time refinement campaign.