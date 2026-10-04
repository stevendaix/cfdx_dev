# 14 — Numerical Analysis

## 1. Truncation error
For continuous operator \(L\) and discrete operator \(L_h\),
\[
\tau_h=L_hu-Lu.
\]
Consistency requires \(\tau_h\to0\) as \(h\to0\) under the regularity assumptions of the scheme.

## 2. Convergence and observed order
If
\[
E(h)=Ch^p,
\]
then
\[
p_{obs}=\frac{\ln(E_h/E_{h/r})}{\ln r}.
\]
This formula is meaningful only for a controlled refinement family with a stable error norm and unchanged mathematical problem.

## 3. Stability
For
\[
\dot u=Lu,
\]
a one-step method gives
\[
u^{n+1}=G(\Delta tL)u^n.
\]
Stability requires the eigenvalues of the amplification operator to lie in the method's stability region.

## 4. Monotonicity/boundedness
A bounded scalar solution should satisfy
\[
\phi_{min}\le\phi_h\le\phi_{max}.
\]
A coefficient condition such as
\[
a_P\ge\sum_Na_{PN}\ge0
\]
is a useful diffusion-dominated diagnostic, but not a universal theorem for every CFD scheme.

## 5. Dissipation/dispersion
For Fourier mode \(e^{ikx}\), a numerical method introduces amplitude and phase errors. Dissipation alters amplitude; dispersion alters phase. These must be studied independently for transport/wave problems.

## 6. Conditioning
\[
\kappa(A)=\|A\|\|A^{-1}\|.
\]
Large \(\kappa\) amplifies perturbations and can explain solver sensitivity without implying a physical modelling error.

## 7. Error budget
Conceptually,
\[
E_{total}\sim E_h+E_t+E_i+E_m+E_g,
\]
for spatial, temporal, iterative, modelling and geometry/input contributions. This is a diagnostic decomposition, not an exact statistical identity.

## 8. Mesh effects
Observed order can be corrupted by changing topology, aspect ratio, skewness or boundary treatment while refining. Each V&V campaign must therefore state its refinement family.

**CFDX paths:** `src/cfdx/core/numerics/`, `src/cfdx/core/linalg/`; reference experiments are under `docs/theory/04-gradients-reconstruction/experiments/`.