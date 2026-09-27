# Turbulence closure equation audit

This document is the equation-level closure contract for the complete 13-model CFDX turbulence catalogue in #461/#473/#474.

The tests in `tests/unit/test_turbulence_qualification_equations.cpp` evaluate the implemented closure against explicit reference expressions. A closure test is verification of the mathematical kernel; it is not physical validation.

## 1. Laminar

No turbulent closure:
\[
\nu_t = 0.
\]

## 2. Standard k-epsilon

Transported variables: \(k,\epsilon\).

Boussinesq closure:
\[
\nu_t=C_\mu\frac{k^2}{\epsilon}.
\]

Production used by the current incompressible transport path:
\[
P_k=\rho\nu_t S^2.
\]

The default CFDX constants are \(C_\mu=0.09\), \(C_1=1.44\), \(C_2=1.92\), \(\sigma_k=1\), \(\sigma_\epsilon=1.3\).

## 3. RNG k-epsilon

\[
\nu_t=C_\mu\frac{k^2}{\epsilon},
\qquad
\eta=S\frac{k}{\epsilon},
\]

with the RNG-modified production coefficient
\[
C_1^*=C_1-
\frac{\eta(1-\eta/\eta_0)}
     {1+\beta\eta^3}.
\]

The current source terms are
\[
P_k=\rho\nu_tS^2,
\quad
S_k=P_k,
\quad
S_\epsilon=C_1^*P_k\frac{\epsilon}{k},
\quad
D_\epsilon=C_2\rho\frac{\epsilon^2}{k}.
\]

## 4. Realizable k-epsilon

\[
\nu_t=C_\mu\frac{k^2}{\epsilon},
\]

where CFDX computes \(C_\mu\) from the strain, rotation and third invariant:
\[
C_\mu=\frac{1}{A_0+A_sU^*k/\epsilon}.
\]

The implemented epsilon coefficient is
\[
C_1=\max\left(0.43,\frac{\eta}{\eta+5}\right),
\qquad
\eta=S\frac{k}{\epsilon}.
\]

The epsilon destruction denominator is
\[
k+\sqrt{\nu\epsilon}.
\]

The implementation deliberately receives the tensor invariants rather than reconstructing them from a scalar strain-rate magnitude.

## 5. k-omega

Baseline Wilcox-type closure:
\[
\nu_t=\frac{k}{\omega}.
\]

The current transport path uses
\[
P_k=\rho\nu_tS^2,
\]
\[
S_k=P_k-\rho\beta^*\omega k,
\]
and the omega production/destruction terms
\[
S_\omega=
\alpha\frac{\omega}{k}P_k
+
\rho\sigma_d\frac{\nabla k\cdot\nabla\omega}{\omega},
\]
\[
D_\omega=\rho\beta_0\omega^2,
\]
with \(\sigma_d=\sigma_{d0}\) only when \(\nabla k\cdot\nabla\omega>0\).

The SST strain-rate limiter is not part of the baseline k-omega closure.

## 6. k-omega SST

The Menter SST eddy viscosity is
\[
\nu_t=
\frac{a_1k}
{\max(a_1\omega,SF_2)}.
\]

The blending is
\[
X=F_1X_1+(1-F_1)X_2
\]
for \(\sigma_k,\sigma_\omega,\beta,\gamma\).

The current transport implementation includes:
- k production with the explicit production limiter;
- \(-\beta^*\rho k\omega\) destruction;
- cross diffusion proportional to
  \[
  2(1-F_1)\rho\frac{\sigma_{\omega2}}{\omega}
  \nabla k\cdot\nabla\omega;
  \]
- sign-split cross-diffusion treatment in the omega equation;
- blended turbulent diffusion.

Wall distance is an explicit input to \(F_1\) and \(F_2\).

## 7. Spalart-Allmaras

Transported variable: \(\tilde\nu\).

\[
\chi=\frac{\tilde\nu}{\nu},
\qquad
f_{v1}=\frac{\chi^3}{\chi^3+c_{v1}^3},
\qquad
\nu_t=f_{v1}\tilde\nu.
\]

The modified vorticity is computed from \(\Omega\), \(f_{v2}\), \(\tilde\nu\), \(\kappa\) and wall distance.

Production:
\[
P= c_{b1}(1-f_{t2})\tilde S\tilde\nu.
\]

Destruction:
\[
D=
\left(
c_{w1}f_w-
\frac{c_{b1}}{\kappa^2}f_{t2}
\right)
\left(\frac{\tilde\nu}{d}\right)^2.
\]

The destruction coefficient uses the standard sixth-power formulation and the \(r\le10\) limiter. The implementation keeps the destruction coefficient signed; it is not artificially clamped to zero.

## 8. Smagorinsky LES

For filter width \(\Delta\):
\[
\nu_t=(C_s\Delta)^2S.
\]

CFDX currently uses the cell-volume cube-root as the default \(\Delta\) supplied to the closure kernel.

This is a closure kernel only; LES transport/temporal integration and physical qualification remain separate gates.

## 9. WALE LES

Using the tensor invariants
\(S^2\) and \(S_d^2\):
\[
\nu_t=
(C_w\Delta)^2
\frac{(S_d^2)^{3/2}}
{(S^2)^{5/2}+(S_d^2)^{5/4}}.
\]

CFDX requires the actual tensor invariants as inputs. It must not manufacture \(S_d^2\) from the scalar strain-rate magnitude.

## 10. Dynamic one-equation LES

The one-equation SGS closure has the form
\[
\nu_t=C_k\Delta\sqrt{k_{sgs}}.
\]

The coefficient \(C_k\) is a dynamic coefficient and therefore requires resolved/test-filter information. CFDX has the coefficient kernel and the algebraic closure kernel, but the complete dynamic transport/filtering path is still PLANNED.

A fixed \(C_k=0.1\) is therefore only a kernel default/test value, not evidence of a completed dynamic model.

## 11. DES

Current CFDX closure:
\[
\ell_{DES}=\min(d,C_{DES}\Delta),
\]
\[
\nu_t=(C_s\ell_{DES})^2S.
\]

The hybrid solver coupling and complete DES base-model path remain KERNEL_ONLY.

## 12. DDES

The shielding function is
\[
f_d=1-\tanh\left[(8r_d)^3\right],
\]
and
\[
\ell_{DDES}
=d-f_d\max(0,d-C_{DES}\Delta).
\]

The closure is
\[
\nu_t=(C_s\ell_{DDES})^2S.
\]

The actual computation and transport of \(r_d\) must be integrated before DDES can be promoted beyond KERNEL_ONLY.

## 13. IDDES

The current CFDX kernel first computes the DDES shielding and combines it with a stress/wake blending parameter \(f_b\):
\[
f_{IDDES}=(1-f_b)+f_bf_d.
\]

The implemented length scale is
\[
\ell_{IDDES}
=d-f_{IDDES}\max(0,d-C_{DES}\Delta),
\]
and
\[
\nu_t=(C_s\ell_{IDDES})^2S.
\]

This is an explicit algebraic kernel contract. Full IDDES requires the physically correct shielding, stress/wake functions, base-model coupling and transient 3-D verification before solver qualification.

## Qualification rule

For every model, the closure test must be independent of the implementation expression being tested. Passing the closure test establishes equation/kernel consistency only.

It does **not** establish:
- correct wall distance;
- correct discretisation of transported equations;
- convergence;
- mesh independence;
- correct y+ regime;
- physical validation;
- DES/DDES/IDDES hybrid-mode behaviour.

References used for cross-checking the model families include the OpenFOAM turbulence-model documentation and the original model literature. OpenFOAM's current documentation lists the corresponding RAS, LES and DES model families and their model-specific source implementations.
