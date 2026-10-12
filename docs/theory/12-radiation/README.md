# 12 — Radiation

## 1. Blackbody and gray emission
\[
E_b=\sigma T^4,
\qquad
E=\varepsilon\sigma T^4.
\]
The \(T^4\) dependence makes radiation strongly nonlinear in temperature.

## 2. Radiosity
For an opaque diffuse-gray surface,
\[
J=\varepsilon\sigma T^4+(1-\varepsilon)G,
\]
where \(G\) is irradiation. Net heat flux into the surface can be written
\[
q=G-J
\]
with the sign convention explicitly declared by the implementation.

**CFDX:** `src/cfdx/physics/radiation_s2s.h`.

## 3. View factors
\[
A_iF_{ij}=A_jF_{ji},
\qquad
\sum_jF_{ij}=1
\]
for a closed enclosure. Reciprocity and closure are geometry-only invariants and are ideal unit/verification tests.

## 4. S2S system
For each surface,
\[
J_i=\varepsilon_i\sigma T_i^4+(1-\varepsilon_i)\sum_jF_{ij}J_j.
\]
Thus
\[
\left[I-(I-\mathcal E)F\right]J=\mathcal E\,E_b,
\]
schematically, where \(\mathcal E\) is the diagonal emissivity matrix.

## 5. P1
A representative gray P1 equation is
\[
-\nabla\cdot\left(\frac{1}{3\beta}\nabla G\right)+\beta G
=4\beta\sigma T^4.
\]
The exact CFDX convention for \(G\), \(\beta\), scattering and boundary closure must be followed.

## 6. DOM
For discrete direction \(s_m\),
\[
s_m\cdot\nabla I_m+\beta I_m=S_m.
\]
Angular quadrature approximates
\[
G=\int_{4\pi}I\,d\Omega
\approx\sum_mw_mI_m.
\]

**CFDX:** `src/cfdx/physics/radiation_streaming.h`, `radiation_models.h`, `radiation_solver.h`.

## 7. Thermal coupling
\[
\rho c_p\frac{DT}{Dt}=\nabla\cdot(k\nabla T)+q_{rad}+q_v.
\]
Radiation-to-energy sign must be tested by a two-surface exchange case.

## 8. Verification
Check Stefan–Boltzmann limits, view-factor reciprocity/closure, energy conservation and deterministic S2S/P1/DOM regression. See `scripts/thermal_radiation_regression.py`.