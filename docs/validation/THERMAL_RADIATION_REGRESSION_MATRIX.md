# Thermal / Radiation Regression Matrix

This matrix is the permanent regression gate for CFDX thermal and radiation physics. It separates component verification from solver-level V&V and requires an independent oracle or conservation invariant for every executable case.

## Thermal coverage

| Family | Models / case | Oracle or invariant | Gate |
|---|---|---|---|
| Scalar property | CONSTANT | exact reference value | executable |
| Scalar property | LINEAR | analytic slope | executable |
| Scalar property | POLYNOMIAL | analytic polynomial | executable |
| Scalar property | POWER_LAW | analytic power law | executable |
| Scalar property | EXPONENTIAL | analytic exponential | executable |
| Scalar property | ARRHENIUS | analytic Arrhenius relation | executable |
| Scalar property | TABLE | interpolation + CLAMP/LINEAR/REJECT | executable |
| Bounds | property limits | bounded output / invalid configuration rejection | executable |
| Density | CONSTANT | exact value | executable |
| Density | LINEAR_TEMPERATURE | analytic law | executable |
| Density | BOUSSINESQ | analytic law + validity-limit rejection | executable |
| Density | IDEAL_GAS | p/(RT) | executable |
| Density | TABLE | interpolation | executable |
| Enthalpy | constant Cp | analytic integral | executable |
| Enthalpy | linear Cp | analytic integral | executable |
| Enthalpy | polynomial Cp | analytic integral | executable |
| Enthalpy | tabulated Cp | piecewise-linear integral | executable |
| Energy solver | 1-D conduction | exact cell-centre temperatures + energy balance | executable |
| Energy solver | volumetric generation | independent conductive balance | executable |

## Radiation coverage

| Family | Model / case | Oracle or invariant | Gate |
|---|---|---|---|
| Blackbody | Stefan-Boltzmann emission | exact sigma T^4 | executable |
| Gray wall | net surface flux | epsilon(Eb-G) | executable |
| P1 | diffusion coefficient | analytic coefficient | executable |
| P1 | absorption/reaction coefficient | analytic coefficient | executable |
| P1 | equilibrium source | zero source at equilibrium | executable |
| P1 solver | isothermal cavity | exact equilibrium irradiation and zero source | executable |
| DOM | isotropic quadrature | zeroth/first/second moments | executable |
| Rosseland | radiative conductivity | analytic 16 sigma T^3/(3 kappa_R) | executable |
| Selector | thin/intermediate/thick regimes | threshold invariants | executable |
| S2S | view-factor bounds/closure | physical bounds | executable |
| S2S | reciprocity | Ai Fij = Aj Fji | executable |
| S2S | gray radiosity | two-surface Stefan-Boltzmann resistance | executable |
| S2S | energy conservation | integrated net power | executable |
| S2S | external irradiation | equilibrium limit | executable |
| S2S Monte Carlo | ray visibility/blocking | statistical visibility + blocker invariant | executable |
| Input guards | invalid radiation states | deterministic exception | executable |

## CI regression policy

The CI campaign runs all thermal/radiation suites and generates:

- individual stdout logs;
- a machine-readable JSON report;
- a Markdown summary;
- residual and energy-balance diagnostics.

A regression is a CI failure when:

1. any executable returns non-zero;
2. any test case fails;
3. required residual/balance diagnostics disappear;
4. a model loses its analytical or conservation gate.

The suite does not use a silently updated numerical baseline. Analytical references are compiled into the tests so a code change cannot make the reference pass by merely changing the stored expected result.

## Current V&V boundary

The regression suite is deliberately broader than the current solver-level V&V, but it does not claim coverage that is not implemented. The next independent solver-level extensions are:

- thermal manufactured-solution verification;
- transient thermal MMS and time-order verification;
- spatial refinement/order for production thermal schemes;
- CHT interface flux conservation and coupled convergence;
- radiation angular convergence (S2/S4/S6/S8);
- radiation/energy coupled conservation;
- spectral multi-band verification where the spectral model is enabled.

These are tracked separately from component-level regression so a passing algebraic test is never reported as a full CFD validation case.
