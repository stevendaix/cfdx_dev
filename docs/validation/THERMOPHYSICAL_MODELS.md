# Thermophysical and turbulence model framework

## Scope

This extension separates thermophysical properties from the energy solver and provides explicit model selection, validation and extrapolation policy.

### Temperature-dependent scalar properties

Supported forms:
- constant;
- linear;
- polynomial;
- power law;
- exponential;
- Arrhenius;
- tabulated piecewise-linear data.

The same scalar model infrastructure can represent thermal conductivity, heat capacity, dynamic viscosity and turbulent Prandtl number. Every evaluation rejects non-finite temperatures and non-positive physical results.

### Density

Supported forms:
- constant;
- linear temperature law;
- Boussinesq;
- ideal gas;
- tabulated data.

The Boussinesq law is intentionally exposed as an algebraic property model. Buoyancy source assembly remains a solver concern.

### Enthalpy

The ThermophysicalProperties enthalpy helper integrates cp(T) independently of the energy equation. This makes the thermodynamic property model verifiable before transient enthalpy transport is coupled.

### Tabulated data

Temperature abscissae must be strictly increasing. Extrapolation is explicit:
- CLAMP: use the nearest endpoint;
- LINEAR: continue the endpoint segment;
- REJECT: throw.

### Turbulence catalogue

The model catalogue covers:
- laminar;
- standard, RNG and realizable k-epsilon families;
- k-omega and SST;
- Spalart-Allmaras;
- Smagorinsky and WALE LES;
- DES, DDES and IDDES hybrid families.

The current PR deliberately distinguishes model closures and selectors from complete transport-equation implementations. A closure helper must not be represented as a fully validated RANS/LES solver. Full equation implementations require their own production, destruction, diffusion, wall-treatment and benchmark validation.

### Turbulent heat transfer

The framework exposes the standard coupling:

alpha_t = nu_t / (rho Pr_t)

with explicit validation of density and turbulent Prandtl number.

## Verification

test_thermophysical_models and test_thermophysical_models_vv check:
- analytic scalar-property models;
- tabulated interpolation and all extrapolation policies;
- Boussinesq and ideal-gas density;
- constant-cp enthalpy;
- representative turbulence closures;
- turbulent thermal diffusivity;
- rejection of invalid physical inputs;
- analytic linear-property gradient verification;
- exact piecewise-linear table interpolation/integration;
- Boussinesq validity-boundary verification;
- DDES/IDDES shielding limits and protected/unprotected length-scale branches.

The V&V executable remains a model-layer verification campaign; wall-resolved heated-channel T+/y+ correlation requires a complete turbulence transport and wall-treatment implementation and is therefore not claimed by this PR.

This test is a model-layer verification gate, not a replacement for solver-level V&V.


## Defensive validity and hybrid shielding

- All property evaluators reject non-finite temperatures and non-positive absolute temperature.
- Boussinesq density has an explicit configurable validity guard on \(|\beta(T-T_0)|\), defaulting to 0.1. Disable the guard only when the caller intentionally accepts the approximation outside this range.
- Spalart–Allmaras accepts a negative working variable by clamping it to zero by default; callers can select explicit rejection with `NegativeNuTildePolicy::REJECT`.
- DDES exposes the standard shielding-function shape \(f_d=1-\tanh((8r_d)^3)\) and a shielded length-scale helper. The helper is an algebraic building block; it is not a substitute for the complete DDES transport equations, wall treatment, or IDDES zonal formulation.
- The turbulent Prandtl number remains a property-model abstraction. An overload accepts local \(y^+\) so a future wall-dependent Prandtl model can be introduced without changing the energy-solver interface.

These checks are deliberately explicit: a closure helper must not silently convert invalid input into a plausible-looking CFD result.
