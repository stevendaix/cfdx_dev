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

test_thermophysical_models checks:
- analytic scalar-property models;
- tabulated interpolation and all extrapolation policies;
- Boussinesq and ideal-gas density;
- constant-cp enthalpy;
- representative turbulence closures;
- turbulent thermal diffusivity;
- rejection of invalid physical inputs.

This test is a model-layer verification gate, not a replacement for solver-level V&V.
