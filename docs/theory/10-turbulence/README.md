# 10 — Turbulence

**Status: REPOSITORY-GROUNDED; qualification remains a V&V question.**

Reynolds decomposition is

$$
u_i=\overline u_i+u_i'.
$$

Averaging introduces the Reynolds stress

$$
-\rho\overline{u_i'u_j'},
$$

which requires closure.

## Eddy-viscosity RANS

A representative closure is

$$
-\rho\overline{u_i'u_j'}
=
2\mu_tS_{ij}-\frac23\rho k\delta_{ij}.
$$

The definition of turbulent viscosity and transported variables is model-specific.

## Spalart–Allmaras

SA transports a modified viscosity variable with production, destruction, diffusion and wall-dependent terms. Current source family: physics/spalart_allmaras.h and turbulence transport/model files.

## k–ω SST

SST transports k and ω and blends near-wall and free-shear behaviour. Current source family: physics/sst_solver.h, turbulence_models.h and turbulence_transport.h.

## Wall distance

Wall distance enters near-wall closure and wall functions. Relevant files are wall_distance.h, wall_distance_turbulence.h and wall_functions.h. The wall-distance algorithm must be verified separately from turbulence-model validation.

## LES and hybrid methods

LES resolves part of the turbulent spectrum and models subgrid stresses. Hybrid methods combine resolved and modelled regions and are strongly mesh/time-resolution dependent.

## Numerical coupling

Turbulence equations use the common reconstruction, convection, diffusion, temporal and linear-algebra contracts. Source stiffness and positivity make their convergence requirements more restrictive in some regimes.

## V&V

Separate equation verification, transport MMS where meaningful, positivity/boundedness, canonical turbulent benchmarks and independent physical validation. Every benchmark records model constants, wall treatment, mesh resolution and inlet turbulence specification.
