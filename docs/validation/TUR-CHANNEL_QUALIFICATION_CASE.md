# TUR-CHANNEL — turbulent channel qualification

## Reference definition

Moser, R. D., Kim, J. & Mansour, N. N. (1999), “Direct numerical simulation of turbulent channel flow up to Re_tau=590”, *Physics of Fluids*, 11(4), 943–945, DOI: 10.1063/1.869966.

Primary data portal: https://turbulence.oden.utexas.edu/MKM_1999.html

The DNS populations are Re_tau = 180, 395 and 590, with profiles and statistical databases including mean velocities, Reynolds stresses, higher moments, spectra, correlations and budgets. The 180 case has stronger low-Reynolds-number effects; do not impose a universal log-law fit across all three cases.

## Canonical setup

- Plane channel with half-height h and walls at y=±h.
- Periodic streamwise and spanwise directions.
- No-slip walls.
- Fully developed flow driven by a controlled mean pressure gradient or body force. The selected forcing convention must match the chosen reference data.
- Re_tau = u_tau*h/nu, where u_tau = sqrt(abs(tau_w)/rho).

## Required quantities

- U+ versus y+.
- u_rms+, v_rms+, w_rms+.
- -u'v'+ and k+.
- Viscous, turbulent and total shear-stress profiles.
- Cf, bulk Reynolds number and pressure-gradient/flow-rate balance.

A linear eddy-viscosity RANS closure can be compared on mean velocity and shear stress, but its Boussinesq hypothesis does not independently reproduce the full DNS Reynolds-stress anisotropy. Use an RSM if the six Reynolds-stress components are a primary validation target. RANS cannot be compared against DNS spectra as if it resolved temporal turbulence.

## Data integrity

Before the first campaign run, freeze each downloaded reference file independently:

- exact source URL and filename;
- retrieval date and source revision/history;
- SHA-256 checksum;
- units, coordinates and normalization;
- interpolation method used by the comparison tool.

Never reconstruct reference profiles from plotted figures.

## Required CFDX evidence

An actual production turbulent CFDX run, statistically stationary averaging interval, wall-normal and homogeneous-direction mesh sensitivity, y+ and resolution metadata, independent mass/momentum balance, bulk-flow stability, quantitative profile errors, runtime and exact software/case revision.

**Status: BLOCKED.** The reference definition is known, but this is not evidence that the current production turbulence path executes or is validated. Reference-data parsing tests alone do not qualify CFDX.
