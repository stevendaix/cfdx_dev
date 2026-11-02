# Thermohydraulic turbulent benchmark contracts

Each benchmark has an independent GitHub Actions workflow and a machine-readable reference contract.

| Case | Reference | Independent CI | Main comparison targets |
|---|---|---|---|
| TH-01 | Moser–Kim–Mansour DNS, nominal Re_tau 180/395/590 | `.github/workflows/th01-channel-dns.yml` | U+, turbulence intensities, Reynolds stresses, total shear-stress balance, Cf, Re_b |
| TH-03 | Driver–Seegmiller straight-wall BFS, AGARD CMP32 | `.github/workflows/th03-backward-facing-step.yml` | x_r/h = 6.26 ± 0.10, Cp, Cf, mean velocity and Reynolds stresses |
| TH-07 | Wagner–Wetzel mixed-convection DNS | `.github/workflows/th07-mixed-convection.yml` | Re_b=4328, Pr=0.71, Gr sweep, derived Ri, wall shear, U/T and turbulence statistics |
| TH-08 | OECD/NEA–Vattenfall T-junction | `.github/workflows/th08-t-junction.yml` | inlet flow split, velocity profiles, T*, T*_rms, 3.5 Hz spectral feature |

## Meaning of CI status

The four workflows are deliberately **reference-contract CIs**: they validate that the frozen reference conditions, units, uncertainties, normalizations and derived dimensionless values are internally consistent. They do not run CFDX and do not claim a physical-validation PASS.

The physical status remains BLOCKED until the production turbulence model is coupled to momentum, energy and (where required) buoyancy, then run against the reference on mesh/time-step refinement levels. A successful reference-contract CI must never be promoted to a solver qualification.

## Comparison and uncertainty policy

- Compare like-for-like quantities, stations, coordinate systems and normalizations.
- Experimental uncertainties are stored with the case contract and must be included in the comparison metrics.
- DNS statistical sampling uncertainty is distinct from discretization error.
- RANS mean-field comparisons are separated from resolved RMS, spectra and Reynolds-stress anisotropy.
- Mesh and time-step uncertainty must be assessed independently from model-form error.
- Do not digitize figures when machine-readable reference profiles exist.
- No pass tolerance may be relaxed after seeing CFDX results.

## Primary references

- Moser et al. (1999): https://turbulence.oden.utexas.edu/MKM_1999.html
- Driver & Seegmiller CMP32: https://torroja.dmt.upm.es/turbdata/agard/docs/AGARD-AR-345.pdf
- Wagner & Wetzel (2022): https://doi.org/10.1002/gamm.202200006
- OECD/NEA Vattenfall (2011): https://www.oecd-nea.org/upload/docs/application/pdf/2021-03/csni-r2011-5.pdf
