# M1-M4 Verification and Validation

CFDX separates verification (mathematical/numerical implementation correctness) from validation (agreement with independent reference data). The implementation stack now contains executable solver paths for M1-M4; this document defines the acceptance gates that remain attached to the final PR.

## M1 — Incompressible

Implemented:
- cell-centred finite-volume momentum assembly;
- continuity and pressure-correction system;
- SIMPLE nonlinear iteration;
- PISO/PIMPLE multi-correction control;
- Rhie-Chow pressure-gradient face-flux correction;
- pressure gauge/reference handling;
- momentum/continuity diagnostics.

Regression coverage:
- Poiseuille analytical oracle;
- constant-state and zero-state fixed-point tests;
- FVM conservation/diagonal checks.

Benchmark target:
- lid-driven cavity, Re=1000, using the Ghia et al. reference family;
- quantitative centerline profiles, mass imbalance and mesh refinement.

The classical Re=1000 lid-driven cavity is a standard benchmark; published benchmark descriptions use a unit square cavity with a moving lid and compare centerline velocity profiles against Ghia et al. (1982). 

## M2 — Turbulence

Implemented:
- k-epsilon transport equations;
- SST k-omega transport equations with F1/F2 blending;
- local turbulent diffusivity with harmonic face interpolation;
- positivity clipping;
- wall epsilon/omega closures;
- existing Smagorinsky/DES eddy-viscosity foundations.

Regression coverage:
- positivity and boundedness;
- production/destruction paths;
- wall-law closure;
- k-epsilon and SST transport solver smoke/verification cases.

Benchmark targets:
- fully developed turbulent channel;
- zero-pressure-gradient flat plate;
- NASA TMR comparisons for the pinned turbulence-model variants.

## M3 — Thermal / CHT

Implemented:
- steady energy equation;
- transient finite-volume energy equation;
- fixed, zero-gradient and spatially varying boundary temperatures;
- cell-dependent conductivity/diffusivity;
- two-region CHT interface matching;
- two-layer thermal-resistance interface temperature;
- interface heat-flux balance.

Regression coverage:
- transient one-cell analytical energy step;
- two-region interface temperature/flux balance.

Benchmark targets:
- 1-D two-material slab;
- conjugate channel;
- multi-region energy conservation.

## M4 — Radiation

Implemented:
- blackbody and gray-surface models;
- two-surface exchange;
- view-factor closure/reciprocity checks;
- discrete ordinates quadrature validation;
- participating-medium directional intensity transport;
- absorption/scattering source terms;
- irradiation and radiative energy source;
- nonlinear radiation/energy outer coupling.

Regression coverage:
- six-direction normalized DOM quadrature;
- isothermal blackbody participating-medium equilibrium;
- radiation source closure.

Benchmark targets:
- gray two-surface enclosure;
- canonical participating-medium enclosure;
- radiation/energy coupled balance.

## Acceptance gates

The final implementation is accepted only after:
1. Release compilation succeeds with warnings enabled.
2. DebugSanitizers compilation and CTest succeed.
3. Analytical and manufactured-solution errors are reported as L2/Linf norms.
4. Conservation errors are reported independently of linear residuals.
5. Mesh-refinement order is measured for the verification cases.
6. Restart/reproducibility is checked.
7. Serial/parallel equivalence is checked where MPI execution is enabled.

The implementation therefore does not use a residual alone as a validation criterion. The solver architecture follows the standard pressure-predictor / pressure-correction / flux-update sequence used by SIMPLE/PISO/PIMPLE implementations, while the thermal/radiation stack follows the same conservative finite-volume assembly principles used by established CFD frameworks.