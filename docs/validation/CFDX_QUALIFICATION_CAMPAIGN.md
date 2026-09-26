# CFDX Engineering Qualification Campaign

Issue #118 is the authoritative qualification programme for CFDX. It is deliberately different from a collection of unit tests or capability checks.

The campaign answers one engineering question:

> Can CFDX execute a reproducible CFD case, conserve the governing quantities, converge the numerical solution, reproduce an independent reference quantity, demonstrate mesh behaviour, and keep that evidence under regression control?

## 1. Qualification matrix

The machine-readable source of truth is `CFDX_QUALIFICATION_REGISTRY.json`. The following matrix is the human-facing programme.

| ID | Domain | Case | Current state | Primary QoI |
|---|---|---|---|---|
| LAM-COUETTE | Laminaire | Couette | PASS | velocity profile, mean velocity |
| LAM-POISEUILLE | Laminaire | Poiseuille | PASS | velocity profile, flow rate |
| INC-GHIA | Incompressible | Ghia cavity | DIAGNOSTIC | centreline velocity, vortices |
| VER-MMS | Conservation | MMS | READY | L2/Linf field error, conservation |
| FORCE-VMFL036 | Force | Sphere Re=100 | READY / #441 | pressure/viscous/total drag, Cd |
| EXT-NACA0012 | External flow | NACA0012 | PLANNED | CL, CD, Cm, Cp |
| INT-CHANNEL | Internal flow | Channel | READY | profile, pressure gradient, friction |
| SEP-BFS | Separation | Backward-facing step | PLANNED | reattachment length, Cp, wall shear |
| TUR-FLATPLATE | Turbulence | Flat plate | PLANNED | Cf, thicknesses, profiles |
| TUR-CHANNEL | Turbulence | Channel | PLANNED | mean velocity, Cf, Reynolds stresses |
| TUR-NACA0012 | Turbulence | NACA0012 | PLANNED | CL, CD, Cm, Cp |
| TH-CONV | Thermal | Conduction/convection | READY | T, heat rate, Nu |
| CHT-INTERFACE | CHT | Solid/fluid interface | READY | interface T and heat flux |
| RAD-P1 | Radiation | P1 | READY | wall heat flux, radiative source |
| RAD-DOM | Radiation | DOM | READY | wall heat flux, radiative source |
| RAD-S2S | Radiation | S2S | READY | view factors, net heat rate |

**Important:** READY means that a case can be implemented/executed; it is not a numerical PASS. DIAGNOSTIC means that an executable exists but the evidence is not yet sufficient for qualification.

## 2. Required evidence for every case

Every promoted case must have a calculation sheet containing, without exceptions:

1. **Reference**
   - publication, analytical solution, DNS/experimental database, or independent code;
   - provenance and version/date where applicable.

2. **Parameters**
   - all dimensional inputs;
   - dimensionless numbers;
   - material properties;
   - operating conditions.

3. **Mesh**
   - topology and geometry;
   - cell/node count;
   - minimum quality metrics;
   - wall resolution/y+ where applicable;
   - coarse/medium/fine definitions;
   - refinement ratio.

4. **Boundary conditions**
   - every named boundary;
   - values and units;
   - pressure reference/gauge treatment;
   - turbulence/thermal/radiation BCs when applicable.

5. **Numerical schemes**
   - gradient reconstruction;
   - convection;
   - diffusion;
   - pressure interpolation;
   - temporal scheme;
   - turbulence/radiation angular discretization where applicable.

6. **Solver settings**
   - algorithm (SIMPLE, SIMPLEC, PISO, PIMPLE, COUPLED, etc.);
   - linear solver/preconditioner;
   - linear tolerances and iteration limits;
   - relaxation/coupling parameters;
   - initialization;
   - stopping policy.

7. **Convergence criteria**
   - nonlinear residuals;
   - normalized residuals;
   - QoI stabilization;
   - statistical stationarity for turbulence;
   - no fixed-iteration PASS.

8. **Conservation criteria**
   - mass imbalance;
   - momentum balance;
   - energy balance;
   - interface balance;
   - radiation closure as applicable.

9. **QoI**
   - precisely defined extraction formula;
   - location/surface/volume;
   - units.

10. **Reference value**
    - scalar, profile, curve, or field;
    - independent provenance.

11. **Error**
    - absolute error;
    - relative error;
    - norm definition;
    - interpolation/sampling procedure.

12. **Observed order**
    - at least three meshes where a formal order is meaningful;
    - error sequence and refinement ratio;
    - observed order;
    - GCI when appropriate.

13. **Runtime**
    - mesh generation time;
    - solver wall time;
    - iteration count;
    - hardware/CI runner;
    - parallelism.

14. **Regression status**
    - executable test;
    - deterministic inputs;
    - CI result;
    - retained logs/artifacts;
    - explicit promotion status.

## 3. Qualification gates

A case is **PASS** only when all gates below are satisfied:

### Gate A — Reference integrity
The reference configuration and provenance are frozen. A convenient value copied from another case is not accepted.

### Gate B — Real CFDX execution
The quantity of interest is produced by CFDX from the numerical solution. An analytical oracle alone is never a CFDX PASS.

### Gate C — Convergence
The solver has converged according to residual and physical-QoI criteria. A finite number after a fixed iteration count is not convergence evidence.

### Gate D — Conservation
The relevant global balances close within the case-specific acceptance contract.

### Gate E — Quantitative agreement
The CFDX QoI is compared against the independent reference using a declared error metric.

### Gate F — Grid verification
Where applicable, at least three systematically refined meshes are used to distinguish discretization error from model/reference error.

### Gate G — Regression
The complete case can be rerun in CI and produces machine-readable diagnostics and retained evidence.

### Gate H — Audit
The calculation sheet records the exact configuration, result, error, observed order, runtime and regression status.

## 4. Verification versus validation

The campaign explicitly separates:

- **Code verification:** MMS, analytical solutions, conservation and numerical-order tests.
- **Solution verification:** iterative convergence, mesh convergence, conservation and sensitivity studies.
- **Validation:** comparison with independent experimental/DNS/reference data.

A case can therefore be numerically verified without being physically validated, and a validation comparison cannot compensate for poor numerical verification.

## 5. Current execution strategy

The first qualification wave should reuse the cases that already have solver-level infrastructure:

- Couette;
- Poiseuille;
- Ghia;
- MMS;
- VMFL036 sphere;
- existing thermal/CHT/radiation executables.

The second wave adds canonical external/internal/separation cases.

The third wave promotes the turbulence matrix only after the turbulence models themselves have passed their model-verification and numerical checks.

This ordering avoids claiming that a turbulence validation case is meaningful while the underlying model, wall treatment, discretization or convergence behaviour is still under qualification.

## 6. No tolerance inflation

The following are explicitly prohibited:

- widening a tolerance to make a case green;
- disabling a failing case;
- replacing a CFDX result by an analytical value;
- hard-coding a reference value as the computed QoI;
- accepting a plausible QoI when conservation is wrong;
- silently changing the mesh/domain/reference configuration;
- interpreting capability presence as validation.

A failure is evidence. It must remain visible until its root cause is understood and corrected.

## 7. References

- Ghia, Ghia & Shin, *High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method*, JCP 48 (1982), 387–411.
- Roache, *Verification of Codes and Calculations*, AIAA Journal 36 (1998), 696–702.
- Armaly, Durst, Pereira & Schönung, *Experimental and theoretical investigation of backward-facing step flow*, JFM 127 (1983), 473–496.
- Ladson (1988), NACA 0012 experimental aerodynamic data; the exact dataset/subcases must be frozen in the NACA0012 calculation sheet before promotion.
- VMFL036 literature configuration as documented in the CFDX validation record and PR #441.

## 8. Relationship with issue #118

#118 remains the authoritative campaign tracker.

This registry is the structured qualification contract. Individual implementation PRs may add or repair solver capabilities, but a case is promoted only from actual CFDX evidence under this programme.

The existing Fluent/VMFL matrix remains useful for compatibility planning, but it must not be confused with the engineering qualification matrix defined here.
