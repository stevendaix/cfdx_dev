# CFDX Verification & Validation

**Status: PARTIAL — governance is established; evidence migration is ongoing.**

The V&V Guide is the operational evidence manual for CFDX. Theory explains the mathematics; Developer explains implementation; V&V defines how claims are demonstrated and retained.

## Evidence vocabulary

- **Code verification:** implementation matches the specified mathematical algorithm.
- **Solution verification:** discretisation and iteration errors are controlled for the declared calculation.
- **Validation:** physical/reference comparison with stated uncertainty.
- **Qualification:** an explicitly bounded capability population satisfies declared acceptance criteria.

A green CI job is not automatically validation or qualification.

## Operational chapters

| # | Chapter | Purpose |
|---|---|---|
| 00 | V&V governance | Claims, evidence and status |
| 01 | Requirements and claims | Requirement-to-test traceability |
| 02 | Code verification | Implementation-level evidence |
| 03 | Solution verification | Iteration/discretisation error |
| 04 | MMS | Manufactured-solution campaigns |
| 05 | Spatial convergence | Mesh refinement and order |
| 06 | Temporal convergence | Time-step refinement |
| 07 | Conservation | Global/local balance and flux closure |
| 08 | Boundedness | Positivity, monotonicity and limiter evidence |
| 09 | Linear solver verification | Krylov/preconditioner evidence |
| 10 | Pressure–velocity verification | Coupling/projection evidence |
| 11 | Gradient verification | GG/LS/WLS/reconstruction evidence |
| 12 | Diffusion verification | Laplacian and non-orthogonal terms |
| 13 | Convection verification | Flux schemes and boundedness |
| 14 | Turbulence verification | Closure-level verification |
| 15 | Thermal verification | Energy/CHT evidence |
| 16 | Radiation verification | S2S/P1/DOM evidence |
| 17 | Multiphysics verification | Coupled residual/Jacobian evidence |
| 18 | Benchmark validation | Independent physical/reference cases |
| 19 | Qualification matrix | Population and acceptance gates |
| 20 | Evidence and reproducibility | Artifact retention and replay |

## Minimum evidence record

Every quantitative claim should identify, where applicable: software revision, case/setup revision, mesh and geometry, physics/model choices, numerical schemes, solver/stopping criteria, metric and oracle, acceptance criterion, execution environment and retained artifacts.

Never substitute a reference result for a CFDX result.

The docs/validation/ tree remains the current evidence store during migration.
