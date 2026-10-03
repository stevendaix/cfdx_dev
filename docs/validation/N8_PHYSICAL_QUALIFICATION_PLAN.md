# N8 physical solver/preconditioner qualification plan

## Purpose

N8 implementation is now largely merged. The remaining work is **physical and production-path qualification**, not another broad linear-algebra implementation pass.

This plan turns the existing CFDX validation assets into a controlled N8 campaign without changing numerical acceptance criteria.

## Existing evidence to reuse

| Evidence | Existing executable | Role in N8 campaign |
|---|---|---|
| Couette | `tests/validation/test_phase9_acceptance.cpp` / `test_couette_quick` | Pressure-velocity / coupled production path |
| Poiseuille | `test_poiseuille_quick`, `test_poiseuille_diagnostics` | Analytic pressure/velocity reference and solver behaviour |
| Ghia Re=100 | `tests/validation/test_ghia_cavity.cpp` / `test_ghia_cavity_quick` | Non-trivial canonical CFD QoI validation |
| Algebraic Schur | `test_n8_schur_production_benchmark` | 4N Schur/AMG production-path evidence |
| Native AMG | `test_amg_preconditioner_qualification` | Multilevel AMG structural/energy evidence |
| Mesh robustness | `test_nonorthogonal_skew_campaign`, `test_convection_polyhedral_campaign` | Controlled difficult-mesh extension |
| Conservation | `test_conservation_boundedness`, `test_transport_conservation`, `test_conservation_assembly` | Independent conservation evidence |

The existing numerical workflow already classifies the N8 algebraic tests and the relevant physical validation tests separately. The campaign must preserve that distinction.

## Campaign matrix

### A. Baseline physical cases

Run, at minimum:

1. **Couette**
   - existing canonical setup;
   - compare the same physical case with the N8-selected linear/preconditioner path;
   - record nonlinear iterations, linear/FGMRES iterations, true residuals, continuity imbalance and velocity QoIs.

2. **Poiseuille**
   - use the existing analytic benchmark;
   - record pressure drop, flow rate, velocity-profile error and solver metrics;
   - repeat on the existing refinement family where available.

3. **Ghia Re=100**
   - use the existing reference implementation and current acceptance oracle;
   - record centreline velocity QoIs, nonlinear/linear convergence and conservation;
   - do not replace the existing Ghia oracle with solver residuals.

4. **Controlled skew/non-orthogonal case**
   - reuse the existing mesh-quality campaign infrastructure;
   - keep the geometry and mesh family fixed and reproducible;
   - measure solver robustness separately from discretisation error.

### B. N8 method matrix

For each case, record explicitly:

- linear solver;
- preconditioner;
- pressure strategy;
- momentum strategy;
- coupled-block strategy;
- Schur approximation (exact / SIMPLE / SIMPLEC / LSC / BFBt where applicable);
- AMG family and hierarchy parameters;
- null-space policy;
- linear tolerance and iteration limit;
- nonlinear convergence criteria;
- mesh/operator characteristics.

No automatic recommendation may silently substitute another method.

### C. Required quantitative outputs

For every run, retain machine-readable evidence for:

- matrix dimensions and NNZ;
- Schur dimensions/NNZ where applicable;
- nonlinear iteration count;
- Krylov/FGMRES iteration count;
- reported residual;
- independently recomputed true residual;
- continuity/global mass imbalance;
- relevant QoIs and reference errors;
- AMG hierarchy/coarse sizes;
- setup cost and solve cost as diagnostics;
- finite-value checks;
- failure reason when a run does not converge.

Timing and memory are diagnostic unless an explicit, reproducible performance gate is later defined.

## LSC/BFBt qualification

The existing #568 harness deliberately leaves approximation quality diagnostic. The physical campaign must therefore establish the acceptance envelope from representative CFD matrices.

For each representative matrix:

1. construct/use the exact Schur reference;
2. run LSC and BFBt with the documented scaling;
3. measure approximation error and conditioning;
4. measure resulting Krylov convergence and true residual;
5. compare against the attainable FP64 floor;
6. retain the complete matrix characteristics.

No single arbitrary absolute Schur-error threshold should be invented before this evidence exists.

## AMG difficult-matrix campaign

Add controlled operator families covering:

- anisotropy;
- strong coefficient scaling;
- increasing problem size;
- representative CFD saddle-point blocks where appropriate.

Separate:

- algebraic validity;
- hierarchy diagnostics;
- convergence;
- failure envelope;
- performance diagnostics.

A failed robustness point must remain visible rather than being hidden by tolerance changes or method substitution.

## Mesh / discretisation separation

The campaign must distinguish:

**solver failure** from **discretisation error**.

The N8 linear solver must converge tightly enough that the linear algebra error does not mask the spatial discretisation error of Couette, Poiseuille, Ghia or the skew/non-orthogonal benchmark.

Conversely, a successful Krylov convergence must not be interpreted as validation of the underlying discretisation.

## MPI boundary

Before claiming N8 production qualification beyond serial:

- either run the same representative production matrices/cases through the MPI path and compare independent true residuals, QoIs and reproducibility;
- or explicitly record MPI as outside the N8 qualification domain and link the follow-up work.

No MPI scaling claim is implied by the serial campaign.

## CI / acceptance sequence

The intended sequence for the implementation PR is:

1. Reuse existing tests and fixtures; do not duplicate existing physical solvers.
2. Add only the missing N8 production-path matrix/reporting layer.
3. Register the campaign in the N8 numerical workflow.
4. Run diagnostics first.
5. Apply hard gates only to independently justified quantities.
6. Run exact-HEAD CFDX CI, numerical maturity and numerical-maturity audit.
7. Run the long validation workflow separately when required.

No tolerance relaxation, disabled validation case, hidden fallback or performance-based method substitution is allowed.

## Definition of done for the next N8 qualification PR

The next PR should be considered complete when it provides:

- [ ] Couette N8 production-path matrix;
- [ ] Poiseuille N8 production-path matrix;
- [ ] Ghia N8 production-path matrix;
- [ ] controlled skew/non-orthogonal N8 case;
- [ ] explicit solver/preconditioner configuration in every run;
- [ ] independently recomputed true residuals;
- [ ] conservation evidence;
- [ ] LSC/BFBt versus exact-Schur evidence on representative CFD matrices;
- [ ] anisotropic/strong-scaling AMG diagnostics;
- [ ] machine-readable report;
- [ ] CI integration without suppressing existing gates.

This PR intentionally stops at the campaign design and evidence wiring. It does not claim N8 qualification before the executable production campaign has been run.

## Executable campaign layer

This PR now contains the first executable qualification layer in addition to the campaign specification.

- `scripts/n8_physical_qualification.py` reuses the existing CTest validation executables; it does not duplicate their physics or solver implementations.
- CMake registers `test_n8_physical_qualification` when Python is available.
- The campaign requires the existing Couette, Poiseuille, Ghia, controlled skew/non-orthogonal, Schur production benchmark, and AMG qualification tests to be present. A missing required test is reported as **INCOMPLETE**, not silently skipped.
- Tests are executed sequentially and the first failing gate stops the campaign while retaining its complete output in the JSON report.
- The report is written to `<build-dir>/n8_physical_qualification.json` and records pass/fail, elapsed time, captured output, and campaign policy metadata.
- No numerical tolerance is changed, no validation case is disabled, and no solver fallback is introduced by the campaign layer.

This is intentionally a first executable layer. The next iterations should extend the machine-readable evidence with explicit solver/preconditioner metadata and quantitative LSC/BFBt-vs-exact-Schur and difficult-matrix AMG evidence, using the existing tests as the numerical oracles.
