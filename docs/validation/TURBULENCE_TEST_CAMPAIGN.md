# CFDX turbulence qualification campaign — preparation

## Purpose and scope

This change prepares a reproducible entry point for the existing turbulence tests and
makes the remaining evidence gaps explicit. It does **not** claim that CFDX already
passes a turbulent physical benchmark. Component tests, transport-equation verification,
solver verification, and physical validation are distinct gates.

Authority: qualification governance issue #118; model inventory #461/N17;
method contract `docs/audits/09-turbulence-rans.md`; model matrix
`docs/validation/TURBULENCE_QUALIFICATION_MATRIX.md`.

## Run the prepared campaign

From a configured CMake build:

```bash
python3 scripts/run_turbulence_validation.py --build-dir build
```

The runner:
- discovers registered CTest names before executing anything;
- runs the equation-closure and hardening tests individually;
- runs `test_channel_qualification` only as a **laminar solver baseline**;
- writes a machine-readable JSON report (default:
  `artifacts/turbulence-qualification/report.json`);
- returns non-zero while required tests are missing, a test fails, or the physical/MMS
  qualification gates remain blocked.

Configure and build using the same project options as the normal developer/CI workflow.
For example, if the build directory is not yet configured:

```bash
cmake -S . -B build -DCFDX_BUILD_TESTS=ON
cmake --build build --parallel
python3 scripts/run_turbulence_validation.py --build-dir build
```

The configure/build commands are generic starting points; use the repository's pinned
compiler, dependencies, presets, and CI options when applicable. Do not treat a build
that silently disables HDF5/MPI-dependent targets as equivalent to the production build.

## Current executable gates

| Gate | Test(s) | What a PASS means |
|---|---|---|
| Closure/equation references | `test_turbulence_qualification_equations` | Tested closure identities agree with the test's frozen expressions |
| Phase 10 hardening | `test_phase10_turbulence_hardening` | The covered hardening regressions pass |
| SA negative-ν-tilde policy | `test_n11_sa_negative_nu_tilde_policy` | The explicitly covered negative-state policy passes |
| Bound reporting | `test_n11_bound_reporting_hardening` | The covered bound-reporting checks pass |\n| Production selection guard | `test_production_turbulence_contract` | The production driver accepts laminar declarations and rejects unsupported declared RANS models instead of silently running laminar equations |
| Laminar baseline only | `test_channel_qualification` | The existing analytical laminar channel regression passes; **not turbulence evidence** |

A missing registered test is reported as `BLOCKED_NOT_REGISTERED`, not skipped or passed.
A test failure is preserved as FAIL. The runner intentionally exits with code 2 even when
all listed component tests pass, until the missing turbulence qualification gates are
implemented and executed.

## Blocking work before a turbulent PASS

### Gate A — implementation audit (source audit completed; coupling remains blocked)

The source audit is recorded in
`docs/validation/TURBULENCE_PRODUCTION_PATH_AUDIT.md`. The case schema and Python setup catalogue
carry a turbulence model, and the C++ transport APIs implement component equations, but the
production driver does not yet call those transport entry points or couple `nu_t` back into
momentum. Therefore this is not yet a production RANS solver.

This PR adds a fail-closed contract at the production-driver boundary: a case declaring
non-laminar turbulence is rejected with an explicit message rather than being solved by the
laminar momentum equations. The test `test_production_turbulence_contract` covers laminar,
RANS, material-level and unknown model declarations. This guard is an interim correctness
barrier, not turbulence implementation.

The canonical C++ and Python catalogues now distinguish `COMPONENT_READY` from
`SOLVER_READY`. The six RANS entries are component-ready, excluded from the production
selector, and have transport MMS, production solver verification and coupled conservation
marked as incomplete in the JSON qualification matrix. Recommendations may still explain
SST/SA applicability, but warn that the current production solver cannot execute them.

The source audit also found and corrected repeated full-field gradient work inside the SST
per-cell loop. `cell_gradient_dot(k, omega, mesh)` is now evaluated once per outer iteration,
avoiding O(N²) work in that stage without changing its formula or convergence controls.

Next implementation work must connect one model end-to-end (recommended first target: SST):
initialization, wall distance, strain, F1/F2, transport, `nu_t`, momentum coupling and outer
convergence. Only then can a transport MMS test the actual production equation path.

### Gate B — isolated transport MMS

Add independent manufactured-solution tests for k-ε, k-ω SST, and SA. Each test must
define the exact transported field(s), analytic source, boundary data, dimensions/units,
discretization, and independent error norms. Use at least three systematically refined
meshes and report L1, volume-weighted L2, Linf, and observed order. Verify positivity and
finite source/flux terms without relying on clipping to hide an error.

The current Phase 13 acceptance document explicitly identifies this as missing:
`docs/validation/PHASE13_ACCEPTANCE.md`, section 13.9.

### Gate C — first physical RANS benchmark: TUR-CHANNEL

Use the frozen Moser–Kim–Mansour DNS reference and selected friction Reynolds number from
`docs/validation/TUR-CHANNEL_QUALIFICATION_CASE.md`. Before running, freeze:
- model and constants (recommend SST as the first wall-distance-dependent model);
- geometry, mesh family, wall-normal cell-centre distribution and target y+;
- inlet/outlet or periodic conditions and pressure-gradient/bulk-flow control;
- numerical schemes, relaxation, convergence and averaging criteria;
- quantitative acceptance tolerances before inspecting CFDX results.

Required outputs: U+, friction velocity, Cf, available Reynolds stresses, mass balance,
pressure-gradient/bulk-flow balance, all momentum and turbulence residual histories,
wall-distance/y+ distributions, iteration/runtime, and mesh sensitivity. A steady RANS
mean profile must be compared to the DNS mean; it must not be presented as resolving DNS
fluctuations.

### Gate D — subsequent physical benchmarks

1. TUR-FLATPLATE: compare Cf(x), boundary-layer thickness/profile quantities against the
   frozen reference and assess inlet turbulence and transition assumptions.
2. TUR-NACA0012: freeze Reynolds number, Mach number, incidence, transition assumption,
   wall treatment and mesh; compare Cp and integrated forces to the matching experimental
   or trusted numerical reference. Do not use a reference at incompatible Reynolds number
   or turbulence/transition regime.

### Gate E — evidence and promotion

Store the exact HEAD SHA, build configuration, commands, exit codes, logs, report JSON,
mesh/config hashes, solver controls, convergence histories, and reference-data version.
Promote a model/case only when its applicable gates pass. No tolerance inflation, fixed-
iteration PASS, oracle substitution, disabled tests, or substitution of a laminar case for
a turbulent one.

## Expected result of this preparation PR

This PR makes the current component-test subset reproducible and reports the missing gates.
It does not implement the production turbulent transport MMS or the physical cases, because
those require an audit of the existing solver APIs and a separately reviewed numerical
implementation. The PR must remain draft until its exact-head CI passes and the generated
report is reviewed.
