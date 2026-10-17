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
| Bound reporting | `test_n11_bound_reporting_hardening` | The covered bound-reporting checks pass |
| Laminar baseline only | `test_channel_qualification` | The existing analytical laminar channel regression passes; **not turbulence evidence** |

A missing registered test is reported as `BLOCKED_NOT_REGISTERED`, not skipped or passed.
A test failure is preserved as FAIL. The runner intentionally exits with code 2 even when
all listed component tests pass, until the missing turbulence qualification gates are
implemented and executed.

## Blocking work before a turbulent PASS

### Gate A — implementation audit

For each model, trace configuration -> initialization -> transport assembly -> linear solve
-> positivity/bounds -> eddy viscosity -> momentum flux assembly -> wall treatment.
Record exact source locations and tests. Confirm the model selected in the case is the
model actually used by the production solver. No silent model fallback is allowed.

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
