# Turbulence production-path audit — PR #941

## Scope and evidence boundary

Audit performed against the PR base revision `120070dedebecd70e9f78f30691e9e5172a71925`.
This is a source audit, not a successful build or physical simulation. Exact-symbol searches,
the production executable entry point, model-selection metadata, solver APIs and existing tests
were inspected. CI must re-run on the final PR head.

## Finding 1 — case metadata exists, but production execution does not consume it

The case schemas carry a `turbulence_model` string:
- `src/cfdx/io/cfdx_io/case_schema.h`
- `src/cfdx/python/cfdx/io/schema.py`

The Python setup catalogue exposes RANS choices in `src/cfdx/python/cfdx/physics_setup.py`.
That establishes a configuration/UI contract, not runtime coupling.

The production executable `apps/production_solver/main.cpp` includes
`steady_incompressible_solver.h` and `production_energy_stage.h`. Its options and help text
expose mesh, restart, iteration/convergence, and velocity/pressure probes; no turbulence
selection or transported turbulence fields are exposed in the production run contract.
Exact-symbol searches for the transport entry points found their definitions and test callers,
but no production-driver caller:
- `solve_kepsilon_transport` in `src/cfdx/physics/turbulence_solver.h`;
- `solve_sst_transport` / `solve_sst_transport_dynamic_blending` in
  `src/cfdx/physics/sst_solver.h`;
- `solve_spalart_allmaras_transport` in `src/cfdx/physics/turbulence_solver.h`.

**Consequence:** the repository does not yet demonstrate that a case selecting SST, SA or
k-epsilon changes the production momentum solve. Until a fail-closed guard or complete coupling
is added, the model catalogue must not be interpreted as proof of a production RANS run. A
case-file field or a standalone transport function is not a coupled RANS capability.

## Finding 2 — transport functions are component-level solvers, not a demonstrated coupled RANS loop

The current APIs solve turbulence scalars on an already supplied mesh and face mass flux:
- k-epsilon and related families assemble scalar FVM equations in
  `turbulence_solver.h`;
- SST assembles k and omega equations in `sst_solver.h`, with F1/F2 supplied by the caller
  or computed by the separate dynamic-blending wrapper;
- SA transports `nu_tilde` with caller-supplied vorticity and wall distance.

The current tests (including `test_phase10_turbulence_hardening`,
`test_m2_m4_solver_validation`, and `test_turbulence_qualification_equations`) test
closures, source identities, guards, and small synthetic meshes. They are useful component
evidence, but do not show an outer iteration coupling turbulence, eddy viscosity, momentum
fluxes, wall treatment and convergence in `cfdx_production_solver`.

The generic k-epsilon production assembly uses hard-coded source linearisation in the solver
function. SST uses the `F1/F2` fields supplied to `solve_sst_transport`; the dynamic wrapper
is a separate entry point. The production caller must explicitly select and update the correct
path, not just serialize a model name.

## Finding 3 — quadratic repeated work in dynamic SST blending

In `solve_sst_transport_dynamic_blending`, the full vector returned by
`cell_gradient_dot(k, omega, mesh)` was computed inside the loop over every cell. Since that
helper returns one value per cell, the operation was repeated N times per outer iteration,
giving O(N²) work for this stage rather than one O(N) gradient pass plus the O(N) cell loop.

PR #941 moves the gradient contraction outside the per-cell loop. This is an algorithmically
equivalent hoist: it does not alter the formula, tolerances, relaxation, or stopping criteria.
CI/build evidence is required before accepting the change.

## Finding 4 — there is no isolated turbulence transport MMS yet

`docs/validation/PHASE13_ACCEPTANCE.md` item 13.9 explicitly lists turbulence MMS as missing.
The existing equation-reference test freezes algebraic closure/source identities, but does not
measure spatial L1, volume-weighted L2 and Linf errors or observed order for the production
transport operator on three or more refined meshes.

A valid MMS campaign cannot be claimed by wrapping a generic scalar-diffusion test around the
turbulence equations. The current transport entry points do not expose a manufactured forcing
term for each transported equation. Implementing a meaningful test requires a reviewed source
injection/test seam (or a manufactured residual assembled through the same production equation
path), independent analytic sources, exact boundary values, and mesh-refinement gates.

## Required implementation order

1. **Fail closed on runtime selection:** until coupled RANS exists, reject non-laminar
   `turbulence_model` in the production driver rather than silently running laminar flow.
2. **Couple one model end to end:** SST is the recommended first target because its eddy
   viscosity and wall-distance/F1/F2 dependencies make the data path explicit. Define where
   k/omega are initialized, how wall distance and strain are computed, how mass flux is
   updated, when nu_t is refreshed, and how outer convergence includes turbulence residuals.
3. **Add a production-path integration test:** case file selecting SST must either produce
   explicit coupled-turbulence evidence or fail closed; it must never silently execute laminar.
4. **Add transport MMS on the actual equation assembly path:** start with a scalar manufactured
   field and source-injection seam; then cover both transported SST equations and SA `nu_tilde`.
   Extend to k-epsilon only after the source/linearisation contract is independently checked.
5. **Run TUR-CHANNEL against frozen Moser–Kim–Mansour DNS** with mesh/y+ sensitivity and
   independent mass and pressure-gradient/bulk-flow balances.
6. **Run TUR-FLATPLATE and TUR-NACA0012** only after the coupled production path and channel
   benchmark pass their predeclared gates.

## Promotion decision

Current state remains **BLOCKED for physical RANS qualification**. No model/case status is
promoted by this audit. No tolerance is relaxed, no test is disabled, and no unit test or
laminar channel case is substituted for a turbulent benchmark.
