# Conservation and boundedness qualification

This document defines the first N11 qualification slice from #461. It adds reusable diagnostics without changing the discretisation or solver algorithm.

## Conservation contract

For a face flux phi_f, the owner receives +phi_f and an internal neighbour receives -phi_f. Internal faces therefore cannot create or destroy a conserved quantity algebraically.

The qualification reports:

- global boundary flux;
- L1/L2/max local cell imbalance;
- number of internal and boundary faces;
- non-finite face fluxes;
- explicit owner/neighbour accounting.

These diagnostics are deliberately independent from a linear-solver residual. A converged linear solve can still expose a physical mass imbalance, while an exactly conservative face assembly can be identified before nonlinear convergence.

## Boundedness contract

For a transported scalar with declared bounds [lower, upper], the diagnostic reports:

- minimum and maximum;
- count of lower/upper violations;
- largest lower/upper violation;
- non-finite values.

For positivity-constrained variables, audit_positive_field() is the zero-lower-bound specialization. Typical applications include k >= 0, epsilon > 0, omega > 0 and non-negative turbulent viscosity. Strictly positive variables should use a model-specific admissibility floor rather than silently clipping negative values.

The diagnostic tolerance is explicit and is not an acceptance-tolerance inflation mechanism: production qualification should normally use zero tolerance or a separately documented floating-point budget.

## Scope

This PR establishes the reusable invariant/diagnostic layer and executable regressions on closed and two-cell finite-volume topologies.

It does not yet claim full solver-level conservation of momentum/energy or boundedness of every production turbulence scalar. Those require attaching the diagnostics to the corresponding assembled face fluxes and post-solve fields in subsequent qualification PRs.

## Follow-up

- attach independent conservation reports to Couette, Poiseuille and cavity;
- add momentum component and energy/scalar balances;
- audit upwind, bounded SOU and TVD/MINMOD for boundedness and conservation;
- add non-orthogonal correction sensitivity;
- add per-cell/per-iteration diagnostics to solver reports;
- add model-specific positivity gates for turbulence and thermal fields.


## Implemented in the current qualification wave

The diagnostic layer is now also attached to the authoritative incompressible solver iteration:

- mass-flux boundary balance and global cell balance;
- local L1/L2/Linf continuity accounting and worst-cell identification;
- explicit non-finite authoritative face-flux detection;
- finite-velocity admissibility check after each nonlinear iteration;
- source and transient-accumulation balance API for transported quantities;
- scalar lower/upper-bound diagnostics with worst-cell identification;
- strict positive-floor diagnostics for turbulence-like variables.

The conservation identity used by the transport audit is

    sum(boundary flux) + sum(cell source) - sum(accumulation) = residual.

The accumulation term is optional for steady equations and must be supplied explicitly for transient equations. No source, flux or state is silently clipped.

## Discretisation qualification status

The solver already exposes UPWIND, bounded SECOND_ORDER_UPWIND and TVD/MINMOD. Existing transport tests cover their algebraic paths. The present conservation wave adds invariant diagnostics around those paths; it does not declare a scheme validated merely because its implementation exists.

Non-orthogonal diffusion/correction paths are already present in the FVM stack and are treated as a separate sensitivity/convergence gate. The acceptance criterion remains physical conservation and measured refinement, not a relaxed residual threshold.

## Remaining validation work

Production qualification still requires running the diagnostics on the complete benchmark matrix (Couette, Poiseuille and Ghia cavity) and, where applicable, energy/scalar and turbulence fields. The helper APIs are in place; benchmark-specific acceptance remains tied to measured results rather than assumed from unit tests.

## Assembly-level audit added

The qualification layer now distinguishes two different conservation contracts:

1. **Field-level face-flux conservation**: one authoritative face flux is applied as +phi_f to the owner and -phi_f to the neighbour.
2. **Assembly-level antisymmetry**: independently assembled owner and neighbour contributions are checked explicitly for owner_contribution + neighbour_contribution = 0 on internal faces, with zero neighbour contribution on boundary faces.

The second diagnostic is intentionally separate. A single stored face flux cannot reveal an upstream assembly bug that produced two inconsistent cell contributions before the flux was collapsed into one value.

reconstruct_cell_balance() also provides an independent post-assembly reconstruction of cell balances from the authoritative face flux. It is suitable for comparing the conservative flux path with independently reconstructed diagnostics without modifying the solution.

## Qualification gates

The current PR implementation is complete for the reusable diagnostic layer. The remaining gates are execution/evidence gates:

- compile and run test_conservation_boundedness;
- compile and run test_transport_conservation;
- compile and run test_conservation_assembly;
- execute Couette, Poiseuille and Ghia with the diagnostics enabled;
- record nonlinear/linear convergence, true residuals, global and local conservation metrics, worst-cell diagnostics and QoIs;
- compare UPWIND, bounded SECOND_ORDER_UPWIND and TVD/MINMOD without changing acceptance tolerances;
- execute the existing non-orthogonal/skew sensitivity campaign and retain quantitative conservation/refinement evidence;
- keep energy/scalar/turbulence positivity as model-specific follow-up gates.

A green unit-test build alone is therefore not a numerical-validation PASS.
