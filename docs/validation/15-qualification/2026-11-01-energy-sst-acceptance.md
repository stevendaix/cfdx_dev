# Independent energy closure and SST input contracts (2026-11-01)

Status: **code and unit tests implemented, physical campaigns not promoted**.

## Why

`solve_energy` previously accepted a small residual of the **assembled**
energy equation as sufficient evidence of convergence. Bounded convection
can subtract `T * div(m_dot)` in the assembled operator. If the carrier mass
flux does not conserve mass in each cell, `A*T=b` can be solved exactly while
the conservative heat/enthalpy transport balance is false. Opposing local
errors can also cancel in the global sum. Neither outcome may be reported
as an energy-converged state.

The independent reconstruction uses the actual accepted temperature and
produces, in watts,

`R_c = sum_faces(q_out,c) + rho*cp*V_c*(T_c - T_old,c)/dt - (Su_c + Sp_c*T_c)*V_c`.

The accepted-energy gate requires all of:

- Successful linear solve and the existing accepted-state and assembled
  backward-error tests.
- Finite face and cell conservation terms (`nonfinite_faces == 0`,
  `nonfinite_cells == 0`).
- `abs(sum_c R_c) / max(1 W, abs(net_boundary_flux), abs(source_power),
  abs(accumulation_power)) <= energy_balance_tolerance`.
- For **every** cell, `abs(R_c) / max(1 W, sum_faces(abs(q_f)) +
  abs(source_c*V_c) + abs(accumulation_c)) <= energy_balance_tolerance`.
  The implemented scale uses the **sum** of these gross contributions, not
  their maximum, so a cell with large throughflow is not incorrectly compared
  to an arbitrary domain-wide flux.

The `1 W` floor belongs to the existing independent audit contract; it is
not a new or relaxed tolerance. For fully relaxed *linear* energy stages,
`solve_energy` may ask Krylov to refine the algebraic solution by a decade
if the independent watts balance still fails. It does not loosen the
requested physical tolerance. If the accepted algebraic solution is exact
or stagnates, the stage terminates **not converged**. Under-relaxed
energy–radiation outer loops continue to converge normally.

`EnergyIteration::independent_conservation_worst_cell` retains its previous
absolute-residual meaning. New fields report the local normalized worst
cell separately, as well as the independent acceptance flag.

## SST hardening

The SST closure paths reject NaN/infinity in blending factors `F1/F2`,
strain, gradients and viscosity inputs **before clamping/max operations**.
They also reject nonrepresentable transport coefficients. This prevents
`std::max(finite, NaN)` from silently manufacturing plausible closure values.
This does not change the model constants, equation sources or near-wall
boundary condition scheme.

## Executable counterexamples

The two new CTest targets are:

- `test_energy_independent_acceptance`: valid balanced flux; single-cell
  nonconservative flux (algebraic PASS, physics FAIL); two cells with cancelling
  global errors (algebraic and global PASS, local physics FAIL); 16-cell
  `T=x^2` manufactured conduction, demanding independent closure.
- `test_sst_closure_input_contract`: healthy SST, nonfinite F1/F2,
  out-of-range F1/F2, negative/nonfinite strain, nonfinite cross-gradient and
  eddy-viscosity inputs.

On the pre-fix source, these counterexample suites failed, as desired.
After the fix both targets pass in Release. Relevant thermal and SST
validation executables also pass in a local `-O0` regression build.

Run from the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target test_energy_independent_acceptance test_sst_closure_input_contract
ctest --test-dir build --output-on-failure -R '^test_(energy_independent_acceptance|sst_closure_input_contract)$'
```

## Limits, integration and next work

- A converged independent energy balance is **necessary but insufficient**
  for physical validation. Accuracy vs reference/MMS, three-mesh GCI,
  conservation of carrier mass, wall heat flux, pressure-velocity coupling
  and performance still need independent campaigns.
- Existing branches for COUPLED/AMG, NACA, SST restart and CHT must remain
  separate. Do not bypass their gates or promote incomplete cases.
- This patch targets the **reported** PR #1155 head `9ced2008` as its
  intended integration base, but was created from a local reproduction of
  that branch. The remote head and CI are **not verified in this environment**.
  Recheck the actual HEAD, patch preimages and CI before creating a PR.
