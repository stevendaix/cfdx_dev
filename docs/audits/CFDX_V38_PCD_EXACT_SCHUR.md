# CFDX V38: exact frozen pressure-Schur audit on a captured BFS state

**Scope.** V37 already established that all six Schur variants solve the same
original saddle-point matrix A,b,x0 in the 152-cell *simplified* BFS case.
V38 adds an opt-in, read-only capture of the pressure operators Mp/Kp/Fp
as built by the real PCD production branch. It independently computes the
Schur complement S=C-D(Auu)^(-1)G on the original saddle-point A, then
compares the two PCD approximation formulae with exact algebraic inverses.
No FV equation, wall condition, accepted nonlinear state, Krylov tolerance,
pressure gauge or default preconditioner is modified.

## Why this matters

The V37 PCD branch is historically a single negative Fp inverse. The actual
PCD lumped-mass commutator is a different, opt-in action:

- historical: P_legacy^(-1) = -Fp^(-1)
- commutator: P_PCD^(-1) = -diag(Mp)^(-1) Fp Kp^(-1)
- exact algebraic saddle Schur: S = C - D(Auu)^(-1)G

The script uses Moore-Penrose inverses to evaluate the ideal pressure actions
when Fp is singular. It does NOT substitute these inverses into CFDX.
The 'oracle-only exact Schur diagonal' is also diagnostic: it is NOT the
production SIMPLE, SIMPLEC or BlockLocal implementation.

## Reproduced evidence (152 cells / 608 coupled unknowns)

- A and b had the same SHA256 as V37:
  - A: `2d56f528ff69ec792d447a503d99ad36db06cec286cb95126d20921cc55b261c`
  - b: `6447a6034677cb224ecc90ee4f8f845b88bf4d0f83ec8e24cd59c4387fe72d61`
- Exact saddle S: numerical rank 152/152; retained condition number about 5.40e4.
- Original PCD Fp: numerical rank 151/152, whereas Kp is full rank 152/152.
- Ideal approximate inverse defect `||S P^(-1) - I||_F / sqrt(n)`:
  - historical negative Fp inverse: 864.9873653
  - lumped pressure commutator: 28.6036456
  - oracle-only diagonal of exact Schur: 0.9001040
  - exact Schur pseudoinverse reference: 2.89e-13
- The actual Schur RHS is also badly approximated: even its best one-scalar
  correction cannot lower the error below ~0.66 (historical) or ~0.85
  (commutator). Scale correction alone is therefore insufficient.
- The optical nullspace/inlet experiments on the same A did not close PCD
  FGMRES within 2000 iterations; other Schur variants passed the linear stage,
  but none qualified the global U,p,k,omega physical residual.

These numbers apply only to the first frozen state of the 152-cell debugging
mesh. They do NOT prove large BFS or NASA agreement.

## Production capture

Explicitly opt into `CFDX_PCD_OPERATORS_CAPTURE_PREFIX_V38=/some/output/pcd`.
The PCD branch writes up to THREE complete triplets per process:
`pcd_state1_mp.cfp`, `pcd_state1_kp.cfp`, `pcd_state1_fp.cfp`, etc.
Each file uses the existing, bounded CFP01 atomic binary writer with ZERO
placeholder RHS and initial vectors. Do not confuse these with genuine
`CFDX_PCD_FP_CAPTURE_PREFIX` action snapshots, which carry actual RHS.
An invalid capture is fatal rather than silently ignored.

## Replay

First capture A,b,x0 with V37 `--capture-original-operator`; the physical
initial state and selected numerical policy must match the PCD run. Then:

```sh
python scripts/diagnose_pcd_schur_v38.py \
  --coupled-A captures/pcd_A.mtx --coupled-b captures/pcd_b.mtx \
  --mp pcd_state1_mp.cfp --kp pcd_state1_kp.cfp --fp pcd_state1_fp.cfp \
  --output pcd_schur_v38.json \
  --max-action-error 0.5 --require-fp-nullity-match
```

V38 returns status 2 after writing a JSON audit when a **requested** strict
algebraic gate fails. That was the expected result on the 152-cell example;
this is not a failed program or relaxed tolerance. An ordinary audit without
these optional gates returns 0 after successfully collecting evidence.
The dense tool is intentionally limited to <=256 cells by default (up to
1024 only on explicit override); it scales cubically and is NOT a production
PCD/preconditioner for BFS 20,540, parallel CPU or GPU.

Optional CTest registration:

```sh
cmake -S . -B build-v38 -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_ENABLE_PCD_V38_DENSE_AUDIT_TESTS=ON
cmake --build build-v38 --target test_pcd_pressure_operator_capture_v38
ctest --test-dir build-v38 -R 'pcd_pressure_operator_capture_v38|diagnose_pcd_schur_v38'
```

The option is OFF by default. Enabling it requires NumPy and SciPy at
configuration time, fail-closed. The production capture feature is OFF by
default independently of the CTest option.

## Next qualification gates

1. Diagnose Schur-model mismatch before increasing Fp Krylov budgets again.
2. Compare PCD pressure BCs, gauge and viscosity scaling against the true
   Rhie-Chow Schur and the other five Schur actions, on identical states.
3. Construct an affordable, *clearly named hybrid* if necessary (e.g., PCD
   with a pressure-Schur coarse/defect correction). Never relabel LSC or
   SIMPLE fallback as pure PCD success.
4. Compare all six methods first on 712 cells, then V4 BFS 20,540 cells.
   Distinguish run success, true U,p,k,omega contraction and NASA validation.
5. Rerun the NACA, SST and thermal suite before changing defaults.

PETSc's separation between nonlinear operator and preconditioning matrix is
relevant: https://petsc.org/release/manual/snes/ ; FGMRES with variable
preconditioning is documented at https://petsc.org/main/manual/ksp/ .
