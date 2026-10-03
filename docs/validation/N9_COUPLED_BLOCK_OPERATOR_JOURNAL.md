# N9 — COUPLED block-operator investigation journal

This is a working record of an unfinished investigation, not a deliverable and not user documentation. It documents why the block `COUPLED` pressure–velocity path in `src/cfdx/physics/steady_incompressible_solver.h` does not reach the same discrete fixed point as the segregated paths, the derivation that identifies the operator defect, the defects found and fixed, the measurements taken, and the ordered next steps.

Two sessions are recorded. Sections 1 to 6 are the first session: the operator analysis and a failed attempt. Section 7 is the second session: the fixes that were applied to the solver header, the measurements that confirmed or refuted each hypothesis, and the failure that is still open. The test file is back to its committed form; `git show HEAD:tests/validation/test_n9_physical_matrix.cpp` is the reference. Branch: `feat/n9-pressure-velocity-closure`. Parent campaign document: `docs/validation/N9_PRESSURE_VELOCITY_CLOSURE.md` (sections 1 and 3 there define the block system and the Rhie–Chow operator this journal works with; nothing here contradicts it, it adds the executable detail).

## 1. Objective and gate

The requirement is that the block `COUPLED` algorithm — momentum and continuity solved together as one `4N x 4N` system — reaches the **same discrete fixed point** as `SIMPLE`, `SIMPLEC`, `PISO`, `PIMPLE` and `FRACTIONAL_STEP`.

The gate is physical cross-algorithm equivalence on the three internal benchmarks in `tests/validation/test_n9_physical_matrix.cpp` (channel Couette, pressure-driven Poiseuille, lid-driven cavity), plus the controlled skew Couette mesh. It is explicitly *not* "each algorithm reports converged". The reason is structural: an incomplete block continuity operator is still a consistent discrete system, so GMRES converges on it, the outer loop's own convergence test can be satisfied, and both runs print "converged" while their velocity fields separate. The cross-algorithm max-difference gates (`tests/validation/test_n9_physical_matrix.cpp:424`) are the only measurement that can distinguish the two situations.

No gate currently passes on `COUPLED`. See section 6 for the measured state.

## 2. Defects found and fixed in the working tree

Each item is stated as: the defect, why it makes the coupled operator differ from the segregated one at the discrete-operator level, and the fix as it now exists in the diff.

1. **Arithmetic average for internal Gauss face pressure in the momentum pressure-coupling block.** `solve_coupled_momentum_continuity` used `coeff = 0.5` on both cells. The segregated path builds the pressure gradient through `gauss_gradient_with_boundary` (`steady_incompressible_solver.h:220`), which interpolates the face value with `w = dn/(dc+dn)` where `dc` is the owner-centre-to-face distance and `dn` the neighbour-centre-to-face distance (`:253-261`), and the gradient is `sum_f Sf.(w p_c + (1-w) p_other)/V` (`:274-279`). An arithmetic average coincides with that operator only when the face is equidistant from both centres. On a skewed mesh it is a different linear form in `p`, i.e. a different momentum operator, so the block fixed point cannot equal the segregated fixed point even when both converge. Fix: `w_own = dn/(dc+dn)`, `w_other = 1 - w_own` with a degenerate-distance throw, at `steady_incompressible_solver.h:833-848`.

2. **`EMPTY` boundary patches in the momentum pressure-coupling block.** `gauss_gradient_with_boundary` skips `PatchType::EMPTY` faces entirely (`:264-266`), so they contribute nothing to the pressure gradient. The block was adding an owner extrapolation on those faces, inventing a gradient component the segregated operator does not have. On any case with empty walls (the internal benchmarks have two empty patches each, `test_n9_physical_matrix.cpp:75-76`) the two operators differ by construction. Fix: `continue` on `PatchType::EMPTY` before any pressure boundary lookup, at `steady_incompressible_solver.h:853-861`; the same skip is mirrored in the `gw` assembly (`:954-956`) and in the continuity boundary branch (`:1110-1112`).

3. **Sign of the deferred non-orthogonal remainder in the continuity row.** The reconstructed flux is

   ```
   phi_f = rho*interp(U)_f . Sf_f - rho*rfn*[(p_N-p_P)/d * A_ortho + gp_f . Snon]
   ```

   (`:386-388` for the pressure-gradient part, `:515` for the assembly). Summing `Sf_c . phi_f` over a cell's faces gives the orthogonal part as `-D*(p_P - p_N)` on the left-hand side of the row (which is how the block assembles it, `:1074-1075`), and the non-orthogonal remainder `+rho*rfn*gp_f.Snon` on the right. So `div(phi)=0` requires the remainder to be **subtracted** from the row. Fix: `b(row) -= nonorth_flux` at `steady_incompressible_solver.h:1107`; it was `+=`.

4. **`SIMPLEC` modified-inverse (`rAtU`) variant removed.** The consistent inverse `rAtU = 1/(A_P - sum|off-diag|)` produced an inconsistent pair: `build_hbya` used `rAU` while the pressure corrector and the Rhie–Chow flux used `rAtU`, so the pressure gradient entered the velocity twice with two different weights. That inverse is also unusable as an operator — `1/(A_P - sum|off-diag|)` diverges on convection-dominated diagonals and reconstructs a divergence-free field that violates momentum. Fix: `build_hbya` now takes only `rAU` (signature at `:1563-1567`), all `rAtU` computation and use is deleted (the `rAU`/`rAtU` selection in `rhie_chow_pressure_flux_internal`, `directional_face_coefficient`, the boundary pressure diagonal and the `U` reconstruction now unconditionally use `rAU`; see `:1917-1919`, `:1960`, `:1976`, `:2003-2006`, `:2110-2116`, `:2124`), and `SIMPLEC` instead re-solves the **already assembled** momentum matrix with the increment

   ```
   rhs_d[c] = (grad p_new - grad p_old)_d[c] * V_c        ->   A (u_new - u_pred) = -V (grad p_new - grad p_old)
   ```

   then subtracts the increment from `U`, and rebuilds `HbyA` and the Rhie–Chow flux from it (`:2127-2193`). Reusing the existing matrix matters: `A_P`, and therefore the pressure operator built from it, is exactly what the predictor left it at, which re-assembling would not preserve. Because that re-solve moves the flux again, `SIMPLEC` is given a second pressure corrector (`:1674-1679`).

5. **`PISO` inner-loop `phiHbyA` must be the convective predictor flux.** The inner loop previously carried `phiHbyA = mass_flux`, the absolute conservative flux. `mass_flux` already contains `rho*rfn*grad p` from `make_rhie_chow_mass_flux` (`:473-474`, `:515`); feeding it back as the predictor adds that operator a second time on top of the pressure equation, so the corrected flux can never reach `div(mass_flux)=0` and the continuity residual sticks at a constant floor. Fix: `phiHbyA = make_mass_flux(mesh, geometry, HbyA, controls.density, velocity_bcs)` at `:2195-2207`, with the momentum matrix still frozen inside the inner loop.

6. **`COUPLED`: segregated corrector skipped, flux rebuilt from `HbyA`.** The block solve *is* the pressure solve, so running the segregated corrector on top of it re-solves the same correction with a different operator and the two iterations limit-cycle against each other. `pressure_correctors = 0` for `COUPLED` at `:1668-1673`. Since the corrector loop is then empty, the authoritative flux has to be rebuilt after it: `mass_flux = make_rhie_chow_mass_flux(mesh, geometry, HbyA, p, rAU, ...)` at `:2210-2217` (and again immediately after the block solve at `:1839-1845`, which also fixed the volume factor: the `rAU` arrays handed to the flux are now `1/A_P`, not `V/A_P`, since `make_rhie_chow_mass_flux` introduces `V` itself, `:360-366`). Two supporting changes belong to this item: `ux/uy/uz` are seeded from the block solution before `build_hbya` (`:1893-1902`), otherwise `HbyA` is built from zeros on the coupled path; and the block increment is accepted without under-relaxation (`:1829-1833`), because the block solve is already a converged linear solve and discarding part of it pins the reported divergence at the discarded fraction.

7. **Per-iteration re-baselining of `result.reference_momentum_residual` removed.** Previously the reference was assigned *after* the iteration's solve, under `if (iter == 1)`, and immediately consumed on the next line. Because the retry controller can re-enter iteration 1, the reference was effectively re-baselined on every attempt at iteration 1, and it was taken from the **corrected** state. For `COUPLED` that is fatal: iteration 1 is already a converged linear solve, so numerator and denominator both sit at the round-off floor and `momentum_equation_residual_relative` can never fall below approximately 1, whatever the accuracy of the state. Fix: the reference is now computed at the top of iteration 1, **before** the solve, from the pre-solve state (`:1776-1794`), and the assignment after the solve is gone. Consequence for the metric: `h.momentum_equation_residual_relative` is now `final_momentum_residual / reference` with `reference` taken from the state the solve started from, so it is a genuine reduction ratio rather than a ratio against itself. For the internal matrix the initial state is the rest state, so the reference is O(1) and the ratio can reach the 1e-10 range. Subtlety not yet addressed: on a retry of iteration 1 the reference is recomputed from whatever state the rollback restored.

8. **Diagnostics: `BLKRES`.** A block under `if (diagnostics.coupled_matrix_summary)` prints `BLKRES continuity=... momentum=...`, the max of `|b - A x|` over the continuity rows and over the momentum rows respectively, computed from the post-gauge finalized matrix and the initial guess `x = (U_old, p_old)`, immediately before `solve_gmres` (`:1394-1421`). Because it is evaluated at the initial guess, it is an absolute measure of how far the assembled system is from being satisfied by the state it was assembled from — its value is as a consistency check on assembly, not as a convergence metric (the post-solve check is `COUPLED_LINEAR_TRUE_RESIDUAL`, `:1433-1476`).

## 3. The central derivation: the block continuity row is incomplete

This is the substantive defect, and it is the reason the `COUPLED` path cannot be equivalent to the segregated paths even after items 1 to 8.

Write the segregated velocity reconstruction. `build_hbya` (`:1563-1590`) forms

```
H[c] = rAU[c] * ( b[c] + V_c * grad(p)_d[c] - sum_{j != c} A_PN(c,j) * u_j[c] )
```

and the pressure corrector reconstructs the corrected cell velocity as

```
u = H - rAU * V * grad(p_new)          (:2110-2116)
```

Now substitute the momentum equation that `b` came from, `b[c] = A_P(c,c) u[c] + sum_{j != c} A_PN(c,j) u_j[c]`, into the expression for `H`:

```
H[c] = rAU[c] * ( V_c * grad(p)_d[c] + A_P(c,c) * u[c] )
     = u[c] + rAU[c] * V_c * grad(p)_d[c]
```

That is an exact algebraic identity, not an approximation. Hence

```
HbyA = U + rAU * V * grad(p)                                        (*)

flux(HbyA)_f = rho * interp(HbyA)_f . Sf_f
             = rho * interp(U)_f . Sf_f  +  rho * interp(rAU*V*grad(p))_f . Sf_f
```

Define the second term

```
T_f = rho * 0.5 * ( G_d[P] + G_d[N] ) . Sf_f ,
G_d[c] := rAU_d[c] * V_c * grad(p)_d[c]
        = rAU_d[c] * sum_g Sf_g . d * ( w_g * p_c + (1-w_g) * p_nb(g) )
```

using the same `0.5` face interpolation that `make_rhie_chow_mass_flux` uses for its cell field (`:472`) and the same `Sf_g` face stencil that `gauss_gradient_with_boundary` uses (`:235-275`).

The segregated pressure equation measures continuity on the flux of `HbyA` (`:1935-1936`, `:1950-1966`). So the segregated continuity operator, written on the unknowns `(U, p)`, is

```
0 = div( rho*interp(U).Sf )  +  T  -  Q_orth(p)  -  Q_non(p)
```

with `Q_orth(p) = sum_f Sf_c . rho*rfn*(p_N-p_P)/d*A_ortho` and `Q_non(p) = sum_f Sf_c . rho*rfn*gp_f.Snon`. The term `T` is the momentum-response term: it is exactly the contribution the velocity unknowns make when they are eliminated to form the pressure equation.

The block continuity row (`:978-1174`) currently assembles only

```
0 = div( rho*interp(U).Sf )  -  D_faces * (p_P - p_N)  +  b_rest  -  Q_non(p)
```

with `D_faces = rho * rfn * area / d` implicit (`:1048`, `:1074-1075`) and the boundary contributions in `b_rest`. **The term `T` is absent.** The six convective entries per internal face (`:1004-1010`) carry `rho*0.5*Sf.d` on the two adjacent velocity unknowns and nothing else.

The same statement in block-elimination terms. Eliminating `u` from

```
[ M   G ] [ u ]   [ b_u ]
[ D   C ] [ p ] = [ b_p ]
```

gives the Schur complement `(C - D M^-1 G) p = b_p - D M^-1 b_u`, and the `-D M^-1 G` term *is* `T` with `M^-1` approximated by `diag(M)^-1 = rAU` — the same approximation the segregated pressure operator makes when it builds `D` from `rAU` (`:1976-1977`, `:1610-1628`). In the block solve that contribution is generated implicitly by the momentum rows, so it does not have to be written down. What *does* have to be written down is the continuity row that defines it: the `D` that the momentum response is contracted with. The row carries the orthogonal face pressure derivative `D_faces`, but it does not carry the matching momentum-response coupling, so the velocity unknowns in the block see a pressure gradient only through their own momentum rows and the continuity row sees only `interp(U)`.

The practical consequence is an algebraic identity, not a heuristic. Let `R(c) = div(mass_flux)` be what the outer loop actually measures (`:2392-2405`, with `mass_flux` built from `HbyA` by `make_rhie_chow_mass_flux`). Using (*):

```
R(c) = div( rho*interp(U).Sf )  +  T  -  Q_orth(p)  -  Q_non(p)
```

Substituting the block row's own fixed point, `div( rho*interp(U).Sf ) = D_faces*(p_P-p_N) + Q_non(p) - b_rest`, and noting that `D_faces*(p_P-p_N)` is the block's approximation to `Q_orth(p)`:

```
R(c) = T  +  ( D_faces*(p_P-p_N) - Q_orth(p) )  -  b_rest
```

The bracket is the residual mismatch between the block's orthogonal pressure derivative and the flux's, which is small on an orthogonal mesh. The leading term is `T`. So the block row can be satisfied to machine precision while the authoritative flux divergence equals `T` and not zero. That is a *fixed point of a different operator*, and it is exactly what the observed stall looks like: the `COUPLED` cavity reaches a state where `velocity_change_inf` and `pressure_change_inf` are both identically zero from iteration to iteration while `continuity_linf` stays at a fixed non-zero value (measured values in section 6). The identification of that floor with `T` is the hypothesis the `HOMOCOUNT` probe exists to confirm or refute; the algebra above is verified from the code, the attribution of the specific number is not.

Omitting `T` still converges, and that is why this was not caught earlier. The incomplete block system is nonsingular and consistent; it just solves a different discrete problem whose fixed point is offset from the segregated one by the momentum response. Both runs report convergence.

### The wiring that is missing

`G_d[c]` is a linear form in `p` (the Gauss gradient is linear in `p`, and `rAU` and the geometry are frozen in the Picard linearisation), which is why the precomputation described in section 5 is possible at all. The contribution is `rho*0.5*Sf_f.d*(G_d[c] + G_d[o])` on **both** continuity rows of the face, i.e.

```
for each internal face f with cells (P,N), outward Sf_f for P:
    for d in 0..2:
        for (col, coef) in gw[d][P]:
            A.push_back(nv+P, col, rho*0.5*Sf_f.d * coef)
            A.push_back(nv+N, col, rho*0.5*Sf_f.d * coef)
        b(nv+P) -= rho*0.5*Sf_f.d * gw_const[d][P]
        b(nv+N) -= rho*0.5*Sf_f.d * gw_const[d][N]
```

with `col` being a global pressure column `nv + cell` and the constant terms moved to the right-hand side because the row is written `A x = b`. The sign is positive on the left-hand side: `T` belongs on the same side as the convective `interp(U)` term, matching `make_rhie_chow_mass_flux`, which forms `phi = rho*interp(HbyA).Sf - correction`.

Boundary faces must use the same `w`-free branches the `gw` assembly already encodes (`EMPTY` skipped, `FIXED_VALUE` into `gw_const`, everything else owner extrapolation), which is what `gauss_gradient_with_boundary` does at `:262-273`.

## 4. Failed first attempt, and the measurement that killed it

The first attempt added the momentum-response term using only the **current** face's area vector, roughly

```
base = 0.5 * rho * (dAU_o + dAU_n) * Sf_f.d^2
```

i.e. the same shape as the convective entries but with `Sf_f.d^2` instead of the cell's face-stencil gradient. Measured, starting `COUPLED` from `SIMPLE`'s converged cavity state and running 2 block solves:

```
HOMOCOUNT dev_after_2_block_solves=4.5143e+01
```

The block did not merely drift, it left the segregated fixed point entirely: a deviation of order 45 in a velocity field whose lid speed is 1. A two-iteration perturbation experiment is the right instrument here — if the coupled operator agrees with the segregated one, starting from the segregated fixed point must leave it at round-off.

Root cause, and it is a category error in the discretisation rather than an indexing or sign slip. `G_d[c]` is

```
G_d[c] := rAU_d[c]*V_c*grad_p_d[c]
        = rAU_d[c] * sum over ALL faces g of cell c of  Sf_g.d * ( w_g * p_c + (1-w_g) * p_nb(g) )
```

— a sum over the whole face stencil of the cell, divided by nothing, because the `V_c` and the `1/V_c` of the gradient cancel. It is emphatically **not** one face. Summing only the current face misrepresents the gradient by construction: the gradient of a cell is a sum over its entire boundary, and replacing it by a single face makes the "gradient" a function of whichever face the assembly happens to be visiting. That the error is large and not a small perturbation is the expected outcome, not a surprise.

The failed code has been removed; the current tree contains the precomputation instead (section 5). The `4.5143e+01` figure therefore cannot be reproduced from the current tree. For reference, the same probe on the **current** tree — which has the precomputation but not the wiring — reads

```
HOMOCOUNT dev_after_2_block_solves=1.0461e-01
```

## 5. Current state of the tree

Everything in this section is verified against the working tree unless explicitly marked otherwise.

### 5.1 The `gw` / `gw_const` precomputation exists and is not yet consumed

`solve_coupled_momentum_continuity` contains, at `steady_incompressible_solver.h:911-976`:

```
std::array<std::vector<std::map<std::size_t,double>>,3> gw;      // linear forms in p, keyed by global pressure column nv + cell
std::array<std::vector<double>,3>                      gw_const; // constant part, per cell and component
```

It is assembled per cell and per component over **all** of that cell's faces, with the same Green-Gauss weight `w = dn/(dc+dn)`, the same owner/neighbour resolution rule, the same `EMPTY` skip and the same boundary treatment as `gauss_gradient_with_boundary`: `EMPTY` faces are skipped (`:954-956`), `FIXED_VALUE` pressure goes into `gw_const` (`:968-971`), any other pressure boundary becomes owner extrapolation onto the `nv + c` column (`:972`). Topology validation throws on invalid neighbour/owner and degenerate-distance cases (`:930-947`), matching the momentum block and the continuity block.

The comment block above it (`:893-910`) states the derivation this journal develops in section 3.

**These forms are not consumed.** A search over the file shows `gw[d][c]` and `gw_const[d][c]` written at `:949`, `:950`, `:970`, `:972` and read nowhere; no continuity-row entry is added from them. The wiring described at the end of section 3 is the immediate next step. Until it lands, the block continuity row is exactly the operator derived in section 3 with `T` missing.

### 5.2 The `coupled_matrix_summary` diagnostic does fire on the current tree

This corrects an earlier working note that the `BLKRES` line did not appear. `coupled_matrix_summary` is declared in `DiagnosticsControls` (`:69-77`), passed through as `controls.diagnostics` to `solve_coupled_momentum_continuity` (`:1807`), and the two `HOMOCOUNT` blocks in the test set it (`:317`, `:388`, `:390`). The build in `build-verify` is current (`cmake --build build-verify --target test_n9_physical_matrix` reports no work to do) and running `./build-verify/test_n9_physical_matrix` does emit the line, twice, from the two block solves of the pre-loop `HOMOCOUNT` probe:

```
COUPLED_MATRIX_SUMMARY n=1024 nnz=17245 M_l2=13.4393 G_l2=1.45371 D_l2=1.36788 C_l2=1.07687 zero_rows=0 zero_cols=0 schur_diag_min=0.0206203 schur_diag_absmax=1 schur_near_zero=0 rhie_chow_schur_vs_algebraic_diag_delta=0.985679 gauge_row=768
COUPLED_PRECONDITIONER name=CoupledBlockSchur-Diagonal restart=512 adaptive_restart=0\nBLKRES continuity=0.0198316 momentum=0.000774408
COUPLED_LINEAR_TRUE_RESIDUAL norm=6.31225e-10 relative=7.89031e-10 gmres_reported=6.31225e-10
...
BLKRES continuity=8.14426e-05 momentum=0.000721001
```

The earlier "did not appear" observation is therefore stale — most plausibly a stale binary or a run of the committed test rather than the instrumented one — and no cause is asserted here. The flag's declaration and propagation are verified as correct as described.

Two cosmetic defects in that diagnostic block, worth cleaning up while it is being touched: `wd` at `:1410` is assigned and never used, and the `rau3` inverse-diagonal array built at `:1395-1401` is never used. Also note `" adaptive_restart=0\\n"` at `:1392` is a literal backslash-n, not a newline, which is why `BLKRES` lands on the `COUPLED_PRECONDITIONER` line and is harder to grep than it should be.

Also worth interpreting before trusting it: `rhie_chow_schur_vs_algebraic_diag_delta=0.985679` is dominated by the gauge row (the diagnostic compares against the pre-overwrite vector, so the reference cell sees algebraic `1.0` against a raw accumulated `~0.014`), not by a per-cell operator defect. **Unverified:** that this fully explains the value.

### 5.3 Cost of the map-based form storage

`gw` is `3 * nc` `std::map`s keyed by pressure column. For each cell and component the entries are the cell itself plus one entry per internal face of that cell (owner column and neighbour column), so the storage is O(nc * degree) map entries per component and the assembly does O(degree) `operator[]` insertions per face per component. Boundary faces contribute at most one entry or one constant each.

For the 256-cell cavity this is irrelevant. For production mesh sizes it is a genuine concern: `std::map` node allocation per entry, per component, per outer iteration, on the pressure-dependent part of the coupled solve. This is stated as a concern, not as a measured cost — no timing has been taken, and no production mesh size has been run through this path. A flat CSR-like structure (offset + column array + value array per component, assembled with a stencil walk) would remove the allocation and the pointer chasing and is the obvious replacement if it turns out to matter.

### 5.4 `coupled_pressure_schur_diagonal` is diagonal-only

The solver accumulates `coupled_pressure_schur_diagonal[c] += D` per internal face and per fixed-pressure boundary face (`:1076`, `:1169`) and passes it to the preconditioner (`:1358`). That accumulation is a diagonal-only approximation of the momentum-weighted face-flux derivative, and it is a remaining approximation, not a resolved item.

Two qualifications, verified in `src/cfdx/core/linalg/preconditioner.h`: the preconditioner builds its own Schur approximation as `S_tilde = C - D diag(M)^-1 G` with the **full** pressure-pressure sparsity pattern (`:328-375`), not from the supplied diagonal; and `set_pressure_schur_diagonal` is consulted only when a Schur row is degenerate, to supply a fallback diagonal (`:399-405`) — and the reference cell's entry is overwritten with `1.0` at `steady_incompressible_solver.h:1346` before the call, so in practice the solver-supplied vector contributes essentially nothing. The vector's real remaining role is the `rhie_chow_schur_vs_algebraic_diag_delta` diagnostic. **Unverified:** whether the diagnostic is worth keeping in this form, and whether a non-diagonal, momentum-weighted Schur quantity is needed anywhere else.

### 5.5 Open items found while reading, not yet triaged

These are code observations, not measured results, and none of them is part of the current diff. They are recorded because they are exactly the kind of thing the skew campaign exposes.

- **Direction used for `rfn`.** `rhie_chow_pressure_flux_internal` builds `rfn` from the centre-to-centre unit vector `e` (`:355`, `:366`), while `directional_face_coefficient` and the block continuity row build it from the face normal `Sf/|Sf|` (`:1617`, `:1626`, `:1016-1018`, `:1043`). These coincide on an orthogonal face and differ on a non-orthogonal one. Both the segregated pressure matrix and the block use the face normal, so the two paths agree with each other; neither reproduces the conservative flux it is supposed to represent.
- **Orthogonal area.** The block and the segregated pressure matrix both use `area/d` for the implicit orthogonal coefficient (`:1048`, `:1978`), while the flux uses `(p_N-p_P)/d * A_ortho` with `A_ortho = Sf . e` (`:367`, `:386-388`). Same situation: cross-algorithm consistent, not flux-consistent on non-orthogonal faces.
- **Volume factor on the boundary pressure diagonal.** `directional_face_coefficient` includes `V` (`:1618-1626`) and the block's fixed-pressure boundary branch includes `V` (`:1155-1163`), but the segregated pressure matrix's fixed-pressure boundary diagonal does not (`:2003-2007`). This is a units-level inconsistency between two operators that the three current benchmarks do not exercise, because `channel_pressure_bc()` is `ZERO_GRADIENT` on every patch (`test_n9_physical_matrix.cpp:176-182`). It would show on any case with a prescribed pressure boundary, Poiseuille-with-outlet-pressure being the obvious internal candidate.

## 6. Verification harness state

`tests/validation/test_n9_physical_matrix.cpp` is currently instrumented with temporary debug output that is not part of the intended test:

- `PROBE_OK` / `PROBE_FAIL` / `HIST` lines in `require_physical_convergence` (`:215-242`), with `HIST` dumping the last three iterations on failure.
- An `EQ` table of cross-algorithm max differences against `SIMPLE` for all three benchmarks (`:347-354`). A second, duplicate `EQ` block exists at `:414-421`.
- `LOC` (`:356-375`): worst cell and component for the `SIMPLE` vs `COUPLED` cavity difference with both values, plus the mean difference and the count of cells above 1e-3.
- `CAV` (`:406-412`): per-algorithm final history line for the cavity.
- A standalone `HOMOCOUNT` block (`:306-325`) that solves the cavity with `SIMPLE`, copies that state, runs 2 `COUPLED` block solves with `coupled_matrix_summary` on, and prints `dev_after_2_block_solves`.
- A second, fuller `HOMOCOUNT` block (`:377-404`) that starts from `cavity[0]` (the converged `SIMPLE` cavity), runs `COUPLED` to its own convergence, and prints the max deviation from `SIMPLE`. It sets `coupled_matrix_summary` twice (`:388`, `:390`).

**Why the first `HOMOCOUNT` block sits before the per-algorithm loop.** The loop over the six algorithms (`:327-345`) calls `require_physical_convergence` on Couette, Poiseuille and cavity for each. At algorithm index `k=5` (`COUPLED`) the cavity case fails the convergence gate, `require_physical_convergence` throws, and everything after the loop is unreachable: the `EQ` tables, the `LOC` block, the second `HOMOCOUNT` block, the `CAV` lines and the cross-algorithm equivalence gate itself (`:424-431`). The pre-loop block is therefore self-contained — it builds its own mesh and its own `SIMPLE` reference run — and is the only way to get the decisive measurement at all while the loop still throws. One correction to the earlier note: it is the **cavity** case that fails at `k=5`, not Couette.

### Measured state, current tree, `build-verify/test_n9_physical_matrix`

Reproduce with:

```
cmake --build build-verify --target test_n9_physical_matrix
./build-verify/test_n9_physical_matrix
```

```
HOMOCOUNT dev_after_2_block_solves=1.0461e-01
PROBE_OK   Couette alg=SIMPLE iters=451 cont=6.077917e-10
PROBE_OK   Poiseuille alg=SIMPLE iters=778 cont=3.664419e-11
PROBE_OK   cavity alg=SIMPLE iters=134 cont=2.693463e-11
... (SIMPLEC, PISO, PIMPLE, FRACTIONAL_STEP all PROBE_OK) ...
PROBE_OK   Couette alg=COUPLED iters=2 cont=2.685963e-11
PROBE_OK   Poiseuille alg=COUPLED iters=2 cont=9.976448e-12
PROBE_FAIL cavity alg=COUPLED iters=1500 cont=1.090042e-02 mom=5.165084e-10
HIST cavity COUPLED it=1500 cont=1.090e-02 mom=5.165e-10 momrel=5.071e-10 contnorm=1.837e-02 vel=0.000e+00 p=0.000e+00 pres=1.837e-02
terminate called after throwing an instance of 'std::runtime_error'
  what():  cavity did not converge
```

Four things in that output matter, and none of them is a gate result:

1. **The `COUPLED` block does not preserve the segregated fixed point.** `1.0461e-01` after 2 block solves started from `SIMPLE`'s converged cavity. The target for this probe is round-off, 1e-12 to 1e-10.
2. **The `COUPLED` cavity stalls at a frozen state, not at a slowly converging one.** `vel=0.000e+00` and `p=0.000e+00` mean the state is bit-identical from one outer iteration to the next, while `continuity_linf` stays at 1.09e-2. The block is at a fixed point of its own operator; the outer loop's continuity measure cannot be driven down because it is measuring something the block does not see. This is the signature predicted by section 3.
3. **The `COUPLED` linear solve itself is fine.** `momentum_residual` 5.17e-10, `momrel` 5.07e-10, and the `COUPLED_LINEAR_TRUE_RESIDUAL` line shows the true residual recomputed independently of GMRES at `relative=7.89e-10` against a requested 1e-9. The defect is in the operator that was assembled, not in the Krylov solve of it. Item 7 of section 2 is what allows `momrel` to reach this range at all.
4. **The linear benchmarks are not evidence.** `COUPLED` reports convergence on Couette and Poiseuille in 2 outer iterations (`converged_now` requires `iter > 1`, `:2525-2526`) with `continuity_linf` at 2.7e-11 and 1.0e-11. Those two cases are linear, the discrete solution is unique, and the block reaches it in one linear solve. This is the clearest illustration of why "each algorithm reports converged" is not the gate: on the linear benchmarks `COUPLED` looks perfect and is still wrong on the nonlinear one.

`EQ`, `LOC`, `CAV` and the second `HOMOCOUNT` block produce no output at all in this run, because the loop throws first.

**All of this instrumentation and both `HOMOCOUNT` blocks are temporary and must be removed or reverted before the work is committed.** The committed form of the file is `git show HEAD:tests/validation/test_n9_physical_matrix.cpp`. Nothing in the solver header may be committed while the test still carries the debug output, and the equivalence gate must not be weakened to accommodate a failing `COUPLED` (see `docs/development-and-validation-rules.md`, principle 1).


## 7. Second session: defects found and fixed

All of the following are in `src/cfdx/physics/steady_incompressible_solver.h`, and each was applied and measured, not merely derived. The reasoning for each is inline at the site.

### 7.1 Sign of the deferred non-orthogonal remainder (regression against `HEAD`)

`b(row) -= nonorth_flux` was **reversed relative to `HEAD`, which had it right**. The continuity row is `A x = b` with `b` holding the explicit part of `div(phi)` moved to the right-hand side, so the `- nonorth_flux` term of the reconstructed flux is added to `b`, not subtracted from it. Restored to `b(row) += nonorth_flux`.

### 7.2 Momentum response of the eliminated velocity unknowns (the section 3 defect)

The `gw` / `gconst` forms are now computed and consumed. Storage was changed from `std::map` per cell to a sorted, merged `std::vector<std::pair<std::size_t,double>>` per cell and component: O(nc * degree) entries, a few MB at production mesh sizes, instead of the hundreds of MB the map node overhead implied.

Coupling added to every continuity row:

- internal face: `rho*0.5*(G_d[owner] + G_d[opposite]).Sf`, both cells of the face, with the constant part produced by fixed-pressure boundary faces moved to `b`;
- boundary face **without** a prescribed velocity: `rho*G_d[owner].Sf`, weight 1, because such a face flux interpolates the owner value alone;
- boundary face **with** a prescribed velocity: nothing, the flux is fully explicit there.

### 7.3 Orientation of the face pressure coefficient

`D` and the Schur diagonal used `|Sf|` and the face normal. The face flux that continuity is measured on uses the centre-to-centre direction `e` and its orthogonal projection `Sf.e` (`rhie_chow_pressure_flux_internal`). The block now uses `Sf.e/d` and `e` squared in `rfn`, so the assembled coefficient is the exact derivative of that flux. The two forms agree only on an orthogonal mesh.

`Sf` must be the vector oriented **outward from the cell whose row is being assembled**, not the raw owner-oriented `face_area_vectors[f]`: on the neighbour side the raw projection is negative and `D` came out negative, which the face-coefficient guard rejects. The same applies to `Snon`.

### 7.4 Segregated pressure matrix was not the Jacobian of its own flux

The same `e` / `Sf.e` correction was applied to the segregated pressure matrix: `directional_face_coefficient` now takes the centre-to-centre unit vector instead of the face area vector, and the fixed-pressure boundary diagonal — which was missing the cell-volume factor entirely — was corrected to match the boundary face flux. The segregated pressure equation is now an exact Newton linearisation, so on a skewed mesh its residual floor comes only from the deferred non-orthogonal remainder instead of from an inconsistent matrix.

This changes nothing on an orthogonal mesh and it changes no converged answer; it removes a Jacobian error.

### 7.5 Equation relaxation was applied on only one path

`relax_momentum_equation` was called inside the segregated branch only, so `COUPLED` assembled its pressure operator from an unrelaxed `A_P` while every other algorithm used the relaxed one. Since `rAU = 1/A_P` feeds the Rhie-Chow operator, the two paths enforced different continuity operators by a factor of `alpha_u` (0.7 in the N9 controls). The relaxation now runs before either solve. The relaxed system shares the unrelaxed fixed point, so the segregated answers are unchanged.

### 7.6 `HbyA` on the coupled path

`mass_flux` was rebuilt from `U` alone right after the block solve, and the flux used for the continuity residual came from `build_hbya` evaluated with the predictor. Both measure continuity on the flux of `U`, a different operator from the one the block enforces. On the coupled path `HbyA` is now the reconstruction `U + rAU*V*grad(p)` at the corrected pressure, which is what the segregated path's `HbyA` equals by construction of its own velocity update.

### 7.7 Initial flux of a restarted solve (new finding, affects every algorithm)

`mass_flux` was initialised before the loop with `make_mass_flux`, the purely **convective** flux. A converged state is however the fixed point of an iteration whose momentum equations are assembled with the **conservative** Rhie-Chow flux. A solve started from a converged field — a restart, or a second algorithm seeded with the first one's answer — therefore began one operator away from the state it was handed. Before the fix, two `COUPLED` block solves from `SIMPLE`'s converged cavity state moved the field by 1.05e-1; after it, by 2.1e-2, and the pre-solve continuity residual dropped from 2.3e-5 to 8.5e-6.

The fix assembles the momentum equations once from the incoming state to obtain a provisional `rAU`, then rebuilds `mass_flux` as the conservative flux. On a cold start the provisional pressure flux vanishes and this reproduces the previous convective predictor exactly, so no cold-start case changes. The loop-invariant body-force and scalar-BC fields were hoisted out of the loop to make the provisional assembly possible.

## 8. What the measurements say now

Recorded on the 16x16 cavity, starting `COUPLED` from `SIMPLE`'s converged state, with the momentum response wired in:

| quantity | before wiring | after |
| --- | --- | --- |
| pre-solve continuity residual at `SIMPLE`'s state | 2.3e-5 | 8.5e-6 |
| pre-solve momentum residual | 7.7e-4 | 3.4e-4 |
| `max abs(U - U_SIMPLE)` after 2 block solves | 1.05e-1 | 2.1e-2 |

One fact is now established firmly, and it constrains everything else: **the assembled continuity row is the exact operator.** An instrumented comparison of `(A x - b)` against `div(make_rhie_chow_mass_flux(HbyA, p, rAU))` with `HbyA = U + rAU*V*grad p`, evaluated per cell over the whole block, matched to round-off (`worst = 1.2e-9`) on every cell except the gauge reference cell, where the eliminated pressure column makes the comparison meaningless by construction. Since both sides are affine in the unknowns, the block is not assembling the wrong thing.

## 9. Test-suite status of this tree

| test | `HEAD` | this tree |
| --- | --- | --- |
| `test_n9_physical_matrix` | aborts: "cavity did not converge" | aborts: "Couette did not converge" |
| `test_couette_quick` | pass | **fail**, `COUPLED/BlockSchur` and `COUPLED/MGR` only |
| `test_dual_time_navier_stokes` | pass | **fail**, physical step rejected after bounded retries |
| the other 145 tests | pass | pass |

The two new failures are confined to the `COUPLED` path and are the same defect as section 10, seen more clearly. `test_couette_quick` reports the signature directly: momentum relative residual 8.3e-11, `dU = 0`, `dp = 0`, `reconstructed_velocity_continuity = 2.5e-12`, and `corrected_flux_continuity = flux_velocity_mismatch = 0.0978`. The velocity field is divergence-free and momentum-converged; only the Rhie-Chow flux carries the residual. The block sits on a state it considers converged while the flux operator says otherwise.

Isolated by bisection: the failures are **not** caused by 7.7, which was disabled and re-tested. They come from one of the changes that affect only `COUPLED` — 7.1, 7.2, 7.3, 7.5 or 7.6.

## 10. Open defect: the block satisfies its system but the measured continuity residual does not

This is the remaining blocker and it is not yet explained.

After a converged block solve the outer loop rebuilds `mass_flux` from `HbyA = U + rAU*V*grad p` and reports `div(mass_flux)` as the continuity residual. That residual stays at 0.04 to 0.10 while the block's own rows are satisfied: `COUPLED_LINEAR_TRUE_RESIDUAL` reports `||b - A x||` relative below the requested 1e-9, in 20 to 34 GMRES iterations, and the post-solve gate re-derives that residual independently.

The two operators were shown to be the same affine function on the pre-solve state (section 8). Two affine functions that agree at one point agree everywhere, so one of the two premises must be wrong. Candidates, with the evidence against each:

- **The gauge reference cell.** `reference_value = 0`, the constraint row pins `p(0) = 0` exactly, and the dropped column contributes nothing. Ruled out for the cavity (all boundary fluxes vanish there, so `sum_c div_c = 0` forces `div_0 = 0` once the others are satisfied), and the worst cell is 249, not 0. Not fully ruled out for a channel with prescribed through-flow, where `sum_c div_c` equals the net boundary mass flux rather than zero.
- **Duplicate `(row, col)` entries.** The momentum response pushes each pressure column twice per face. `SparseMatrix::finalize` keeps duplicates unsummed, but every consumer sums them, and the pre-solve comparison was exact.
- **The frozen non-orthogonal remainder.** Identically zero on the cavity (orthogonal mesh), so it cannot explain a post-solve difference there.
- **`rAU`.** Both the block and the outer loop take it from the same relaxed diagonal.
- **The momentum rows.** The pre-solve momentum residual is 3.4e-4, and `momentum_equation_residual_relative` is measured against `result.reference_momentum_residual`, captured from the *pre-solve* state at iteration 1. A gate that compares a residual against itself can report convergence while the absolute residual is orders of magnitude above tolerance. If the converged `SIMPLE` state satisfies the momentum equations only to 3.4e-4 absolute, a block that solves them exactly *must* move the field, and the 2.1e-2 deviation is that move rather than an assembly error. That would make section 10 and the `HOMOCOUNT` gate a question about what "the same fixed point" means for an algorithm whose momentum solve is not under-relaxed, not about the continuity operator.

## 11. Ordered next steps

1. Settle section 10 before touching anything else. The cheapest discriminator: run the same `COUPLED` case with 7.5 disabled (no momentum relaxation in the block). If the measured residual goes to zero, the momentum rows and the pressure operator are inconsistent with each other inside the block and 7.5 was masking it rather than causing it. If it does not, the disagreement is in the flux reconstruction and the post-solve `(A x - b)` versus `div(phi)` comparison has to be run at the block's solution rather than at its input.
2. Report the reference cell explicitly. `corrected_flux_continuity_linf` currently mixes a constrained maximum with one unconstrained cell. Either exclude `pressure_reference_cell` from the reported continuity measures, or enforce its divergence by moving the gauge constraint onto the pressure level instead of onto a continuity row. Decide this before any convergence gate is trusted.
3. Make `coupled_matrix_summary` usable: drop the unused `rau3` and `wd`, fix the literal `\n` in the `COUPLED_PRECONDITIONER` line, and drop `rhie_chow_schur_vs_algebraic_diag_delta` — the supplied diagonal is overwritten with 1.0 at the reference cell before the preconditioner call, so that metric is contaminated by the gauge row.
4. Re-run the skew campaign. Neither `docs/validation/PHASE3_6_SKEW_CAMPAIGN.md` (a Laplacian operator sweep, no pressure-velocity coverage) nor `docs/validation/LEVEL_C_COUPLED.md` (turbulence/radiation/CHT) covers `COUPLED`. The only skew exposure in the repository is `make_channel_mesh(12,16,0.25)` inside the N9 test. Sections 7.3 and 7.4 are exactly what skewed meshes expose.
5. Once `COUPLED` converges, restore the full cross-algorithm matrix and check the segregated paths for regressions from 7.4, the skew cases in particular, since that is the only place the Jacobian change bites.

Status: in progress. Sections 7.1 to 7.7 are fixed and measured; section 10 is open. The test suite has two new `COUPLED` failures, section 9. `tests/validation/test_n9_physical_matrix.cpp` is back to its committed form.
