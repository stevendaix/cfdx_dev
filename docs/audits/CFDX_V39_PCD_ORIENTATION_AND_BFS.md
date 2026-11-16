# CFDX V39 — PCD pressure commutator orientation, Schur mismatch, production A/B

**11 November 2026.** The exact `cfdx_dev-physics-bfs-v6-hilo-rhie-chow-production.zip` remains the **cumulative-patch base**, via V37 → V38 → V39. V39 is a **preconditioner-only, OFF-by-default** research extension to V38. No FV equation, Rhie–Chow flux, SST closure, conservation gate, physical tolerance or pressure reference is changed.

## 1. Primary finding on captured matrices

V38 recorded the **actual** production-coupled `A,b,x0` and PCD `Mp,Kp,Fp` from the **same** first linear solve, on a 152-cell debug BFS mesh (608 unknowns). It formed the frozen Schur complement

\[
S=C-D A_{uu}^{-1}G,\qquad r_S=b_p-D A_{uu}^{-1}b_u.
\]

The operator `S` is NOT the full Jacobian of `U,p,k,omega`, and the 152-cell synthetic/debug BFS is NOT the experimental 20,540-cell V4 mesh.

| Read-only algebraic inverse action | Frobenius `||SP-I||_F/sqrt(n)` | Best-scaled physical-state Schur-RHS residual |
| --- | ---: | ---: |
| Historical single inverse, `-Fp^+` | 864.987 | 0.659912 |
| Existing opt-in `-Mp^{-1} Fp Kp^{-1}` | 28.6036 | 0.850988 |
| **V39 opt-in** `-Kp^{-1} Fp Mp^{-1}` | 8.15950 | 0.871716 |
| Sign-altered reverse **oracle only** `+Kp^{-1} Fp Mp^{-1}` | 7.79885 | 0.871716 |
| Exact dense inverse of `C-D diag(Auu)^{-1}G` **oracle only** | 1.31262 | 0.570238 |
| Exact dense inverse of `C-D blockdiag_3x3(Auu)^{-1}G` **oracle only** | 1.31262 | 0.570238 |
| True `S^+` **reference only** | 2.894e-13 | ~1e-13 |

The two PCD orders are not commutative, and **smaller Frobenius defect does not imply a smaller true production linear residual**. The diagonal and 3×3 surrogate pressure inverses are diagnostic, never substituted for the actual production BlockLocal/SIMPLE/SIMPLEC algorithms. The sign comparison is diagnostic only: alternate sign is **not** exposed as a production flag.

**Rank mismatch:** `rank(Fp)=151/152`, but `rank(S)=152/152`. Any product containing the singular `Fp` as a central direct multiplication has rank at most 151. Thus a full-rank Schur inverse cannot be recovered simply by making the `Kp` inversion more accurate. This is necessary context, not an excuse to relax inner/outer tolerances.

**Dominant mismatch rows:** `||S_i||_2/||Fp_i||_2` is ~3696.4 at cell 0 and ~2307.2 at cell 12. Recovered from native case geometry, these are near the downstream face of the step: their representative centers are both `x=0.0211667 m` and `y=0.003175 m` / `y=0.009525 m` in this debug grid. This suggests a high-priority investigation of near-wall/inflow auxiliary pressure rows and gauge treatment, **not** a validated causal attribution to a particular boundary condition.

## 2. Controlled production A/B on 152 cells

V39 adds **one** new `PcdSchurActionMode` enum value and a case flag:

```
pcd_reverse_lumped_mass_v39_experimental=true
```

The original default, legacy `-Fp^{-1}`, remains untouched; the existing `pcd_lumped_mass_commutator_experimental=true` remains a distinct, mutually exclusive flag. The new flag is accepted only for COUPLED/PCD; simultaneous commutator flags are rejected. Both compositions use a **strictly positive diagonal** mass and check the CURRENT `Kp` inverse residual at `<= 1e-9 * ||rhs|| + 1e-12`, without loosening the outer coupled linear criterion. The production `Kp` nested tolerance is tightened to `min(original,1e-10)` for **both** commutator variants. On solver rejection the pressure output is transactional and the coupled outer step remains rejected.

| Metric, same frozen production `A,b,x0` | Original opt-in commutator | V39 reverse opt-in |
| --- | ---: | ---: |
| SHA-256 `A`, `b`, `x0` | Identical | Identical |
| FGMRES iterations | 2000 | 2000 |
| Independent `||Ax-b||/||b||` | **0.006230732** | **0.006346720** |
| Outer state accepted | No | No |
| Physical global `U,p,k,omega` converged | No | No |

The reverse mode **does not improve the production coupled solve** on this test; it must remain OFF by default. A preliminary reverse run failed before its first Krylov step because its Kp tolerance had not been tightened like the original commutator; that integration defect was corrected, and the table above reports the corrected run.

## 3. Tests and reproducibility

New C++ regression `test_pcd_reverse_commutator_v39` checks the noncommuting operator order, setup validation, positive lumped-mass contract, output transactionality, and independently detects a stale `Kp` callback after numerical updates. Python `test_diagnose_pcd_orientations_v39.py` verifies exact Schur and both orders on synthetic matrices, strict optional algebraic gates, valid/invalid HDF5 geometry capture, and no physical-qualification claims. Production build and selected surrounding CTests are explicitly recorded in the delivery report.

The new bounded diagnostic (SciPy and NumPy) is independently reproducible:

```bash
python scripts/diagnose_pcd_orientations_v39.py \
  --coupled-A capture_old_A.mtx --coupled-b capture_old_b.mtx \
  --mp pressure_state1_mp.cfp \
  --kp pressure_state1_kp.cfp --fp pressure_state1_fp.cfp \
  --mesh bfs_pcd_old.cfdx.h5 \
  --output pcd_orientation_v39.json
```

Dense operations are bounded to **256 pressure cells by default**, hard maximum **1024** with explicit override. The tool neither writes matrices back to production nor uses the `S` oracle as a solver preconditioner. The CTest Python audit requires `CFDX_ENABLE_PCD_V38_DENSE_AUDIT_TESTS=ON` plus NumPy and SciPy; the new C++ regression remains part of ordinary tests. The optional mesh-center enrichment requires `h5py`.

## 4. Next numerical priorities

1. Compare actual `S` near-step rows 0/12 with the **same-row** `Fp`, `Kp`, and `Mp` stencils; classify exact physical face groups and auxiliary BC treatment. Avoid changes to the NASA physical BC based on a debug mesh alone.
2. Evaluate scalable approximate Schur operators `C-D diag(Auu)^{-1}G` and `C-D blockdiag_3x3(Auu)^{-1}G` in independent oracle tests on different meshes before implementing a **clearly named** hybrid. The dense inverse numbers above do not represent a production AMG solve.
3. Obtain successive **accepted** nonlinear states of `U,p,k,omega` (V34 only observed rejected SST states), with a physical residual recomputed after every candidate step; use identical matrices for all six Schur variants.
4. Replay 712 then 20,540 V4 BFS. Preserve separate levels: (a) execution completed, (b) global residuals converged with no threshold relaxation, (c) NASA Driver–Seegmiller validation of `x_r/h`, `Cf/Cp` and velocity/Reynolds-stress profiles.
5. Complete OpenFOAM Foundation native wall-function oracle (`updateCoeffs`, `manipulateMatrix`, `G`, `nut`, `omega`) independently of any Schur acceleration.

## Bibliographic controls

- Elman & Tuminaro (2009), *Boundary conditions in approximate commutator preconditioners for the Navier–Stokes equations*, ETNA 35, 257–280: https://etna.ricam.oeaw.ac.at/volumes/2001-2010/vol35/abstract.php?pages=257-280
- Firedrake `PCDPC` documentation describes an `K^{-1} F_p M^{-1}` arrangement, while explicitly warning of limitations for inflow/outflow boundary conditions: https://www.firedrakeproject.org/_modules/firedrake/preconditioners/pcd.html
- PETSc `PCFIELDSPLIT` Schur factorization and user-controlled pressure preconditioner: https://petsc.org/main/manual/ksp/

These are **methodological comparisons**, not proofs that a finite-element-oriented PCD formula should replace the finite-volume CFDX default. No GitHub push; no production-default change; full BFS convergence remains open.
