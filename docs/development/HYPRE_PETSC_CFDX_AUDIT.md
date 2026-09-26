# HYPRE and PETSc methods applicable to CFDX

## Purpose and boundary

This audit identifies algorithms and architecture patterns from HYPRE and
PETSc that can improve CFDX without adding either project as a mandatory
dependency. It is an implementation roadmap, not a claim that the native CFDX
classes are API-compatible replacements.

The audit follows three rules:

1. Reimplement documented numerical methods behind CFDX interfaces; do not
   copy large upstream modules or preserve their internal API.
2. Keep setup, apply, convergence, memory location and execution policy
   explicit. A requested GPU path must never silently run on the CPU.
3. Accept an optimization only with true-residual, conservation, numerical
   equivalence and measured-performance evidence.

The source review is pinned to:

- HYPRE commit `5b582610a113ede9162a5744d33873b860498a11`;
- PETSc commit `74bbc03299ffb37f3e3b88d943fcbbe2cdac1b1b`.

These projects remain the authority for their own implementations. CFDX uses
the algorithms and design lessons, not their names as evidence of equivalent
capability.

## Current CFDX baseline

CFDX already has CSR matrices, matrix-free operators, CG, BiCGStab, restarted
GMRES, Jacobi, Gauss-Seidel, ILU(0), Chebyshev smoothing, mixed-precision
helpers, a native Galerkin AMG V-cycle, native two-field splitting, MPI
reductions and halo plans, CUDA gradient/divergence kernels, and a CUDA
Jacobi-CG Poisson path.

The important gaps are integration and lifecycle rather than missing class
names:

- most pressure/momentum solves rebuild temporary state instead of reusing a
  symbolic pattern, AMG hierarchy or device allocation;
- the steady incompressible path does not select a preconditioner by equation
  type and does not expose a complete solver context;
- MPI reductions are not systematically fused or overlapped with useful work;
- the CUDA Poisson solve allocates and transfers on each call, so it is not a
  device-resident CFD solve;
- local CSR pressure systems can use an explicit orthonormal null-space
  projector and projected CG, but existing CFD algorithms still use a
  fixed-cell gauge and distributed pure-Neumann Poisson is not yet connected;
- the coupled velocity-pressure preconditioner has no scalable LSC/MGR-style
  Schur approximation.

## HYPRE inventory

| HYPRE facility | CFDX relevance | Decision |
|---|---|---|
| BoomerAMG | Highest for pressure Poisson, diffusion and scalar transport | Continue the native implementation: parallel coarsening, interpolation truncation, complexity controls, hierarchy reuse and GPU kernels |
| MGR | Highest for coupled velocity-pressure, energy and multiphysics block systems | Implement a native reduction preconditioner over explicit variable groups; use AMG on the reduced pressure/scalar operator |
| ILU: ILU(k), ILUT, RAS and Schur variants | High for nonsymmetric momentum/transport and local MPI subdomains | Extend ILU(0) first with ILUT and local RAS; avoid a monolithic port |
| FSAI | High on GPU because apply is sparse matrix-vector work and avoids triangular solves | Implement after persistent GPU CSR; use an adaptive sparse pattern and SPD local dense solves during setup |
| Hybrid | High as a policy: try cheap diagonal Krylov, monitor convergence factor, build AMG only when needed | Implement as a solver-selection policy, not a separate matrix type |
| FlexGMRES | High when AMG, inner solves or nonlinear preconditioners vary per iteration | Add native FGMRES before variable MGR/Schur preconditioning |
| LGMRES | Medium-high for repeated slowly converging nonsymmetric solves | Add augmentation/recycling after stable GMRES workspace reuse |
| GMRES, PCG, BiCGSTAB | Already present in CFDX | Harden convergence reasons, breakdown checks, true residual and MPI reduction policy |
| ParaSails | Algorithmically useful sparse approximate inverse, but upstream recommends FSAI | Do not port; implement FSAI instead |
| Euclid and PILUT | Useful historical parallel ILU approaches | Mine ordering/drop-policy ideas only; native ILUT/RAS has better scope |
| SMG and PFMG | Excellent for logically rectangular diffusion grids | Defer until CFDX has a structured-grid specialization; unsuitable as the generic unstructured path |
| SysPFMG, Split and SSAMG | Semi-structured block-grid specialization | Defer; MGR and FieldSplit fit the current unstructured storage better |
| AMS | Solver for edge-element H(curl) Maxwell systems | Out of scope until CFDX adds an electromagnetic edge-element model |
| ADS | Solver for face-element H(div) systems | Out of scope for the current cell-centred finite-volume equations |
| LOBPCG | Block eigensolver | Defer until stability, modal or spectral-analysis requirements exist |
| IJ/ParCSR interface | Strong model for distributed ownership and off-process assembly | Adopt the ownership/assembly concepts in CFDX CSR, without reproducing the API |
| Memory/execution policy | Essential for portable CPU/GPU behavior | Unify CFDX host/device ownership, explicit policy and no-fallback behavior |
| Runtime mixed precision | Useful for GPU throughput and AMG smoothers | Use FP32 operator/preconditioner work with FP64 reductions and periodic true-residual refresh |

### HYPRE source findings that matter

The reviewed BoomerAMG setup separates graph strength, C/F selection,
interpolation, Galerkin products, smoothing and the coarse solve. CFDX should
preserve those stage boundaries so each can be tested and moved to the GPU
independently. The current native direct interpolation is a valid first stage,
but it is not HMIS/PMIS plus extended+i interpolation.

`par_mgr_setup.c` builds levels from user-supplied variable groups, configures
F-point relaxation, transfer operators and a separate final coarse solver.
`par_mgr_solve.c` applies residual restriction, F-relaxation and coarse-grid
correction. This maps naturally to a CFDX coupled matrix split into velocity,
pressure, energy and species groups.

`par_ilu_setup.c` makes the useful design split explicit: local factorization,
optional interior/interface permutation, and a Schur solve for interface
unknowns. CFDX should first implement local ILUT/RAS and reuse its sparsity
pattern; a distributed Schur level comes later.

`par_fsai_setup.c` grows a sparse inverse factor pattern adaptively and solves
small SPD dense systems during setup. The expensive setup is acceptable only
when reused across several nonlinear iterations or time steps. Its application
is attractive on GPUs because it avoids serial triangular substitution.

`amg_hybrid.c` monitors the average residual convergence factor before paying
for AMG setup. That policy is particularly useful for changing CFD systems in
which some momentum or scalar equations are already strongly diagonally
dominant.

## PETSc inventory

| PETSc facility | CFDX relevance | Decision |
|---|---|---|
| KSP solver context | Highest: uniform setup/solve/reuse/monitor contract | Add a native reusable linear-solver context with explicit operator and preconditioner generation counters |
| FGMRES | Highest for varying AMG/MGR/Schur applications | Implement before advanced block preconditioners |
| PIPEFGMRES/PIPEGCR/PIPECG | Highest at MPI scale | Introduce asynchronous fused reductions and overlap them with preconditioner/SpMV work; retain a stable non-pipelined fallback |
| LGMRES/GCRODR-style recycling | High for successive pressure and momentum systems | Start with LGMRES augmentation; add subspace recycling only after deterministic restart tests |
| PCFieldSplit | Already partially implemented | Add multiplicative and symmetric multiplicative modes, scalable inner solvers and nested block operators |
| Schur factorization | Already has diag/lower/upper/full forms | Replace dense inner block solves for production systems and expose Schur approximation policy |
| PCLSC | Highest for coupled incompressible Navier-Stokes | Implement a native least-squares commutator pressure Schur approximation using existing block operators |
| PCGAMG/PCMG | High | Add smoothed aggregation as a second AMG family, near-null-space input, interpolation reuse, complexity limits and process reduction on coarse levels |
| ASM/GASM/RAS and block Jacobi | High for MPI-local smoothing and additive domain decomposition | Implement overlap-aware RAS over the existing halo/partition plan |
| BDDC/HPDDM | Powerful but large domain-decomposition frameworks | Defer; do not recreate without a demonstrated scaling need and dedicated design |
| Jacobi/SOR/ILU/ICC/Chebyshev | Mostly present | Add block Jacobi, ILUT/ICC where matrix properties permit, and robust eigenvalue estimation for Chebyshev |
| MatShell | Already reflected by CFDX matrix-free operators | Complete solver/preconditioner support for matrix-free operators and explicit state invalidation |
| MatNest | Very relevant to coupled multiphysics | Add a lightweight native block-operator container; avoid extracting monolithic submatrices |
| MatNullSpace/NearNullSpace | Highest for pressure and AMG robustness | Add constant pressure projection and optional near-null-space vectors instead of relying only on a pinned row |
| COO GPU assembly | Highest for repeated device assembly | Precompute row/column locations and update only values on GPU; retain CSR for solve operations |
| VecScatter/PetscSF concepts | High for halo exchange | Represent communication as a reusable graph and overlap begin/end exchange with owned-cell work |
| SNES Newton line-search/trust-region | High for future fully coupled nonlinear solves | First add inexact Newton and backtracking around CFDX residual/Jacobian operators; do not port all SNES methods |
| Anderson/nonlinear GMRES | Medium-high for segregated fixed-point acceleration | Evaluate after SIMPLE/PIMPLE residual semantics are reliable |
| TS BDF/ARKIMEX/adaptivity | High for transient stiff multiphysics | Reuse method concepts in the native time integrator later; current steady/M1 work remains first |
| DM/DMPlex/DMStag | Broad mesh and discretization framework | Do not reproduce; CFDX already owns mesh/FVM data structures. Adopt only local/global ownership and ghost-update contracts |
| Options, monitors and convergence reasons | High operational value | Add typed runtime solver configuration, monitors and explicit converged/diverged reasons |
| Logging/stages/events | High for optimization evidence | Instrument assembly, halo, SpMV, preconditioner, reductions, transfers and kernels separately |
| TAO and adjoints | Useful for optimization, not current solver closure | Defer until sensitivity/optimization requirements are defined |

### PETSc source findings that matter

`fieldsplit.c` applies block additive, multiplicative and Schur factorizations
through independently configured inner solvers. CFDX has the factorization
order but must replace the current small dense solve before calling it
scalable.

`lsc.c` constructs or reuses the pressure Laplacian-like product and applies
the inverse Schur approximation as two pressure solves around momentum and
gradient/divergence products. This is a better target for the coupled 4N
solver than a diagonal Schur approximation.

`pipefgmres.c` and `pipegcr.c` begin a multi-dot reduction, launch a nonblocking
global reduction, perform preconditioner and matrix work, then finish the
reduction. CFDX's existing reduction packet is only the first half of this
design: useful computation must occur while the collective is in flight.

`gamg.c` reuses interpolation and Galerkin product structure when the nonzero
pattern is unchanged, and can reduce the number of active MPI ranks on coarse
levels. `agg.c` preserves near-null-space vectors while filtering the
prolongator. Both ideas are important for repeated pressure solves and large
MPI runs.

PETSc's GPU design keeps matrix and vector objects resident in the selected
memory space and treats host/device synchronization as object state. Efficient
GPU assembly uses a preallocated COO pattern. CFDX should adopt those lifecycle
rules rather than adding isolated CUDA entry points that allocate and transfer
for every solve.

## Pressure-velocity algorithms

The upstream methods improve each CFDX coupling algorithm differently:

| CFDX algorithm | Native methods to apply |
|---|---|
| SIMPLE | Reuse momentum/pressure symbolic setup and AMG hierarchy; hybrid cheap-to-AMG policy; pressure null-space projection; equation-specific convergence reasons |
| SIMPLEC | Same as SIMPLE, while preserving the consistent velocity and face-flux correction; use FGMRES when the effective preconditioner changes |
| PISO | Reuse one momentum predictor across genuine pressure/neighbor corrections; choose corrector count from continuity progress; avoid rebuilding unchanged sparse structure |
| PIMPLE | Combine outer nonlinear monitoring/inexact linear tolerances with inner PISO corrections; consider line search or pseudo-transient continuation only after residual consistency |
| Fractional step | Use a strong reusable pressure Poisson AMG, explicit projection/null-space treatment and device-resident gradient/divergence/update kernels |
| Coupled 4N | Nested block operator plus FieldSplit/MGR; AMG-preconditioned velocity blocks; LSC or pressure-Laplacian Schur approximation; FGMRES outer iteration |

The three scalar momentum component solves are independent on the segregated
CPU path and may execute concurrently above a measured size threshold. On a
GPU, three host threads launching independent full solves are not automatically
faster: the preferred design is persistent data plus streams or batched/fused
component kernels, selected from profiling.

Adaptive pressure-corrector counts may reduce work for PISO/PIMPLE/fractional
step, but must be disabled by default until tests prove that final momentum,
continuity and analytical errors are unchanged.

## GPU plan

GPU support is an execution architecture, not a CUDA version of one Poisson
function. The implementation order is:

1. Introduce persistent device mirrors for CSR pattern, values, vectors and
   solver workspace. Track host/device generations and synchronize only stale
   data.
2. Provide device vector kernels, fused dot/norm reductions, CSR SpMV and
   Jacobi/Chebyshev smoothers. Use vendor sparse primitives only behind an
   optional backend selected by benchmark.
3. Keep the AMG hierarchy resident and add restriction, prolongation, residual
   and smoothing kernels. Build the hierarchy on CPU initially, copy it once,
   then evaluate GPU PMIS/coarse setup only when setup time dominates.
4. Add preallocated COO/value-only assembly for repeated FVM operators and
   persistent pressure-velocity fields.
5. Move Rhie-Chow flux, velocity correction, pressure projection and residual
   evaluation to the device so each nonlinear iteration does not bounce data
   through the host.
6. Add FSAI as the GPU-friendly nonsymmetric/local preconditioner. Keep ILU
   optional because triangular solves can expose less GPU parallelism.
7. Overlap CUDA streams, halo transfer and MPI collectives. CUDA-aware MPI is
   an explicit capability; stage through pinned host buffers otherwise.
8. Add FP32 preconditioner/operator modes with FP64 accumulation, reliable
   updates and periodic FP64 true-residual recomputation.

Required GPU evidence includes transfer bytes, allocation count, kernel and
wall time, setup/apply split, iteration count, true residual, mass imbalance,
and CPU/GPU solution error. Small problems should remain on CPU when transfer
and launch overhead wins.

## Prioritized implementation slices

### Implemented foundation

- reusable CG/GMRES contexts with explicit full-setup, numeric-refresh and
  unchanged-operator counters;
- native AMG numeric hierarchy refresh that preserves interpolation and the
  Galerkin sparsity structure when the CSR pattern is unchanged;
- native smoothed-aggregation AMG with deterministic aggregates, a
  constant-mode tentative prolongator, damped-Jacobi interpolation, Galerkin
  levels and the same reusable numeric-refresh lifecycle;
- one pressure CG/preconditioner context per segregated incompressible solve,
  shared by SIMPLE, SIMPLEC, PISO, PIMPLE and fractional-step corrections;
- typed equation-specific Krylov/preconditioner selection;
- explicit constant or multi-vector null-space projection for local CSR
  systems, including basis orthonormalization, operator-mode validation,
  incompatible-RHS rejection and a canonical projected-CG gauge;
- `NullSpaceModel::Constant` selection for pressure/diffusion problems, with
  conservative Jacobi policy until AMG near-null-space propagation is
  implemented.

The null-space API is currently host-local. MPI and GPU capability flags remain
disabled; the distributed pressure solver and CUDA Poisson path must be wired
and validated before those capabilities are advertised.

### P0: remaining correctness and reusable lifecycle

- near-null-space input and propagation through AMG levels;
- quiet-by-default monitors, explicit convergence/divergence reasons and true
  residual refresh;
- reusable contexts for momentum, scalar transport and coupled Schur solves;
- replace explicit CSR snapshots with matrix mutation generations once direct
  coefficient access is encapsulated.

### P1: scalable segregated solvers

- safe concurrent CPU momentum-component solves;
- adaptive corrector policy for PISO/PIMPLE/fractional step behind an explicit
  control;
- FGMRES and LGMRES augmentation;
- fused/nonblocking MPI reductions and halo overlap;
- ILUT plus overlap-aware local RAS.

### P2: coupled CFD preconditioning

- nested block operator without submatrix extraction;
- native LSC pressure Schur approximation;
- MGR-style variable reduction with AMG coarse solve;
- scalable inner Krylov methods in FieldSplit;
- Couette, Poiseuille, cavity and manufactured coupled validation.

### P3: device-resident execution

- persistent CUDA sparse/vector/workspace objects;
- GPU SpMV, vector reductions and AMG transfer/smoothing;
- device-side pressure-velocity update chain;
- optional FSAI and mixed precision;
- CUDA-aware MPI/pinned staging with measured overlap.

### P4: specialized extensions

- smoothed-aggregation AMG with near-null-space preservation;
- structured semicoarsening only for a future structured-grid backend;
- inexact Newton with line search and nonlinear acceleration;
- transient BDF/ARKIMEX methods when stiff transient validation is ready.

## Explicit non-goals

- No wholesale source copy from HYPRE or PETSc.
- No claim that `NativeBoomerAMGPreconditioner` is BoomerAMG-equivalent.
- No mandatory PETSc, HYPRE, CUDA, HIP, SYCL or vendor sparse dependency.
- No GPU success status after a CPU fallback.
- No port of AMS, ADS, LOBPCG, BDDC, HPDDM, DM, SNES or TS as complete
  frameworks without a concrete CFDX requirement.
- No performance claim without reproducible hardware and numerical evidence.

## Acceptance matrix

Every implementation slice must include:

| Area | Evidence |
|---|---|
| Numerical | independent `b - Ax`, finite values, breakdown reason, iteration history |
| CFD | momentum residual, continuity, boundary/internal flux balance, pressure gauge invariance |
| Validation | Couette, Poiseuille, cavity and a manufactured or matrix-level test appropriate to the slice |
| Parallel | serial/thread/MPI/GPU equivalence at documented tolerances |
| Performance | setup and apply time, reductions, halo time, transfers, memory, speedup and efficiency |
| Robustness | sanitizer build, deterministic mode where promised, no hidden fallback |
| Traceability | source commit, equations used, code/test links and known limitations |

## Upstream references

- [HYPRE solvers and preconditioners](https://hypre.readthedocs.io/en/latest/ch-solvers.html)
- [HYPRE BoomerAMG](https://hypre.readthedocs.io/en/stable/solvers-boomeramg.html)
- [HYPRE MGR](https://hypre.readthedocs.io/en/latest/solvers-mgr.html)
- [HYPRE ILU](https://hypre.readthedocs.io/en/latest/solvers-ilu.html)
- [HYPRE mixed precision](https://hypre.readthedocs.io/en/latest/ch-mprecision.html)
- [PETSc KSP and multigrid](https://petsc.org/main/manual/ksp/)
- [PETSc FieldSplit](https://petsc.org/main/manualpages/PC/PCFIELDSPLIT/)
- [PETSc LSC](https://petsc.org/main/manualpages/PC/PCLSC/)
- [PETSc matrices and GPU COO assembly](https://petsc.org/main/manual/mat/)
- [PETSc nonlinear solvers](https://petsc.org/main/manual/snes/)
- [PETSc time integration](https://petsc.org/main/manual/ts/)
- [PETSc GPU programming model](https://petsc.org/main/manual/getting_started/)

The work is tracked under issues #34, #390, #407 and #416. Each implementation
PR should close a bounded checklist item instead of treating this audit as one
large port.
