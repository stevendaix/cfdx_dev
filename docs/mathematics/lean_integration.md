# CFDX Lean 4 mathematical foundations

## Purpose
CFDX uses Lean 4 and Mathlib as a separate formal-methods project for
machine-checked mathematical specifications and proofs. The production solver
remains C++ (with Python for workflows and analysis). Lean is not a CMake,
runtime, MPI, or HPC dependency.

The intended development chain is:

1. Continuous physical and mathematical formulation.
2. Explicit discrete numerical specification.
3. Lean definitions, hypotheses, and machine-checked theorems.
4. C++ implementation of the specified operators.
5. Bounded conformance tests linking selected C++ outputs to the specification.
6. Independent numerical verification and CFD qualification.

A Lean theorem about an abstract operator does not, by itself, prove that the
production C++ implementation constructs that operator. Tests establish
evidence on tested inputs; they are not a universal proof of all C++ executions.

## Layout

- `formal/lean/lean-toolchain`: pinned Lean compiler version.
- `formal/lean/lakefile.lean`: pinned Mathlib release.
- `formal/lean/lake-manifest.json`: committed lockfile pinning exact revisions
  of Mathlib and all transitive Lake dependencies.
- `formal/lean/CFDX/`: CFDX-owned formal definitions and theorems.
- `scripts/check_lean_placeholders.py`: rejects explicit `sorry`, `admit`,
  and local `axiom` declarations in CFDX-owned Lean source.
- `.github/workflows/lean-verification.yml`: independent GitHub Actions job.

## Local build

Install the Lean toolchain manager (elan), then run:

```bash
cd formal/lean
lake exe cache get
lake build
cd ../..
python3 scripts/check_lean_placeholders.py
```

The Lean toolchain, Mathlib release, and transitive dependency revisions are
pinned. Ordinary builds must use the committed `lake-manifest.json` and must
not run `lake update`. A deliberate dependency update must regenerate and
review the manifest, update the relevant pins, and pass the complete Lean
workflow.

## Initial mathematical contract

`CFDX.LinearAlgebra.SPD` defines a finite-dimensional real matrix contract
with symmetry and strict positivity of its quadratic form for every nonzero
vector. The theorem `SPDMatrix.quadratic_form_pos` exposes the positivity
property as a kernel-checked result.

This first theorem is deliberately modest: it validates the formal project
and establishes a vocabulary for later specifications. It does **not** prove
that any CFDX production matrix is SPD, that CG converges for a production
operator, or that a finite-precision distributed solve converges.

## CI policy and limitations

The Lean workflow:

- builds the standalone Lean project against the committed Lake lockfile;
- checks that the lockfile remains unchanged during the build;
- checks CFDX-owned Lean files for explicit proof placeholders and local
  `axiom` declarations;
- does not alter or replace the existing C++/Python/CFD validation workflows.

The source guard is a narrow policy check, not a complete Lean security audit.
Imported Mathlib theorems may rely on foundational axioms accepted by Lean's
logical environment. Future work should add an explicit audit of theorem
assumptions (for example with `#print axioms`) for designated critical
results, with a documented allowlist.

## First finite-volume conservation theorem

`CFDX.FiniteVolume.OrientedMesh` represents a finite cell/face incidence
relation using signed real coefficients. Its explicit
`internal_face_balance` hypothesis states that, for each internal face, the
sum of its signed incidence coefficients over all cells is zero.

`CFDX.FiniteVolume.internal_face_fluxes_cancel` proves that summing the cell
balances over all cells cancels the contribution of every internal face.
The proof is exact and algebraic: it exchanges the finite cell/face sums,
factors out the one shared face flux, and applies the incidence cancellation
hypothesis.

This theorem deliberately has a narrow scope:

- The same scalar face flux is shared by every cell incident to a face.
  Independently evaluated left/right fluxes need an additional equality
  hypothesis or a conservative numerical flux construction.
- Boundary faces are excluded. Their contribution is not zero; it belongs in
  the global boundary-flux balance.
- The theorem does not yet derive incidence cancellation from a concrete
  mesh data structure, face orientation, or cell adjacency representation.
- It does not prove consistency, accuracy, stability, positivity, pressure-
  velocity coupling, or convergence of a CFDX discretization.

The next formal step is to derive the incidence hypothesis from the actual
CFDX mesh topology and face ownership/orientation rules, then compare the
formal cell-face balance with the C++ flux assembly on deterministic meshes.

## Integration roadmap

1. **Infrastructure:** pinned Lean/Mathlib, committed transitive lockfile,
   independent CI, initial algebraic contract, and documentation.
2. **Finite-volume foundations:** oriented faces, incidence signs, fluxes,
   internal-face cancellation, and global conservation under explicit mesh
   assumptions.
3. **Discrete operators:** actual CFDX gradient/divergence pairings, boundary
   conditions, adjointness, symmetry, positivity, and null spaces.
4. **Linear and block solvers:** conditional contracts for CG/PCG, GMRES,
   Schur complements, and preconditioner properties.
5. **Implementation correspondence:** trace mathematical definitions to C++
   assembly routines and add deterministic small-mesh conformance tests.
6. **Research:** anisotropic-mesh operators and preconditioners, pressure-
   velocity coupling, error estimation, and separate numerical qualification.

Each stage must state its hypotheses and distinguish proved properties,
tested properties, assumptions, and open questions. No solver tolerance,
physical model, convergence gate, or qualification result may be weakened to
make a formal or numerical check pass.
