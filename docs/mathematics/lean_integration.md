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

The owner/neighbour model also proves
`CFDX.FiniteVolume.global_cell_flux_balance_eq_boundary_flux`: summing the
complete owner/neighbour face balances over all cells leaves exactly the
outward owner-oriented boundary flux. Internal faces cancel; boundary faces
are retained once. This is the discrete global balance identity corresponding
to the C++ conservation audit, subject to the topology and shared-flux
assumptions above.

The next formal step is to derive the incidence hypothesis from the actual
CFDX mesh topology and face ownership/orientation rules, then compare the
formal cell-face balance with the C++ flux assembly on deterministic meshes.

## Mapping to the current CFDX mesh and flux code

A source audit of the current C++ implementation identifies a direct
owner/neighbour representation:

- `src/cfdx/core/mesh/ownership.h`: each face has one owner and either one
  neighbour or the boundary sentinel `-1`. Validation rejects out-of-range
  cells, unknown negative sentinels, and self-connected internal faces.
- `src/cfdx/core/mesh/mesh.h`: topology validation checks owner/neighbour
  consistency, the expected number of cell-face references, and that each
  face is referenced by its owner and (for internal faces) its neighbour.
- `src/cfdx/core/numerics/divergence.h`: cell divergence traverses the
  cell-face CSR and adds the face field for the owner, subtracting it for the
  neighbour, then divides by cell volume.
- `src/cfdx/core/numerics/conservation.h`: the independent face-balance
  reconstruction traverses each stored face once, adds one face flux to its
  owner, subtracts that same value from its neighbour, and accumulates
  boundary flux separately.

The new `OwnerNeighbourTopology` Lean structure models the owner, optional
neighbour, internal-face set, and distinct-owner/neighbour invariant. It also
makes the internal/boundary classification explicit: every internal face has
a neighbour, and every face outside the internal set has no neighbour. Its
incidence coefficient is +1 for the owner and -1 for the neighbour; the
boundary case has no neighbour contribution. The new theorem derives the
zero-sum incidence condition from that representation, rather than taking
the condition as an unexplained axiom.

This is a mathematical correspondence at the level of the owner/neighbour
data model, not yet a proof that the C++ validators and CSR traversal satisfy
the Lean structure for every runtime mesh. An existing small-mesh numerical counterpart is
`tests/unit/test_conservation_assembly.cpp`: its two-cell/one-internal-face
case checks owner +flux, neighbour -flux, and a zero global sum; it also tests
mismatched independently assembled contributions and the boundary-neighbour
zero contribution. This is valuable implementation evidence, but it is not a
machine-checked cross-language equivalence test. The remaining bridge must
establish (and document the exact assumptions of) the mapping from validated
C++ mesh storage to the finite Lean topology, including exactly-once cell-face
references, boundary classification, and shared face-flux storage.

## Gauss-gradient contract: constant preservation

`CFDX.FiniteVolume.GaussGradient` defines a finite-volume Gauss-gradient
operator matching the algebraic structure of
`compute_gradient_gauss` in `src/cfdx/core/numerics/gradient.h`:

[
(G_h \phi)_c =
\frac{1}{V_c}
\sum_f \phi_f\,\mathbf S_{cf},
\qquad
\phi_f =
\begin{cases}
\phi_o, & f\text{ boundary},\\
(1-w_f)\phi_o+w_f\phi_n, & f\text{ internal}.
\end{cases}
]

Here \(\mathbf S_{cf}\) is the face area vector oriented outward from cell
\(c\), and \(V_c>0\). The formal geometry requires the cell closure identity
\(\sum_f\mathbf S_{cf}=0\), with non-incident faces represented by zero
vectors. It also records the production interpolation range
\(0\leq w_f\leq1\).

The kernel-checked theorem
`gauss_gradient_const_eq_zero` proves that a constant cell field has zero
Gauss gradient under this closure identity. The proof is deliberately
conditional: it does not yet prove that every geometry cache built by CFDX
satisfies the closure identity or that every CSR cell-face traversal contains
exactly the intended oriented faces. That is the subject of the subsequent
C++/Lean conformance bridge.

**Important limitation:** this theorem is not affine exactness. CFDX computes
the internal interpolation weight by a clamped projection of the face centre
onto the owner-neighbour centre line. On a skew mesh, that interpolation need
not reproduce the exact affine field value at the face centre. Therefore the
current theorem does not justify a claim of second-order or linear-exact
Gauss gradients on arbitrary unstructured meshes. The next mathematical
extension must state the face-value exactness and geometric first-moment
conditions needed for affine exactness, and quantify the interpolation defect
when those conditions fail.

A separate sufficient-condition theorem,
`affine_gauss_numerator_exact`, states the affine-exactness contract for
exact face-centre values. It assumes (i) geometric area closure and (ii) the
first-moment identity
\[
\sum_f (\mathbf g\cdot\mathbf C_f) S_{cf,d}=V_c g_d.
\]
Under those hypotheses it proves that the Gauss numerator is
\(V_c g_d\). This is deliberately a conditional theorem: it does not yet
prove the moment identity for every CFDX mesh, nor that the current
`GAUSS_TWO_POINT` interpolation supplies exact affine face values. The
condition is relevant to point-reconstructed Gauss and other schemes only
after their reconstruction and geometric moment hypotheses have been
verified against the implementation.

The formal module also includes a kernel-checked counterexample,
`two_point_gauss_skew_face_defect`, matching the C++ test
`two_point_gauss_face_interpolation_has_affine_skewness_defect`. For
\(\phi(x,y,z)=2x+3y+1\), owner centre \((0,0,0)\), neighbour centre
\((1,0,0)\), and skew face centre \((1/2,1/2,0)\), the projection weight is
\(w=1/2\). Owner/neighbour interpolation gives \(2\), while the exact
affine face value is \(3.5\), hence the interpolation defect is \(-1.5\).
A companion theorem confirms exactness for the aligned face centre
\((1/2,0,0)\). This is a concrete counterexample, not a general error bound;
it demonstrates why two-point interpolation cannot be claimed affine-exact
on arbitrary skew faces.

The formal operator is an algebraic specification of the production formula,
not a proof of C++ implementation equivalence. The relevant implementation
path is the interpolation loop in `compute_gradient_gauss`, followed by the
cell-face CSR accumulation with owner-positive and neighbour-negative
`face_Sf`, and division by the positive cell volume.

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
