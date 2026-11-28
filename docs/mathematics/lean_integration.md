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

The reference-cell geometry is now checked explicitly for a unit
tetrahedron. Lean proves exact closure of its four outward-oriented face-area
vectors and the first-moment tensor identity
\[
\sum_f C_{f,i} S_{f,j}=V\,\delta_{ij},\qquad V=1/6.
\]
This verifies the geometric hypothesis for this canonical simplex, with
rational coordinates and exact arithmetic. It does not establish the same
identity for arbitrary CFDX polyhedra, nor validate the production face
centres, area vectors, cell volumes, or CSR incidence. Those remain the
runtime conformance target.

The geometric identity is also substituted into the generic affine
Gauss theorem in `reference_tetra_affine_gauss_numerator_exact`. For any
affine intercept and gradient, and exact values at the four face centroids,
the reference tetrahedron's Gauss numerator is proved to be exactly
\(V\mathbf g\). This closes the formal implication for the reference cell;
it still does not make a claim about the production interpolation on
arbitrary skew meshes or arbitrary geometry-cache inputs.

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

The next conformance fixture extends the reference simplex to two adjacent
unit cubes sharing one internal face. Lean's
`two_cell_first_moment` proves the first-moment identity separately for
both cells. The two cells use local face slots, so the shared face is
represented as cell 0's high-x face (slot 1) and cell 1's low-x face (slot 0).
The additional theorems `two_cell_shared_face_centre_matches` and
`two_cell_shared_face_area_opposes` prove that these local slots have the
same geometric centre and exactly opposite outward-oriented area vectors.
The matching C++ extension evaluates the production geometry cache and
traverses each cell's CSR face list, applying the owner/neighbour sign
convention used by `compute_gradient_gauss`. It checks all nine tensor
components for both cells. This specifically tests the shared-face
orientation bridge; it remains a structured two-cell fixture, not a proof
for arbitrary polyhedral meshes or malformed CSR.

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

## CSR cell-face storage invariants

The Lean module `CFDX.Geometry.CsrTopology` formalizes the storage contract
of the CSR cell-face structure on top of the owner/neighbour topology:

- a cell never stores the same face twice (`CsrCellFaceStorage.nodup`);
- every face is stored by its owner, and every internal face is also stored
  by its neighbour (`owner_stores`, `neighbour_stores`);
- no third cell ever stores a face (`only_owner_neighbour_stores`).

Derived results: an internal face is stored by exactly two cells
(`internal_face_storing_count`), a boundary face by exactly one
(`boundary_face_storing_count`), and the CSR-stored signed incidence of an
internal face cancels over all cells
(`csr_internal_face_incidence_cancels`), which connects the storage
contract directly to the conservation theorem. This is an abstract storage
contract over finite types, not yet a proof that every runtime mesh imported
by CFDX satisfies it; the C++ conformance tests on deterministic meshes
remain the bridge.

## Skew-tetrahedron generalization

The geometric moment identities were so far checked only on the axis-aligned
reference tetrahedron. The new module `CFDX.Geometry.SkewTetrahedron`
repeats the verification on a genuinely skewed cell with vertices
(0,0,0), (1,0,0), (0,2,0) and (1,1,1), whose edge directions `p1 - p0` and
`p3 - p0` (respectively `p2 - p0` and `p3 - p0`) have nonzero dot
products, kernel-checked by `skew_tetra_edges_p1_p3_not_orthogonal` and
`skew_tetra_edges_p2_p3_not_orthogonal`. Lean proves exactly:

- the edge-determinant volume relation `6 V = 2`, i.e. `V = 1/3`
  (`skew_tetra_volume_edge_determinant`);
- exact closure of the four outward-oriented face-area vectors
  (`skew_tetra_area_closure`);
- the nine first-moment components
  `sum_f C_{f,i} S_{f,j} = V delta_{ij}`
  (`skew_tetra_first_moment`);
- exact affine Gauss numerators for exact face-centre values
  (`skew_tetra_affine_gauss_numerator_exact`).

The C++ counterpart in `tests/unit/test_oriented_numerics.cpp` builds the
same tetrahedron in the production `Mesh`/`GeometryCache` structures and
checks the volume, the four oriented face-area vectors, oriented surface
closure, and all nine first-moment components through the runtime cache.

Together with the two-point interpolation counterexample, this isolates the
source of the affine defect: skewness breaks the clamped two-point face
interpolation, while the Gauss operator itself remains affine-exact on the
skewed cell whenever the supplied face values are the exact affine values at
the face centres. This remains a reference-cell result on one skewed simplex;
arbitrary polyhedral cells and runtime CSR conformance are still open.

## Unit-hexahedron generalization

Beyond simplices, the module `CFDX.Geometry.UnitHexahedron` checks the same
identities on the unit cube [0,1]^3 with its six quadrilateral faces: exact
area closure (`unit_hexa_area_closure`), the nine first-moment components
`sum_f C_{f,i} S_{f,j} = V delta_{ij}` with V = 1
(`unit_hexa_first_moment`), and the affine Gauss numerator for exact
face-centre values (`unit_hexa_affine_gauss_numerator_exact`). This is a
non-simplex reference cell; runtime polyhedral conformance remains open.

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


## Discrete energy transfer: conservation is not stability

The module `CFDX.FiniteVolume.EnergyTransfer` is included in the root
`CFDX.lean` build. It defines a two-cell discrete velocity, two equal-and-opposite force
patterns, and three velocity-force pairings. The kernel-checked results establish:

- a conservative internal force can have positive discrete power (+2);
- another conservative force can have negative discrete power (-2);
- the same conservative force can have zero power when paired with a
  uniform velocity.

The theorem `conservation_does_not_imply_nonpositive_power` is a constructive
counterexample to the implication from zero global force balance to
non-positive discrete power. The theorem
`conservative_power_can_have_either_sign` makes the sign-indefiniteness
explicit.

This is a finite-dimensional logical counterexample, not a discretization
of the CFDX Navier–Stokes operator and not a turbulence closure. Its purpose
is to prevent an invalid inference: pairwise conservation of internal face
forces does not establish kinetic-energy stability.

The Euclidean inner product is used only in this minimal example. A production
energy contract must use CFDX's actual cell volumes, density convention,
time-discrete kinetic-energy definition, pressure/velocity coupling,
interpolation, non-orthogonal corrections and boundary work. It must also
allow physically legitimate backscatter. Therefore no blanket
non-positive-power requirement is imposed on every turbulence closure.

The mesh-level vector-work identity and positive-volume weighting have
since been formalized in `CFDX.FiniteVolume.MeshVectorWork`. A bounded C++
audit fixture now compares the volume-weighted work from the production
transpose-viscous force divergence with the corresponding owner/neighbour
face work on a deterministic mesh containing internal and boundary faces.
The remaining implementation bridge is broader: compare the full accepted
momentum residual, including pressure, convection, Rhie–Chow coupling,
non-orthogonal diffusion and temporal terms. Neither full residual
correspondence nor solver stability is claimed yet.


### Face-jump work identity (two-cell prototype)

The same Lean module now defines `twoCellInternalFaceForce q`: the owner
receives the integrated face force (+q), and the neighbour receives
(-q). The theorem `two_cell_power_equals_face_jump_work` proves exactly
[
  sum_{c\in\{P,N\}} U_c f_c
  = (U_P-U_N)q.
]
This is the discrete face-jump work identity for one internal face. It is
the simplest prototype of the summation-by-parts identity needed for a
finite-volume stress divergence: internal face work depends on the jump of
the cell velocity and the shared face traction, not on the global force sum.

The `volumeWeightedPower` definition was refined to distinguish an
integrated face-force residual from a per-volume force density:
[
  \sum_c V_c U_c (F_c/V_c).
]
For strictly positive cell volumes, the same conservative two-cell example
still gives positive power. This removes the dimensional ambiguity of
multiplying an integrated face force directly by a cell volume.

This remains a two-cell scalar prototype. The next extension must handle
vector traction, density-weighted kinetic energy, nonuniform cell volumes,
boundary-face work, and the actual momentum residual including pressure,
viscous stress, Rhie–Chow and temporal terms. It is not yet a theorem about
the C++ `momentum_face_budget` implementation.


## Vector-valued internal and boundary face work

The module `CFDX.FiniteVolume.VectorFaceWork` generalizes the scalar face
identity to three velocity components and an integrated vector traction.
For one internal face, with owner receiving +T and neighbour receiving -T,
`internal_face_vector_work_equals_velocity_jump` proves

\[
\mathbf U_P\cdot\mathbf T_f-\mathbf U_N\cdot\mathbf T_f
=(\mathbf U_P-\mathbf U_N)\cdot\mathbf T_f.
\]

The theorem `internal_face_vector_scatter_is_conservative` checks the
owner/neighbour cancellation component by component. Concrete examples show
that the same scatter convention admits positive or negative work depending
on the velocity orientation. Zero velocity gives zero work. This is
intentional: conservation does not select the sign of energy transfer.

Boundary work is represented separately as
\[
\mathbf U_P\cdot\mathbf T_{f,b},
\]
where the traction orientation carries the boundary sign convention. It is
not paired with a fictitious neighbour velocity.

**Scope:** these are face-local algebraic identities on a three-component
vector. They do not yet sum over a mesh, define outward-normal conventions
for every boundary patch, include density/time weights, or establish that
CFDX's C++ stress and momentum assembly implements the same traction. In
particular, the current C++ transpose-viscous traction must not be relabelled
as a Reynolds-stress closure. The next step is a mesh-level incidence/scatter
contract and a deterministic vector-valued C++ conformance test.


### Finite collection of vector internal faces

The vector module also introduces `InternalFace3` and
`sum_internal_face_work_equals_sum_velocity_jump_work`. The latter lifts
the single-face identity to an arbitrary finite list of internal-face
records by summing the face-local equalities. This is a useful intermediate
lemma for a mesh proof, but the list abstraction does not encode mesh
incidence, owner/neighbour uniqueness, boundary classification, or whether
each geometric face is included exactly once. Those remain explicit proof
obligations before calling this a mesh-global finite-volume theorem.


## Mesh-level vector face-work specification

The module `CFDX.FiniteVolume.MeshVectorWork` now builds on the existing
`OwnerNeighbourTopology` and `ownerNeighbourIncidence`; it does not define
a competing mesh representation. It defines one shared integrated vector
traction per face, the owner/neighbour incidence-scattered cell force, and
the corresponding cell-wise velocity/force work. The theorem
`mesh_vector_work_eq_sum_incidence_face_work` exchanges the finite cell
and face sums, expressing total cell work as a sum of per-face incidence work.

The module now proves that each per-face incidence sum reduces to the
owner-minus-neighbour velocity-jump work for internal faces and to owner-only
work for boundary faces. The global theorem combines these reductions with
the finite cell/face sum interchange. This is an abstract finite-topology
contract, not yet a proof of runtime C++ mesh conformance.

The existing C++ transpose-viscous traction path remains a candidate for a
bounded implementation bridge, but it is a viscous stress operator, not a
Reynolds-stress closure. A conformance test must exercise the production
face scatter and compare integrated traction/work, including boundary faces,
without duplicating the implementation formula in the test.


### Owner/neighbour vector-work reduction

The follow-up module now adds the reduction lemmas
`internal_face_incidence_work_eq_velocity_jump` and
`boundary_face_incidence_work_eq_owner_work`. These use the actual
`OwnerNeighbourTopology` incidence definition: owner coefficient +1,
neighbour coefficient -1, and no neighbour coefficient on boundary faces.
The derived `incidence_face_work_eq_owner_neighbour_work` selects the
neighbour only for faces classified as internal, and the global theorem
`mesh_vector_work_eq_owner_neighbour_face_work` combines this reduction
with the finite cell/face sum interchange.

This is still an abstract algebraic theorem over finite types. It does not
prove runtime mesh validity, face orientation of geometry vectors, density
or time weighting, nor correspondence to the production C++ assembly. The
exact-head CI must compile these new lemmas before they are treated as
machine-checked results.


### Volume-weighted vector power

The mesh vector-work contract now also defines
`cellVectorFaceForceDensity` as integrated face force divided by cell volume,
and `meshVolumeWeightedVectorWork` as the cell-volume-weighted velocity
pairing. The theorem
`mesh_volume_weighted_work_eq_integrated_face_work` states that, for
strictly positive cell volumes, this power pairing is exactly the integrated
face-work sum. This keeps the dimensions distinct: traction integrated over
a face is a force (N), its cell divergence is a force density (N/m^3), and
multiplication by cell volume recovers force before taking velocity work (W).

This theorem does not insert density automatically. The correct density
factor depends on whether the transported momentum equation and closure
traction are expressed in dynamic or kinematic form. That convention must be
fixed explicitly when the Reynolds-stress traction interface is designed.


## Backward-Euler kinetic-energy work identity

The new module `CFDX.FiniteVolume.BackwardEulerEnergy` formalizes the
physical backward-Euler momentum increment paired with the new velocity.
For density \\(\\rho>0\\), time step \\(\\Delta t>0\\), and
\\(\\delta\\mathbf U=\\mathbf U^{n+1}-\\mathbf U^n\\), it proves

\\[
\\rho\\,\\mathbf U^{n+1}\\cdot
\\frac{\\mathbf U^{n+1}-\\mathbf U^n}{\\Delta t}
=
\\frac{\\rho}{\\Delta t}
\\left[
K(\\mathbf U^{n+1})-K(\\mathbf U^n)
+\\frac12\\lVert\\delta\\mathbf U\\rVert^2
\\right],
\\qquad K(\\mathbf U)=\\frac12\\lVert\\mathbf U\\rVert^2.
\\]

The quadratic defect is non-negative for positive density and time step.
This identifies the implicit time-discretization contribution to the
kinetic-energy balance. It does **not** prove stability of the complete
pressure-velocity algorithm, nor does it justify treating pseudo-time,
under-relaxation, or nonlinear globalization terms as physical time work.
The next bridge is to map the exact accepted transient accumulation in
`steady_incompressible_solver.h` to this identity, including variable
cell volumes, density, BDF2 history and boundary/source work.


### Constant-step BDF2 G-stability identity

The same module now adds `constant_step_bdf2_g_stability_identity`. For
constant \\(\\Delta t\\), it expresses the BDF2 time derivative paired
with \\(\\mathbf U^{n+1}\\) as the difference of a two-level G-energy
plus the non-negative second-difference defect. The G-energy is the sum
\\[
\\lVert\\mathbf U^{n+1}\\rVert^2+
\\lVert 2\\mathbf U^{n+1}-\\mathbf U^n\\rVert^2.
\\]
The defect is proportional to
\\(
\\lVert\\mathbf U^{n+1}-2\\mathbf U^n+\\mathbf U^{n-1}\\rVert^2
\\).

**Critical limitation:** CFDX's production BDF2 path calls
`bdf2_coefficients(dt_previous, dt, ...)` and therefore supports
variable-step coefficients. The new theorem is only the constant-step
identity; it must not be used to claim G-stability for arbitrary step ratios.
A variable-step energy estimate with its precise admissible step-ratio
condition is a separate open task. The first-step backward-Euler fallback
and constant-step BDF2 are now covered by separate abstract contracts;
neither yet proves the complete C++ temporal assembly.


#### Variable-step BDF2: prior-art correction and exact scope

The literature audit rejects the simplistic statement that variable-step BDF2
is either unconditionally energy-stable or requires one universal ratio cap.
The answer depends on the stability notion and operator class:

- Akrivis, Chen, Han, Yu & Zhang (2024), *The variable two-step BDF method
  for parabolic equations*, BIT Numerical Mathematics 64:14,
  DOI [10.1007/s10543-024-01007-y](https://doi.org/10.1007/s10543-024-01007-y),
  prove a parabolic stability estimate for a self-adjoint positive elliptic
  operator under adjacent step ratios \\(r_n\\le 1.9398\\), using a
  two-multiplier energy technique. The estimate includes a mesh-variation
  factor; it is not a universal CFL-like limit for every BDF2 application.
- Li & Liao (2022), *Stability of Variable-Step BDF2 and BDF3 Methods*,
  SIAM Journal on Numerical Analysis 60(4), 2253–2272,
  DOI [10.1137/21M1462398](https://doi.org/10.1137/21M1462398),
  establish BDF2 stability on arbitrary time grids using discrete
  orthogonal convolution kernels and an elliptic-type matrix norm. This is
  a different global stability framework, not the same local G-energy
  identity as the constant-step formula above.

Therefore the research target is **not** to invent a blanket step-ratio cap.
It is to specify which estimate CFDX needs (local G-energy telescoping,
global DOC-kernel stability, or a solver-specific energy bound), then check
the actual variable-step coefficients and operator hypotheses. Neither
published result automatically transfers to CFDX's coupled, nonlinear,
convection-dominated, pressure-constrained momentum operator. The constant-
step Lean theorem remains a verified algebraic identity only.
