import Mathlib

/-!
# Gauss-gradient contract matching CFDX's cell-based implementation

The production Gauss gradient interpolates a scalar cell field to faces,
multiplies by the cell-oriented face area vector, sums over each cell's faces,
and divides by cell volume.

This module isolates one exact contract of that implementation:
constant preservation, conditional on geometric closure of each cell's
oriented face-area vectors. It deliberately does not claim affine exactness
for the clamped centre-line interpolation used by CFDX on general skew meshes.
-/

namespace CFDX.FiniteVolume

/-- Three-dimensional vector, represented componentwise over the reals. -/
abbrev Vec3 := Fin 3 → ℝ

/-- Geometry and owner/neighbour interpolation data used by the Gauss gradient.

orientedArea c f is the face area vector oriented outward from cell c;
it is zero when face f is not incident to c. weight f is the geometric
owner-to-neighbour interpolation weight. Boundary faces use the owner value
and do not use the weight.
-/
structure GaussGradientGeometry (Cell Face : Type*)
    [Fintype Cell] [Fintype Face] where
  owner : Face → Cell
  neighbour : Face → Option Cell
  weight : Face → ℝ
  weight_nonneg : ∀ f, 0 ≤ weight f
  weight_le_one : ∀ f, weight f ≤ 1
  orientedArea : Cell → Face → Vec3
  volume : Cell → ℝ
  volume_pos : ∀ c, 0 < volume c
  cell_area_closure : ∀ c, ∑ f : Face, orientedArea c f = 0

/-- CFDX-style interpolation of a cell scalar onto one face.

Internal faces interpolate owner and neighbour values. Boundary faces use
the owner value (zero-gradient boundary treatment in compute_gradient_gauss).
-/
def gaussFaceValue
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (geometry : GaussGradientGeometry Cell Face)
    (cellValue : Cell → ℝ) (f : Face) : ℝ :=
  match geometry.neighbour f with
  | none => cellValue (geometry.owner f)
  | some n =>
      (1 - geometry.weight f) * cellValue (geometry.owner f) +
        geometry.weight f * cellValue n

/-- The face interpolation preserves every constant cell field. -/
theorem gauss_face_value_const
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (geometry : GaussGradientGeometry Cell Face) (a : ℝ) (f : Face) :
    gaussFaceValue geometry (fun _ => a) f = a := by
  unfold gaussFaceValue
  cases h : geometry.neighbour f <;> simp [h] <;> ring

/-- Unscaled cell Gauss-gradient numerator. -/
def gaussGradientNumerator
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (geometry : GaussGradientGeometry Cell Face)
    (cellValue : Cell → ℝ) (c : Cell) : Vec3 :=
  ∑ f : Face,
    gaussFaceValue geometry cellValue f • geometry.orientedArea c f

/-- A constant scalar field has a zero Gauss-gradient numerator on a
geometrically closed cell. -/
theorem gauss_gradient_numerator_const
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (geometry : GaussGradientGeometry Cell Face) (a : ℝ) (c : Cell) :
    gaussGradientNumerator geometry (fun _ => a) c = 0 := by
  unfold gaussGradientNumerator
  calc
    (∑ f : Face,
        gaussFaceValue geometry (fun _ => a) f • geometry.orientedArea c f)
      = ∑ f : Face, a • geometry.orientedArea c f := by
          apply Finset.sum_congr rfl
          intro f hf
          rw [gauss_face_value_const]
    _ = a • ∑ f : Face, geometry.orientedArea c f := by
          rw [← Finset.smul_sum]
    _ = 0 := by simp [geometry.cell_area_closure c]

/-- Volume-normalized cell Gauss gradient, matching the production scaling. -/
def gaussGradient
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (geometry : GaussGradientGeometry Cell Face)
    (cellValue : Cell → ℝ) (c : Cell) : Vec3 :=
  (1 / geometry.volume c) • gaussGradientNumerator geometry cellValue c

/-- Constant-preservation theorem. Under cell geometric closure, the
volume-normalized Gauss gradient annihilates constant fields. -/
theorem gauss_gradient_const_eq_zero
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (geometry : GaussGradientGeometry Cell Face) (a : ℝ) (c : Cell) :
    gaussGradient geometry (fun _ => a) c = 0 := by
  unfold gaussGradient
  rw [gauss_gradient_numerator_const]
  simp


/-! ## Explicit skew-face counterexample for two-point interpolation

The next definitions encode the same affine field and skew face used by
`two_point_gauss_face_interpolation_has_affine_skewness_defect` in the C++
unit tests. The owner-neighbour centre line is the x-axis, while the face
centre is displaced in y. The clamped projection weight is 1/2, so the
two-point interpolation samples the affine field at (1/2, 0, 0), not at the
actual face centre (1/2, 1/2, 0).
-/

/-- The affine scalar field phi(x,y,z) = 2x + 3y + 1. -/
def skewCounterexamplePhi (x y z : ℝ) : ℝ :=
  2 * x + 3 * y + 1

/-- The exact counterexample used by the C++ skewness-defect unit test. -/
theorem two_point_gauss_skew_face_defect :
    (1 - (1 / 2 : ℝ)) * skewCounterexamplePhi 0 0 0 +
        (1 / 2 : ℝ) * skewCounterexamplePhi 1 0 0 -
        skewCounterexamplePhi (1 / 2) (1 / 2) 0 = -(3 / 2 : ℝ) := by
  norm_num [skewCounterexamplePhi]

/-- The same interpolation is exact when the face centre lies on the
owner-neighbour centre line. This isolates skew displacement as the source
of the counterexample above. -/
theorem two_point_gauss_orthogonal_face_exact :
    (1 - (1 / 2 : ℝ)) * skewCounterexamplePhi 0 0 0 +
        (1 / 2 : ℝ) * skewCounterexamplePhi 1 0 0 =
        skewCounterexamplePhi (1 / 2) 0 0 := by
  norm_num [skewCounterexamplePhi]


/-! ## Sufficient moment condition for affine exactness

This contract is intentionally stated separately from the two-point face
interpolation above. It applies when the face values supplied to Gauss are
the exact affine values at the geometric face centres (for example, if a
reconstruction has established that property). It does not assert that the
current two-point interpolation supplies those values on skew meshes.
-/

/-- Euclidean dot product on the three stored spatial components. -/
def dot3 (u v : Vec3) : ℝ :=
  ∑ i : Fin 3, u i * v i

/-- Affine scalar field with intercept a and gradient g. -/
def affineValue3 (a : ℝ) (g x : Vec3) : ℝ :=
  a + dot3 g x

/-- Gauss numerator component using exact affine values at face centres. -/
def affineGaussNumeratorComponent
    {Cell Face : Type*} [Fintype Face]
    (faceCentre : Face → Vec3)
    (orientedArea : Cell → Face → Vec3)
    (a : ℝ) (g : Vec3) (c : Cell) (d : Fin 3) : ℝ :=
  ∑ f : Face,
    affineValue3 a g (faceCentre f) * orientedArea c f d

/-- Affine exactness of the Gauss numerator under geometric closure and
the first-moment identity. The first-moment hypothesis is the componentwise
geometric condition
sum_f (g dot Cf) S_cf = V_c g_d for the requested gradient direction.
It must be verified for the actual face centres, oriented areas and volume.
-/
theorem affine_gauss_numerator_exact
    {Cell Face : Type*} [Fintype Face]
    (faceCentre : Face → Vec3)
    (orientedArea : Cell → Face → Vec3)
    (volume : Cell → ℝ)
    (a : ℝ) (g : Vec3) (c : Cell) (d : Fin 3)
    (area_closure :
      ∑ f : Face, orientedArea c f d = 0)
    (first_moment :
      ∑ f : Face, dot3 g (faceCentre f) * orientedArea c f d =
        volume c * g d) :
    affineGaussNumeratorComponent faceCentre orientedArea a g c d =
      volume c * g d := by
  unfold affineGaussNumeratorComponent affineValue3
  calc
    (∑ f : Face,
        (a + dot3 g (faceCentre f)) * orientedArea c f d)
      = (∑ f : Face, a * orientedArea c f d) +
          ∑ f : Face, dot3 g (faceCentre f) * orientedArea c f d := by
            rw [← Finset.sum_add_distrib]
            apply Finset.sum_congr rfl
            intro f hf
            ring
    _ = 0 + volume c * g d := by
          rw [← Finset.mul_sum, area_closure, first_moment]
    _ = volume c * g d := by ring


/-! ## A concrete tetrahedral geometric moment identity

The conditional affine-exactness theorem above is useful only when its
geometric first-moment hypothesis is established. This section verifies that
identity exactly for the reference tetrahedron with vertices
(0,0,0), (1,0,0), (0,1,0), and (0,0,1). Face area vectors are outward
oriented and face centres are triangle centroids. This is a mathematical
reference-cell result; it is not yet a conformance proof for CFDX's runtime
geometry cache on arbitrary polyhedra.
-/

/-- Centroid of each face of the unit reference tetrahedron. -/
def referenceTetraFaceCentre (f : Fin 4) : Vec3 :=
  fun d =>
    if f.val = 0 then 1 / 3
    else if f.val = 1 then (if d.val = 0 then 0 else 1 / 3)
    else if f.val = 2 then (if d.val = 1 then 0 else 1 / 3)
    else (if d.val = 2 then 0 else 1 / 3)

/-- Outward-oriented area vector of each reference tetrahedron face. -/
def referenceTetraOrientedArea (f : Fin 4) : Vec3 :=
  fun d =>
    if f.val = 0 then 1 / 2
    else if f.val = 1 then (if d.val = 0 then -(1 / 2) else 0)
    else if f.val = 2 then (if d.val = 1 then -(1 / 2) else 0)
    else (if d.val = 2 then -(1 / 2) else 0)

/-- The reference tetrahedron has volume 1/6. -/
def referenceTetraVolume : ℝ := 1 / 6

/-- Exact closure of the four outward face-area vectors. -/
theorem reference_tetra_area_closure (d : Fin 3) :
    (∑ f : Fin 4, referenceTetraOrientedArea f d) = 0 := by
  fin_cases d <;>
    norm_num [referenceTetraOrientedArea, Fin.sum_univ_succ]

/-- The first geometric moment tensor of the reference tetrahedron is
exactly V times the identity: sum_f C_f,i * S_f,j = V δ_ij. -/
theorem reference_tetra_first_moment (i j : Fin 3) :
    (∑ f : Fin 4,
      referenceTetraFaceCentre f i * referenceTetraOrientedArea f j) =
      if i = j then referenceTetraVolume else 0 := by
  fin_cases i <;> fin_cases j <;>
    norm_num [referenceTetraFaceCentre, referenceTetraOrientedArea,
      referenceTetraVolume, Fin.sum_univ_succ]


/-- The generic conditional affine Gauss theorem instantiated on the
reference tetrahedron: its exact face-centre values produce the exact
constant gradient for every affine field. -/
theorem reference_tetra_affine_gauss_numerator_exact
    (a : ℝ) (g : Vec3) (d : Fin 3) :
    affineGaussNumeratorComponent
      referenceTetraFaceCentre
      (fun _ f => referenceTetraOrientedArea f)
      a g () d =
      referenceTetraVolume * g d := by
  have hmoment :
      (∑ f : Fin 4,
        dot3 g (referenceTetraFaceCentre f) *
          referenceTetraOrientedArea f d) =
        referenceTetraVolume * g d := by
    fin_cases d <;>
      norm_num [dot3, referenceTetraFaceCentre,
        referenceTetraOrientedArea, referenceTetraVolume,
        Fin.sum_univ_succ] <;> ring
  exact affine_gauss_numerator_exact
    (faceCentre := referenceTetraFaceCentre)
    (orientedArea := fun _ f => referenceTetraOrientedArea f)
    (volume := fun _ => referenceTetraVolume)
    a g () d
    (reference_tetra_area_closure d)
    hmoment

end CFDX.FiniteVolume
