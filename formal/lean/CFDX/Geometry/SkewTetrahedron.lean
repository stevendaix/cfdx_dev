import Mathlib
import CFDX.FiniteVolume.GaussGradient

/-!
# Skew-tetrahedron geometric moment identities

The reference tetrahedron in `CFDX.FiniteVolume.GaussGradient` is axis
aligned. This module verifies the same geometric identities on a genuinely
skewed tetrahedron whose edges are not mutually orthogonal.

The tetrahedron has vertices
`p0 = (0,0,0)`, `p1 = (1,0,0)`, `p2 = (0,2,0)`, `p3 = (1,1,1)`.
Face vertices are defined from these points. Their order is outward:
face 0 is (p1,p2,p3), face 1 is (p0,p3,p2), face 2 is (p0,p1,p3),
and face 3 is (p0,p2,p1). In particular, the last ordering is reversed
from the naive (p0,p1,p2) order so its area vector points outward (-z).

Face centres and oriented area vectors below are derived from those vertices.
The explicit coordinate tables are retained as named reference values only
and are proved equal to the derived geometry.

Scope: this is a second mathematical reference cell, not a conformance proof
for CFDX's runtime geometry cache on arbitrary polyhedra.
-/

namespace CFDX.FiniteVolume

/-- Vertex `p0` of the skew reference tetrahedron. -/
def skewTetraP0 : Vec3 := fun _ => 0

/-- Vertex `p1` of the skew reference tetrahedron. -/
def skewTetraP1 : Vec3 := fun d => if d.val = 0 then 1 else 0

/-- Vertex `p2` of the skew reference tetrahedron. -/
def skewTetraP2 : Vec3 := fun d => if d.val = 1 then 2 else 0

/-- Vertex `p3` of the skew reference tetrahedron. -/
def skewTetraP3 : Vec3 := fun _ => 1

/-- Componentwise difference of two three-dimensional vectors. -/
def skewTetraSub (u v : Vec3) : Vec3 := fun d => u d - v d

/-- Three-dimensional cross product, defined componentwise. -/
def skewTetraCross (u v : Vec3) : Vec3 := fun d =>
  if d.val = 0 then u 1 * v 2 - u 2 * v 1
  else if d.val = 1 then u 2 * v 0 - u 0 * v 2
  else u 0 * v 1 - u 1 * v 0

/-- The oriented edge determinant of the declared tetrahedron. -/
def skewTetraEdgeDeterminant : ℝ :=
  dot3
    (skewTetraCross
      (skewTetraSub skewTetraP1 skewTetraP0)
      (skewTetraSub skewTetraP2 skewTetraP0))
    (skewTetraSub skewTetraP3 skewTetraP0)

/-- The four outward-wound triples of vertices, indexed by opposite vertex.
The final face uses (p0,p2,p1), giving the outward -z orientation. -/
def skewTetraFaceVertices (f : Fin 4) : Vec3 × Vec3 × Vec3 :=
  if f.val = 0 then (skewTetraP1, skewTetraP2, skewTetraP3)
  else if f.val = 1 then (skewTetraP0, skewTetraP3, skewTetraP2)
  else if f.val = 2 then (skewTetraP0, skewTetraP1, skewTetraP3)
  else (skewTetraP0, skewTetraP2, skewTetraP1)

/-- Face centroid derived from the three declared face vertices. -/
def skewTetraFaceCentreDerived (f : Fin 4) : Vec3 :=
  let a := (skewTetraFaceVertices f).1
  let b := (skewTetraFaceVertices f).2.1
  let c := (skewTetraFaceVertices f).2.2
  fun d => (a d + b d + c d) / 3

/-- Outward-oriented area vector derived as one half of the cross product
of the two wound triangle edges. -/
def skewTetraOrientedAreaDerived (f : Fin 4) : Vec3 :=
  let a := (skewTetraFaceVertices f).1
  let b := (skewTetraFaceVertices f).2.1
  let c := (skewTetraFaceVertices f).2.2
  fun d =>
    (1 / 2 : ℝ) *
      skewTetraCross (skewTetraSub b a) (skewTetraSub c a) d

/-- The edges p1-p0 and p3-p0 are not orthogonal. -/
theorem skew_tetra_edges_p1_p3_not_orthogonal :
    dot3 (skewTetraSub skewTetraP1 skewTetraP0)
      (skewTetraSub skewTetraP3 skewTetraP0) = 1 := by
  norm_num [dot3, skewTetraSub, skewTetraP0, skewTetraP1,
    skewTetraP3, Fin.sum_univ_succ]

/-- The edges p2-p0 and p3-p0 are not orthogonal. -/
theorem skew_tetra_edges_p2_p3_not_orthogonal :
    dot3 (skewTetraSub skewTetraP2 skewTetraP0)
      (skewTetraSub skewTetraP3 skewTetraP0) = 2 := by
  norm_num [dot3, skewTetraSub, skewTetraP0, skewTetraP2,
    skewTetraP3, Fin.sum_univ_succ]

/-- Explicit reference values for the face centroids. -/
def skewTetraFaceCentre (f : Fin 4) : Vec3 :=
  fun d =>
    if f.val = 0 then
      (if d.val = 0 then 2 / 3 else if d.val = 1 then 1 else 1 / 3)
    else if f.val = 1 then
      (if d.val = 0 then 1 / 3 else if d.val = 1 then 1 else 1 / 3)
    else if f.val = 2 then
      (if d.val = 0 then 2 / 3 else if d.val = 1 then 1 / 3 else 1 / 3)
    else
      (if d.val = 0 then 1 / 3 else if d.val = 1 then 2 / 3 else 0)

/-- Explicit reference values for the outward-oriented face area vectors. -/
def skewTetraOrientedArea (f : Fin 4) : Vec3 :=
  fun d =>
    if f.val = 0 then
      (if d.val = 0 then 1 else if d.val = 1 then 1 / 2 else -(1 / 2))
    else if f.val = 1 then
      (if d.val = 0 then -1 else if d.val = 1 then 0 else 1)
    else if f.val = 2 then
      (if d.val = 0 then 0 else if d.val = 1 then -(1 / 2) else 1 / 2)
    else
      (if d.val = 0 then 0 else if d.val = 1 then 0 else -1)

/-- The explicit centroid table is exactly the centroid derived from vertices. -/
theorem skew_tetra_face_centre_matches_vertices (f : Fin 4) :
    skewTetraFaceCentre f = skewTetraFaceCentreDerived f := by
  funext d
  fin_cases f <;> fin_cases d <;>
    norm_num [skewTetraFaceCentre, skewTetraFaceCentreDerived,
      skewTetraFaceVertices, skewTetraP0, skewTetraP1,
      skewTetraP2, skewTetraP3]

/-- The explicit area table is exactly the outward cross-product area derived
from the declared, outward-wound face vertices. -/
theorem skew_tetra_area_matches_vertices (f : Fin 4) :
    skewTetraOrientedArea f = skewTetraOrientedAreaDerived f := by
  funext d
  fin_cases f <;> fin_cases d <;>
    norm_num [skewTetraOrientedArea, skewTetraOrientedAreaDerived,
      skewTetraFaceVertices, skewTetraSub, skewTetraCross,
      skewTetraP0, skewTetraP1, skewTetraP2, skewTetraP3]

/-- The skew tetrahedron volume is one sixth of the absolute edge determinant. -/
def skewTetraVolume : ℝ := 1 / 3

/-- The edge determinant is computed from the declared vertices, not a
separately hard-coded scalar. -/
theorem skew_tetra_edge_determinant_value :
    skewTetraEdgeDeterminant = 2 := by
  norm_num [skewTetraEdgeDeterminant, dot3, skewTetraCross,
    skewTetraSub, skewTetraP0, skewTetraP1, skewTetraP2,
    skewTetraP3, Fin.sum_univ_succ]

/-- The declared volume agrees with |det(p1-p0,p2-p0,p3-p0)|/6. -/
theorem skew_tetra_volume_from_vertices :
    skewTetraVolume = |skewTetraEdgeDeterminant| / 6 := by
  rw [skew_tetra_edge_determinant_value]
  norm_num [skewTetraVolume]

/-- The determinant formula gives 6 V = 2 for the declared vertices. -/
theorem skew_tetra_volume_edge_determinant :
    6 * skewTetraVolume = 2 := by
  rw [skew_tetra_volume_from_vertices, skew_tetra_edge_determinant_value]
  norm_num

/-- Exact closure of the four outward face-area vectors of the skew tetrahedron. -/
theorem skew_tetra_area_closure (d : Fin 3) :
    (∑ f : Fin 4, skewTetraOrientedArea f d) = 0 := by
  fin_cases d <;>
    norm_num [skewTetraOrientedArea, Fin.sum_univ_succ]

/-- The first geometric moment tensor is V times the identity. -/
theorem skew_tetra_first_moment (i j : Fin 3) :
    (∑ f : Fin 4,
      skewTetraFaceCentre f i * skewTetraOrientedArea f j) =
      if i = j then skewTetraVolume else 0 := by
  fin_cases i <;> fin_cases j <;>
    norm_num [skewTetraFaceCentre, skewTetraOrientedArea,
      skewTetraVolume, Fin.sum_univ_succ]

/-- The generic conditional affine Gauss theorem instantiated on the skew
tetrahedron: for every affine field, exact face-centre values yield the exact
Gauss numerator. -/
theorem skew_tetra_affine_gauss_numerator_exact
    (a : ℝ) (g : Vec3) (d : Fin 3) :
    affineGaussNumeratorComponent
      skewTetraFaceCentre
      (fun _ f => skewTetraOrientedArea f)
      a g () d =
      skewTetraVolume * g d := by
  have hmoment :
      (∑ f : Fin 4,
        dot3 g (skewTetraFaceCentre f) *
          skewTetraOrientedArea f d) =
        skewTetraVolume * g d := by
    fin_cases d <;>
      norm_num [dot3, skewTetraFaceCentre, skewTetraOrientedArea,
        skewTetraVolume, Fin.sum_univ_succ] <;> ring
  exact affine_gauss_numerator_exact
    (faceCentre := skewTetraFaceCentre)
    (orientedArea := fun _ f => skewTetraOrientedArea f)
    (volume := fun _ => skewTetraVolume)
    a g () d
    (skew_tetra_area_closure d)
    hmoment

end CFDX.FiniteVolume
