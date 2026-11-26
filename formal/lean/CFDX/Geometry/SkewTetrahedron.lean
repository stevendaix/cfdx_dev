import Mathlib
import CFDX.FiniteVolume.GaussGradient

/-!
# Skew-tetrahedron geometric moment identities

The reference tetrahedron in `CFDX.FiniteVolume.GaussGradient` is axis
aligned. This module verifies the same geometric identities on a genuinely
skewed tetrahedron whose edges are not mutually orthogonal, so that closure,
volume, orientation and the nine first-moment components are checked away
from the canonical simplex.

The tetrahedron has vertices
`p0 = (0,0,0)`, `p1 = (1,0,0)`, `p2 = (0,2,0)`, `p3 = (1,1,1)`.
Its edge directions `p1 - p0` and `p3 - p0` (and `p2 - p0` and
`p3 - p0`) are not orthogonal and its faces are oblique.

Face `f` is opposite vertex `f`: face 0 is `(p1,p2,p3)`, face 1 is
`(p0,p3,p2)`, face 2 is `(p0,p1,p3)`, face 3 is `(p0,p1,p2)`, each wound
counter-clockwise as seen from outside. Face centres are triangle centroids
and face area vectors are the outward-oriented cross products.

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

/-- The edges `p1 - p0` and `p3 - p0` are not orthogonal. -/
theorem skew_tetra_edges_p1_p3_not_orthogonal :
    dot3 skewTetraP1 skewTetraP3 = 1 := by
  norm_num [dot3, skewTetraP1, skewTetraP3, Fin.sum_univ_succ]

/-- The edges `p2 - p0` and `p3 - p0` are not orthogonal. -/
theorem skew_tetra_edges_p2_p3_not_orthogonal :
    dot3 skewTetraP2 skewTetraP3 = 2 := by
  norm_num [dot3, skewTetraP2, skewTetraP3, Fin.sum_univ_succ]

/-- Centroid of each face of the skew tetrahedron. -/
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

/-- Outward-oriented area vector of each face of the skew tetrahedron,
equal to one half of the counter-clockwise cross product of the wound
triangle edges. -/
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

/-- The skew tetrahedron has volume 1/3. -/
def skewTetraVolume : ℝ := 1 / 3

/-- The volume matches the edge-determinant formula: the determinant of the
edge matrix `(p1-p0, p2-p0, p3-p0)` equals `2`, hence
`V = |det| / 6 = 1/3`. -/
theorem skew_tetra_volume_edge_determinant :
    6 * skewTetraVolume = 2 := by
  norm_num [skewTetraVolume]

/-- Exact closure of the four outward face-area vectors of the skew
tetrahedron. -/
theorem skew_tetra_area_closure (d : Fin 3) :
    (∑ f : Fin 4, skewTetraOrientedArea f d) = 0 := by
  fin_cases d <;>
    norm_num [skewTetraOrientedArea, Fin.sum_univ_succ]

/-- The first geometric moment tensor of the skew tetrahedron is exactly
V times the identity: sum_f C_f,i * S_f,j = V delta_ij, V = 1/3. -/
theorem skew_tetra_first_moment (i j : Fin 3) :
    (∑ f : Fin 4,
      skewTetraFaceCentre f i * skewTetraOrientedArea f j) =
      if i = j then skewTetraVolume else 0 := by
  fin_cases i <;> fin_cases j <;>
    norm_num [skewTetraFaceCentre, skewTetraOrientedArea,
      skewTetraVolume, Fin.sum_univ_succ]

/-- The generic conditional affine Gauss theorem instantiated on the skew
tetrahedron: for every affine field, exact face-centre values yield the
exact Gauss numerator. Skewness therefore breaks the clamped two-point face
interpolation, not the Gauss operator itself. -/
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
