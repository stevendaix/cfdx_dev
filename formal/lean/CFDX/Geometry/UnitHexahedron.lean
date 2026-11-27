import Mathlib
import CFDX.FiniteVolume.GaussGradient

/-!
# Unit-hexahedron geometric moment identities

The reference cells checked so far are simplices (axis-aligned and skewed
tetrahedra). This module verifies the same geometric identities on the unit
cube [0,1]^3, the simplest polyhedral cell with quadrilateral faces, so that
closure, volume, orientation and the nine first-moment components are
checked on a non-simplex reference cell.

Faces are indexed 0..5: bottom (z=0), top (z=1), front (y=0), back (y=1),
left (x=0), right (x=1). Face centres are quadrilateral centroids and face
area vectors are the outward unit normals (each face has unit area).

Scope: this is a non-simplex mathematical reference cell, not a conformance
proof for CFDX's runtime geometry cache on arbitrary polyhedra.
-/

namespace CFDX.FiniteVolume

/-- Centroid of each face of the unit hexahedron. -/
def unitHexaFaceCentre (f : Fin 6) : Vec3 :=
  fun d =>
    if f.val = 0 then
      (if d.val = 2 then 0 else 1 / 2)
    else if f.val = 1 then
      (if d.val = 2 then 1 else 1 / 2)
    else if f.val = 2 then
      (if d.val = 1 then 0 else 1 / 2)
    else if f.val = 3 then
      (if d.val = 1 then 1 else 1 / 2)
    else if f.val = 4 then
      (if d.val = 0 then 0 else 1 / 2)
    else
      (if d.val = 0 then 1 else 1 / 2)

/-- Outward-oriented area vector of each face of the unit hexahedron. -/
def unitHexaOrientedArea (f : Fin 6) : Vec3 :=
  fun d =>
    if f.val = 0 then
      (if d.val = 2 then -1 else 0)
    else if f.val = 1 then
      (if d.val = 2 then 1 else 0)
    else if f.val = 2 then
      (if d.val = 1 then -1 else 0)
    else if f.val = 3 then
      (if d.val = 1 then 1 else 0)
    else if f.val = 4 then
      (if d.val = 0 then -1 else 0)
    else
      (if d.val = 0 then 1 else 0)

/-- The unit hexahedron has volume 1. -/
def unitHexaVolume : ℝ := 1

/-- Exact closure of the six outward face-area vectors of the unit
hexahedron. -/
theorem unit_hexa_area_closure (d : Fin 3) :
    (∑ f : Fin 6, unitHexaOrientedArea f d) = 0 := by
  fin_cases d <;>
    norm_num [unitHexaOrientedArea, Fin.sum_univ_succ]

/-- The first geometric moment tensor of the unit hexahedron is exactly
V times the identity: sum_f C_f,i * S_f,j = V delta_ij, V = 1. -/
theorem unit_hexa_first_moment (i j : Fin 3) :
    (∑ f : Fin 6,
      unitHexaFaceCentre f i * unitHexaOrientedArea f j) =
      if i = j then unitHexaVolume else 0 := by
  fin_cases i <;> fin_cases j <;>
    norm_num [unitHexaFaceCentre, unitHexaOrientedArea,
      unitHexaVolume, Fin.sum_univ_succ]

/-- The generic conditional affine Gauss theorem instantiated on the unit
hexahedron: for every affine field, exact face-centre values yield the exact
Gauss numerator. -/
theorem unit_hexa_affine_gauss_numerator_exact
    (a : ℝ) (g : Vec3) (d : Fin 3) :
    affineGaussNumeratorComponent
      unitHexaFaceCentre
      (fun _ f => unitHexaOrientedArea f)
      a g () d =
      unitHexaVolume * g d := by
  have hmoment :
      (∑ f : Fin 6,
        dot3 g (unitHexaFaceCentre f) *
          unitHexaOrientedArea f d) =
        unitHexaVolume * g d := by
    fin_cases d <;>
      norm_num [dot3, unitHexaFaceCentre, unitHexaOrientedArea,
        unitHexaVolume, Fin.sum_univ_succ] <;> ring
  exact affine_gauss_numerator_exact
    (faceCentre := unitHexaFaceCentre)
    (orientedArea := fun _ f => unitHexaOrientedArea f)
    (volume := fun _ => unitHexaVolume)
    a g () d
    (unit_hexa_area_closure d)
    hmoment

end CFDX.FiniteVolume
