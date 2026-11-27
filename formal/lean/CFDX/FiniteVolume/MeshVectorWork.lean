import Mathlib
import CFDX.FiniteVolume.Conservation
import CFDX.FiniteVolume.VectorFaceWork

/-!
# Mesh-level vector face work on owner/neighbour topology

This module lifts vector traction work from one face to a finite
owner/neighbour mesh. Internal faces are paired with the owner-minus-neighbour
velocity jump; boundary faces retain owner velocity and are not given a
fictitious neighbour.

The theorem is conditional on the abstract finite topology. It does not
establish that arbitrary runtime C++ mesh storage satisfies this structure.
-/

namespace CFDX.FiniteVolume

/-- Integrated vector force assembled in a cell from the shared face
traction and owner/neighbour incidence. -/
def cellVectorFaceForce
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (traction : Face → Fin 3 → ℝ) (cell : Cell) : Fin 3 → ℝ :=
  fun i => ∑ face : Face,
    ownerNeighbourIncidence topology cell face * traction face i

/-- Global cell-wise work of the assembled face-force field. -/
def meshCellVectorWork
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (traction : Face → Fin 3 → ℝ) : ℝ :=
  ∑ cell : Cell,
    dot3 (velocity cell) (cellVectorFaceForce topology traction cell)

/-- Face work obtained by summing the cell work contributions associated
with that face's signed incidence. -/
def incidenceVectorFaceWork
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (traction : Face → Fin 3 → ℝ) (face : Face) : ℝ :=
  ∑ cell : Cell,
    dot3 (velocity cell)
      (fun i => ownerNeighbourIncidence topology cell face * traction face i)

/-- The mesh work equals the sum of the work of its individual shared faces.
This is the finite-sum interchange underlying the discrete integration by
parts identity; no face is duplicated by the algebraic definition. -/
theorem mesh_vector_work_eq_sum_incidence_face_work
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (traction : Face → Fin 3 → ℝ) :
    meshCellVectorWork topology velocity traction =
      ∑ face : Face, incidenceVectorFaceWork topology velocity traction face := by
  classical
  unfold meshCellVectorWork cellVectorFaceForce incidenceVectorFaceWork
  simp_rw [dot3]
  simp_rw [Finset.sum_add_distrib]
  simp_rw [Finset.mul_sum]
  simp_rw [Finset.sum_comm]

/-- For an internal face, incidence work is exactly the owner-minus-neighbour
velocity jump dotted with the shared traction. -/
theorem internal_face_incidence_work_eq_velocity_jump
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (traction : Face → Fin 3 → ℝ) (face : Face) (neighbour : Cell)
    (hneighbour : topology.neighbour face = some neighbour) :
    incidenceVectorFaceWork topology velocity traction face =
      dot3 (velocity (topology.owner face) - velocity neighbour) (traction face) := by
  classical
  unfold incidenceVectorFaceWork ownerNeighbourIncidence
  simp [dot3, hneighbour, Finset.sum_add_distrib, Finset.sum_ite_eq'] <;> ring

/-- For a boundary face, incidence work is owner work only. -/
theorem boundary_face_incidence_work_eq_owner_work
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (traction : Face → Fin 3 → ℝ) (face : Face)
    (hneighbour : topology.neighbour face = none) :
    incidenceVectorFaceWork topology velocity traction face =
      dot3 (velocity (topology.owner face)) (traction face) := by
  classical
  unfold incidenceVectorFaceWork ownerNeighbourIncidence
  simp [dot3, hneighbour, Finset.sum_add_distrib, Finset.sum_ite_eq']

/-- The canonical work contribution of a face, using the topology's
internal/boundary partition. The neighbour is selected only for internal
faces, whose existence is part of OwnerNeighbourTopology. -/
def ownerNeighbourFaceWork
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (traction : Face → Fin 3 → ℝ) (face : Face) : ℝ :=
  if h : face ∈ topology.internalFaces then
    let neighbour := Classical.choose (topology.internal_face_has_neighbour face h)
    dot3 (velocity (topology.owner face) - velocity neighbour) (traction face)
  else
    dot3 (velocity (topology.owner face)) (traction face)

/-- Every face's incidence work reduces to its physical owner/neighbour
velocity-jump or owner-only boundary work. -/
theorem incidence_face_work_eq_owner_neighbour_work
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (traction : Face → Fin 3 → ℝ) (face : Face) :
    incidenceVectorFaceWork topology velocity traction face =
      ownerNeighbourFaceWork topology velocity traction face := by
  classical
  unfold ownerNeighbourFaceWork
  split_ifs with h
  · let neighbour := Classical.choose (topology.internal_face_has_neighbour face h)
    have hneighbour : topology.neighbour face = some neighbour :=
      Classical.choose_spec (topology.internal_face_has_neighbour face h)
    rw [internal_face_incidence_work_eq_velocity_jump
      topology velocity traction face neighbour hneighbour]
    rfl
  · have hneighbour : topology.neighbour face = none :=
      topology.boundary_face_has_no_neighbour face h
    exact boundary_face_incidence_work_eq_owner_work
      topology velocity traction face hneighbour

/-- Mesh work is the sum of physical face work: internal faces contribute
the owner-minus-neighbour velocity jump, while boundary faces contribute
owner velocity dotted with the outward-oriented traction. -/
theorem mesh_vector_work_eq_owner_neighbour_face_work
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (traction : Face → Fin 3 → ℝ) :
    meshCellVectorWork topology velocity traction =
      ∑ face : Face, ownerNeighbourFaceWork topology velocity traction face := by
  rw [mesh_vector_work_eq_sum_incidence_face_work]
  exact Finset.sum_congr rfl (fun face _ =>
    incidence_face_work_eq_owner_neighbour_work topology velocity traction face)

/-- Internal face work is the owner-minus-neighbour velocity jump dotted
with the one shared integrated traction. -/
def ownerNeighbourInternalVectorWork
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (traction : Face → Fin 3 → ℝ) (face : Face)
    (neighbour : Cell) : ℝ :=
  dot3 (velocity (topology.owner face) - velocity neighbour) (traction face)

/-- Boundary face work is retained on its sole owner cell. -/
def ownerBoundaryVectorWork
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (traction : Face → Fin 3 → ℝ) (face : Face) : ℝ :=
  dot3 (velocity (topology.owner face)) (traction face)

end CFDX.FiniteVolume
