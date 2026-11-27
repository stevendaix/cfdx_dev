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

/-- Integrated vector traction stored once per face, in owner orientation. -/
def meshVectorTraction
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (_topology : OwnerNeighbourTopology Cell Face)
    (traction : Face → Fin 3 → ℝ) : Face → Fin 3 → ℝ :=
  traction

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
  simp_rw [Finset.sum_comm]

/-- Internal face work is the owner-minus-neighbour velocity jump dotted
with the one shared integrated traction. -/
def ownerNeighbourInternalVectorWork
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (traction : Face → Fin 3 → ℝ) (face : Face)
    (neighbour : Cell) : ℝ :=
  dot3 (velocity (topology.owner face)) (traction face) -
    dot3 (velocity neighbour) (traction face)

/-- Boundary face work is retained on its sole owner cell. -/
def ownerBoundaryVectorWork
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (traction : Face → Fin 3 → ℝ) (face : Face) : ℝ :=
  dot3 (velocity (topology.owner face)) (traction face)

end CFDX.FiniteVolume
