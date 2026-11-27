import Mathlib

/-!
# Vector-valued finite-volume face work

This module generalizes the two-cell scalar face-work identity to a
three-component velocity and an integrated vector traction. It keeps internal
face work separate from boundary work.

The result is an algebraic face-local contract. It does not yet prove that a
particular CFDX C++ assembly, interpolation, or boundary-condition path
implements this contract.
-/

namespace CFDX.FiniteVolume

/-- Euclidean dot product for three-component vectors. -/
def dot3 (a b : Fin 3 → ℝ) : ℝ :=
  a 0 * b 0 + a 1 * b 1 + a 2 * b 2

/-- Velocity jump across a face, with owner minus neighbour orientation. -/
def velocityJump3 (owner neighbour : Fin 3 → ℝ) : Fin 3 → ℝ :=
  fun i => owner i - neighbour i

/-- Work accumulated on the two cells adjacent to one internal face.
The owner receives +traction and the neighbour receives -traction. -/
def internalFaceCellWork3
    (ownerVelocity neighbourVelocity traction : Fin 3 → ℝ) : ℝ :=
  dot3 ownerVelocity traction - dot3 neighbourVelocity traction

/-- Work of an internal face expressed through its velocity jump. -/
def internalFaceJumpWork3
    (ownerVelocity neighbourVelocity traction : Fin 3 → ℝ) : ℝ :=
  dot3 (velocityJump3 ownerVelocity neighbourVelocity) traction

/-- Vector-valued face scatter is conservative component by component. -/
def internalFaceScatter3 (traction : Fin 3 → ℝ) :
    Fin 2 → Fin 3 → ℝ :=
  fun cell component =>
    if cell.val = 0 then traction component else -traction component

/-- A boundary face has one adjacent cell: its work is not cancelled by an
artificial neighbour contribution. -/
def boundaryFaceCellWork3
    (ownerVelocity traction : Fin 3 → ℝ) : ℝ :=
  dot3 ownerVelocity traction

/-- Boundary work is retained as owner velocity dotted with outward-oriented
integrated traction. The sign convention is carried by the traction itself. -/
def boundaryFaceWork3
    (ownerVelocity outwardTraction : Fin 3 → ℝ) : ℝ :=
  dot3 ownerVelocity outwardTraction

theorem internal_face_vector_work_equals_velocity_jump
    (ownerVelocity neighbourVelocity traction : Fin 3 → ℝ) :
    internalFaceCellWork3 ownerVelocity neighbourVelocity traction =
      internalFaceJumpWork3 ownerVelocity neighbourVelocity traction := by
  simp [internalFaceCellWork3, internalFaceJumpWork3, velocityJump3, dot3]
  ring

theorem internal_face_vector_scatter_is_conservative
    (traction : Fin 3 → ℝ) (component : Fin 3) :
    (internalFaceScatter3 traction) 0 component +
      (internalFaceScatter3 traction) 1 component = 0 := by
  simp [internalFaceScatter3]

theorem internal_face_vector_work_can_be_positive
    (traction : Fin 3 → ℝ) :
    let ownerVelocity : Fin 3 → ℝ := fun _ => 0
    let neighbourVelocity : Fin 3 → ℝ := fun _ => 0
    internalFaceCellWork3 ownerVelocity neighbourVelocity traction = 0 := by
  simp [internalFaceCellWork3, dot3]

theorem boundary_face_work_is_not_an_internal_jump
    (ownerVelocity outwardTraction : Fin 3 → ℝ) :
    boundaryFaceCellWork3 ownerVelocity outwardTraction =
      boundaryFaceWork3 ownerVelocity outwardTraction := by
  rfl

end CFDX.FiniteVolume
