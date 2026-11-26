import Mathlib

/-!
# Finite-volume conservation from oriented face incidence

The theorem in this module is a topological conservation statement. Every
internal face has one shared face flux, and its signed incidence over all
adjacent cells sums to zero. The result does not prescribe how the flux is
computed (convective, diffusive, pressure correction, etc.).

Boundary faces are deliberately excluded: their contribution is a physical
boundary flux and must appear in the global balance rather than cancel.
-/

namespace CFDX.FiniteVolume

/-- An oriented finite-volume mesh represented by signed cell-face incidence.

For an internal face, incidence c f is the signed coefficient with which
that face contributes to cell c's balance. The internal-face balance
condition states that the oriented contributions cancel over all cells.
-/
structure OrientedMesh (Cell Face : Type*) [Fintype Cell] [Fintype Face] where
  incidence : Cell → Face → ℝ
  internalFaces : Finset Face
  internal_face_balance :
    ∀ f ∈ internalFaces, ∑ c : Cell, incidence c f = 0

/-- The contribution of the selected internal faces to one cell's balance. -/
def cellInternalFluxBalance
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (mesh : OrientedMesh Cell Face) (flux : Face → ℝ) (c : Cell) : ℝ :=
  ∑ f ∈ mesh.internalFaces, mesh.incidence c f * flux f

/-- Summing all cell balances cancels every internal face exactly.

The face flux is a single value flux f shared by all cells incident to face
f. This shared-flux property is essential: independently reconstructed
left/right face fluxes would not satisfy this theorem without an additional
equality/conformity hypothesis.
-/
theorem internal_face_fluxes_cancel
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (mesh : OrientedMesh Cell Face) (flux : Face → ℝ) :
    ∑ c : Cell, cellInternalFluxBalance mesh flux c = 0 := by
  classical
  unfold cellInternalFluxBalance
  rw [Finset.sum_comm]
  apply Finset.sum_eq_zero
  intro f hf
  rw [← Finset.sum_mul]
  rw [mesh.internal_face_balance f hf]
  simp

end CFDX.FiniteVolume
