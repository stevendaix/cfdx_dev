import Mathlib

/-!
# Finite-volume conservation from oriented face incidence

The theorem in this module is a topological conservation statement. Every
internal face has one shared face flux, and its signed incidence over all
adjacent cells sums to zero. The result does not prescribe how the flux is
computed (convective, diffusive, pressure correction, etc.).

Boundary faces are deliberately excluded: their contribution is a physical
boundary flux and must appear in the global balance rather than cancel.

This is a topological conservation theorem only. It does not prove that a
chosen flux is consistent with the PDE, bounded, monotone, energy stable,
coercive, or convergent. Those properties depend on the flux reconstruction,
diffusion coefficients, boundary treatment, discrete inner products and
linear/nonlinear solver; they require separate contracts.
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

/-- Owner/neighbour topology matching CFDX's face ownership convention.

A missing neighbour represents a boundary face. Every internal face has a
neighbour distinct from its owner, matching FaceOwnership::is_consistent.
-/
structure OwnerNeighbourTopology (Cell Face : Type*)
    [Fintype Cell] [Fintype Face] where
  owner : Face → Cell
  neighbour : Face → Option Cell
  internalFaces : Finset Face
  internal_face_has_neighbour :
    ∀ f ∈ internalFaces, ∃ n, neighbour f = some n
  internal_face_owner_neighbour_distinct :
    ∀ f ∈ internalFaces, ∀ n, neighbour f = some n → owner f ≠ n
  boundary_face_has_no_neighbour :
    ∀ f, f ∉ internalFaces → neighbour f = none

/-- Signed incidence induced by one owner and an optional neighbour.

The owner contributes +1; an existing neighbour contributes -1. A boundary
face has no neighbour contribution.
-/
def ownerNeighbourIncidence
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (c : Cell) (f : Face) : ℝ :=
  (if topology.owner f = c then 1 else 0) +
  (if topology.neighbour f = some c then -1 else 0)

/-- The signed incidence of every internal owner/neighbour face sums to zero.

This derives the abstract incidence hypothesis from the owner/neighbour
representation rather than assuming it directly.
-/
theorem owner_neighbour_incidence_cancels
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    {f : Face} (hf : f ∈ topology.internalFaces) :
    ∑ c : Cell, ownerNeighbourIncidence topology c f = 0 := by
  classical
  obtain ⟨n, hneighbour⟩ := topology.internal_face_has_neighbour f hf
  rw [Finset.sum_add_distrib]
  simp [ownerNeighbourIncidence, hneighbour]

/-- Convert owner/neighbour topology to the general oriented-mesh contract. -/
def OwnerNeighbourTopology.toOrientedMesh
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face) : OrientedMesh Cell Face where
  incidence := ownerNeighbourIncidence topology
  internalFaces := topology.internalFaces
  internal_face_balance := by
    intro f hf
    exact owner_neighbour_incidence_cancels topology hf

/-- Global internal-face conservation derived from owner/neighbour topology. -/
theorem owner_neighbour_internal_fluxes_cancel
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face) (flux : Face → ℝ) :
    ∑ c : Cell,
      cellInternalFluxBalance topology.toOrientedMesh flux c = 0 := by
  exact internal_face_fluxes_cancel topology.toOrientedMesh flux

/-- The complete face-flux contribution to one cell, including boundaries. -/
def cellFluxBalance
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (flux : Face → ℝ) (c : Cell) : ℝ :=
  ∑ f : Face, ownerNeighbourIncidence topology c f * flux f

/-- The global sum of cell face balances equals the outward boundary flux.

Internal faces cancel pairwise. A boundary face contributes once, through
its owner, because it has no neighbour cell.
-/
theorem global_cell_flux_balance_eq_boundary_flux
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face) (flux : Face → ℝ) :
    ∑ c : Cell, cellFluxBalance topology flux c =
      ∑ f : Face, if f ∉ topology.internalFaces then flux f else 0 := by
  classical
  unfold cellFluxBalance
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro f hf
  rw [← Finset.sum_mul]
  by_cases hInternal : f ∈ topology.internalFaces
  · rw [owner_neighbour_incidence_cancels topology hInternal]
    simp
  · have hBoundary := topology.boundary_face_has_no_neighbour f hInternal
    simp [ownerNeighbourIncidence, hBoundary]

end CFDX.FiniteVolume
