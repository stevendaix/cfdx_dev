import Mathlib
import CFDX.FiniteVolume.Conservation
import CFDX.FiniteVolume.TwoPointLaplacian

/-!
# Discrete integration by parts for the cell divergence

The divergence implemented by CFDX on face fluxes is the cell balance
\`cellFluxBalance\` (owner +1, internal neighbour -1). This module proves
the exact discrete integration-by-parts identity that the
\`OwnerNeighbourTopology\` incidence justifies, with the boundary term
explicit:

- \`cellFluxBalance_pairing\` : the divergence pairing equals the face sum
  of the flux times the incidence-weighted cell field (owner difference on
  internal faces, owner value on boundary faces) -- the adjointness of the
  incidence matrix and its transpose for the plain Euclidean pairings;
- \`cellFluxBalance_pairing_split\` : the same identity with the boundary
  term isolated;
- \`cellFluxBalance_pairing_internal\` : for a face field vanishing on the
  boundary, the pure adjoint relation with the owner-to-neighbour
  differences.

Faithfulness note (per the Stage 3 audit): the production Gauss gradient
does NOT use the incidence transpose -- it uses face interpolation of cell
gradients. This identity therefore does not prove adjointness of the
production gradient/divergence pair; it isolates exactly the boundary term
and the internal pairing that any such adjointness statement must extend.
-/

namespace CFDX.FiniteVolume

variable {Cell Face : Type*} [Fintype Cell] [Fintype Face]

/-- Boundary version: pairing a field with the incidence of a boundary face
returns the owner value only (a boundary face has no neighbour). -/
theorem sum_weighted_incidence_boundary
    (topology : OwnerNeighbourTopology Cell Face) (φ : Cell → ℝ)
    {f : Face} (hf : f ∉ topology.internalFaces) :
    ∑ c : Cell, φ c * ownerNeighbourIncidence topology c f
      = φ (topology.owner f) := by
  classical
  have hnone := topology.boundary_face_has_no_neighbour f hf
  simp only [ownerNeighbourIncidence, mul_ite, mul_one, mul_zero]
  rw [Finset.sum_add_distrib]
  simp [hnone, Finset.sum_ite_eq, Finset.mem_univ]

/-- Discrete integration by parts (universal form): the cell-divergence
pairing equals the face sum of the flux times the incidence-weighted cell
field. This is the exact adjointness of the incidence matrix and its
transpose for the plain Euclidean pairings on faces and cells. -/
theorem cellFluxBalance_pairing
    (topology : OwnerNeighbourTopology Cell Face) (q : Face → ℝ)
    (φ : Cell → ℝ) :
    ∑ c : Cell, cellFluxBalance topology q c * φ c
      = ∑ f : Face, q f * (if f ∈ topology.internalFaces then
          φ (topology.owner f) - φ (someNeighbour topology f)
        else φ (topology.owner f)) := by
  classical
  simp only [cellFluxBalance, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun f _ => ?_
  rw [← Finset.sum_mul]
  by_cases hf : f ∈ topology.internalFaces
  · rw [sum_weighted_incidence topology φ hf]
    ring
  · rw [sum_weighted_incidence_boundary topology φ hf]
    ring

/-- Discrete integration by parts with the boundary term explicit: internal
faces pair the flux with the owner-to-neighbour difference of the cell
field; boundary faces pair it with the owner value alone. -/
theorem cellFluxBalance_pairing_split
    (topology : OwnerNeighbourTopology Cell Face) (q : Face → ℝ)
    (φ : Cell → ℝ) :
    ∑ c : Cell, cellFluxBalance topology q c * φ c
      = (∑ f ∈ topology.internalFaces,
            q f * (φ (topology.owner f) - φ (someNeighbour topology f)))
        + (∑ f ∈ topology.internalFacesᶜ,
            q f * φ (topology.owner f)) := by
  classical
  rw [cellFluxBalance_pairing]
  rw [← Finset.sum_add_sum_compl topology.internalFaces
        (fun f => q f * (if f ∈ topology.internalFaces then
            φ (topology.owner f) - φ (someNeighbour topology f)
          else φ (topology.owner f)))]
  have h1 : ∑ f ∈ topology.internalFaces,
      q f * (if f ∈ topology.internalFaces then
          φ (topology.owner f) - φ (someNeighbour topology f)
        else φ (topology.owner f))
    = ∑ f ∈ topology.internalFaces,
        q f * (φ (topology.owner f) - φ (someNeighbour topology f)) :=
    Finset.sum_congr rfl (fun f hf => by rw [if_pos hf])
  have h2 : ∑ f ∈ topology.internalFacesᶜ,
      q f * (if f ∈ topology.internalFaces then
          φ (topology.owner f) - φ (someNeighbour topology f)
        else φ (topology.owner f))
    = ∑ f ∈ topology.internalFacesᶜ, q f * φ (topology.owner f) :=
    Finset.sum_congr rfl (fun f hf => by rw [if_neg hf])
  rw [h1, h2]

/-- For a face field vanishing on the boundary, the divergence pairing is
exactly the flux/difference pairing over internal faces. -/
theorem cellFluxBalance_pairing_internal
    (topology : OwnerNeighbourTopology Cell Face)
    (q : Face → ℝ) (hq : ∀ f ∉ topology.internalFaces, q f = 0)
    (φ : Cell → ℝ) :
    ∑ c : Cell, cellFluxBalance topology q c * φ c
      = ∑ f ∈ topology.internalFaces,
          q f * (φ (topology.owner f) - φ (someNeighbour topology f)) := by
  classical
  rw [cellFluxBalance_pairing_split]
  rw [Finset.sum_eq_zero (fun f hf => by rw [hq f hf, zero_mul])]
  ring

end CFDX.FiniteVolume
