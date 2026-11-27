import Mathlib
import CFDX.FiniteVolume.Conservation

/-!
# Two-point Laplacian as a weighted graph operator

CFDX's `compute_laplacian` with the `ORTHOGONAL`/`UNCORRECTED` schemes
assembles, on each internal face, the conductance
\$$\\alpha_f = (\\mathbf S_f \\cdot \\mathbf d_f)/|\\mathbf d_f|^2$$\$
and the two-point difference \$$\\phi_N - \\phi_P$$\$. After volume weighting
and a sign change, the resulting matrix is a weighted graph Laplacian
\$$L = B^\\top W B$$\$ on the owner/neighbour incidence.

This module makes that identification explicit on the abstract
`OwnerNeighbourTopology` and proves the structural contracts that the
linear-solver layer may rely on, each under explicit hypotheses:

- the bilinear pairing is symmetric (`laplacian_pairing_symmetric`);
- the quadratic form is \$$\\sum_f \\alpha_f (\\phi_P - \\phi_N)^2$$\$
  (`laplacian_quadratic_form`);
- positive semidefiniteness under nonnegative conductances
  (`laplacian_psd`);
- constants are in the kernel (`twoPointLaplacian_const`);
- under positive conductances and an explicit connectedness assumption on
  the internal-face adjacency, every kernel element is constant
  (`laplacian_kernel_eq_const`), so an anchor in each component makes the
  operator definite (`laplacian_definite_with_anchor`).

Non-goals (faithful to the implementation audit): the corrected/limited
non-orthogonal variants are not covered; no adjointness claim between the
Gauss gradient and this operator is made; no convergence claim is made.
-/

namespace CFDX.FiniteVolume

variable {Cell Face : Type*} [Fintype Cell] [Fintype Face]

/-- Neighbour cell of a face, chosen classically; on an internal face this
is the actual neighbour, on a boundary face it degenerates to the owner
(the difference term then vanishes in every theorem that ranges over
internal faces). -/
noncomputable def someNeighbour
    (topology : OwnerNeighbourTopology Cell Face) (f : Face) : Cell :=
  if h : ∃ n, topology.neighbour f = some n then h.choose else topology.owner f

/-- On an internal face, `someNeighbour` is the actual neighbour. -/
theorem someNeighbour_eq
    (topology : OwnerNeighbourTopology Cell Face)
    {f : Face} (hf : f ∈ topology.internalFaces) :
    topology.neighbour f = some (someNeighbour topology f) := by
  have hex : ∃ n, topology.neighbour f = some n :=
    topology.internal_face_has_neighbour f hf
  unfold someNeighbour
  rw [dif_pos hex]
  exact hex.choose_spec

/-- The CFDX two-point Laplacian on internal faces: for cell \$$c$$\$, the
conductance-weighted sum of signed incidences times the owner-to-neighbour
difference of \$$\\phi$$\$. -/
def twoPointLaplacian (topology : OwnerNeighbourTopology Cell Face)
    (α : Face → ℝ) (φ : Cell → ℝ) (c : Cell) : ℝ :=
  ∑ f ∈ topology.internalFaces,
    ownerNeighbourIncidence topology c f * α f *
      (φ (topology.owner f) - φ (someNeighbour topology f))

/-- Pairing a field with the incidence of one internal face returns the
owner-to-neighbour difference of the field. -/
theorem sum_weighted_incidence
    (topology : OwnerNeighbourTopology Cell Face) (φ : Cell → ℝ)
    {f : Face} (hf : f ∈ topology.internalFaces) :
    ∑ c : Cell, φ c * ownerNeighbourIncidence topology c f
      = φ (topology.owner f) - φ (someNeighbour topology f) := by
  classical
  have hn := someNeighbour_eq topology hf
  simp only [ownerNeighbourIncidence, mul_ite, mul_one, mul_zero, hn,
    Option.some.injEq]
  rw [Finset.sum_add_distrib,
    Finset.sum_ite_eq Finset.univ (topology.owner f) φ,
    Finset.sum_ite_eq Finset.univ (someNeighbour topology f) (fun c => -φ c)]
  simp [Finset.mem_univ]

/-- The pairing \$$\\langle \\varphi, L \\psi \\rangle$$\$ reduces to the
edge sum \$$\\sum_f \\alpha_f (\\varphi_P - \\varphi_N)(\\psi_P - \\psi_N)$$\$. -/
theorem laplacian_bilinear
    (topology : OwnerNeighbourTopology Cell Face) (α : Face → ℝ)
    (φ ψ : Cell → ℝ) :
    ∑ c : Cell, φ c * twoPointLaplacian topology α ψ c
      = ∑ f ∈ topology.internalFaces,
          α f * (φ (topology.owner f) - φ (someNeighbour topology f)) *
            (ψ (topology.owner f) - ψ (someNeighbour topology f)) := by
  classical
  simp only [twoPointLaplacian, Finset.mul_sum, mul_assoc]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun f hf => ?_
  rw [← Finset.sum_mul, sum_weighted_incidence topology φ hf]
  ring

/-- The two-point Laplacian pairing is symmetric. -/
theorem laplacian_pairing_symmetric
    (topology : OwnerNeighbourTopology Cell Face) (α : Face → ℝ)
    (φ ψ : Cell → ℝ) :
    ∑ c : Cell, φ c * twoPointLaplacian topology α ψ c
      = ∑ c : Cell, ψ c * twoPointLaplacian topology α φ c := by
  rw [laplacian_bilinear, laplacian_bilinear]
  refine Finset.sum_congr rfl fun f _ => ?_
  ring

/-- The quadratic form is the conductance-weighted squared difference
across each internal face. -/
theorem laplacian_quadratic_form
    (topology : OwnerNeighbourTopology Cell Face) (α : Face → ℝ)
    (φ : Cell → ℝ) :
    ∑ c : Cell, φ c * twoPointLaplacian topology α φ c
      = ∑ f ∈ topology.internalFaces,
          α f * (φ (topology.owner f) - φ (someNeighbour topology f))^2 := by
  rw [laplacian_bilinear]
  refine Finset.sum_congr rfl fun f _ => ?_
  ring

/-- Positive semidefiniteness under nonnegative face conductances. -/
theorem laplacian_psd
    (topology : OwnerNeighbourTopology Cell Face) (α : Face → ℝ)
    (hα : ∀ f ∈ topology.internalFaces, 0 ≤ α f) (φ : Cell → ℝ) :
    0 ≤ ∑ c : Cell, φ c * twoPointLaplacian topology α φ c := by
  rw [laplacian_quadratic_form]
  exact Finset.sum_nonneg
    (fun f hf => mul_nonneg (hα f hf) (sq_nonneg _))

/-- Constant fields are in the kernel. -/
theorem twoPointLaplacian_const
    (topology : OwnerNeighbourTopology Cell Face) (α : Face → ℝ)
    (μ : ℝ) (c : Cell) :
    twoPointLaplacian topology α (fun _ => μ) c = 0 := by
  simp only [twoPointLaplacian]
  exact Finset.sum_eq_zero (fun f _ => by simp)

/-- Adjacency of two cells across one internal face. -/
def cellAdjacent (topology : OwnerNeighbourTopology Cell Face)
    (u v : Cell) : Prop :=
  ∃ f ∈ topology.internalFaces,
    (topology.owner f = u ∧ topology.neighbour f = some v) ∨
    (topology.owner f = v ∧ topology.neighbour f = some u)

/-- Under strictly positive conductances and connectedness of the internal
face adjacency, every element of the two-point Laplacian kernel is constant.
This is the exact hypothesis under which a pressure/scalar reference (an
anchor per connected component) makes the operator definite. -/
theorem laplacian_kernel_eq_const
    (topology : OwnerNeighbourTopology Cell Face) (α : Face → ℝ)
    (hα : ∀ f ∈ topology.internalFaces, 0 < α f)
    (hconn : ∀ c₁ c₂ : Cell, EqvGen (cellAdjacent topology) c₁ c₂)
    (φ : Cell → ℝ) (hzero : ∀ c, twoPointLaplacian topology α φ c = 0) :
    ∀ c₁ c₂ : Cell, φ c₁ = φ c₂ := by
  classical
  have hq : ∑ c : Cell, φ c * twoPointLaplacian topology α φ c = 0 := by
    refine Finset.sum_eq_zero fun c _ => ?_
    rw [hzero c, mul_zero]
  rw [laplacian_quadratic_form] at hq
  have hface_eq : ∀ u v : Cell, cellAdjacent topology u v → φ u = φ v := by
    intro u v hadj
    obtain ⟨f, hf, hcase⟩ := hadj
    have hdiff : φ (topology.owner f) = φ (someNeighbour topology f) := by
      have hnonneg : ∀ g ∈ topology.internalFaces,
          (0:ℝ) ≤ α g *
            (φ (topology.owner g) - φ (someNeighbour topology g))^2 :=
        fun g hg => mul_nonneg (le_of_lt (hα g hg)) (sq_nonneg _)
      have hterm :=
        (Finset.sum_eq_zero_iff_of_nonneg hnonneg).mp hq f hf
      have hd2 :
          (φ (topology.owner f) - φ (someNeighbour topology f))^2 = 0 := by
        rcases mul_eq_zero.mp hterm with h0 | hs
        · exact absurd h0 (ne_of_gt (hα f hf))
        · exact hs
      linarith [sub_eq_zero.mp (sq_eq_zero_iff.mp hd2)]
    rcases hcase with ⟨hou, hneigh⟩ | ⟨hov, hneigh⟩
    · rw [hou] at hdiff
      have hns := someNeighbour_eq topology hf
      rw [hneigh] at hns
      have hv : v = someNeighbour topology f := Option.some.inj hns
      rw [hv]
      exact hdiff
    · rw [hov] at hdiff
      have hns := someNeighbour_eq topology hf
      rw [hneigh] at hns
      have hu : u = someNeighbour topology f := Option.some.inj hns
      rw [hu]
      exact hdiff.symm
  intro c₁ c₂
  induction hconn with
  | rel u v h => exact hface_eq u v h
  | refl x => rfl
  | symm y x _ ih => exact ih.symm
  | trans x y z _ ih₁ ih₂ => exact ih₁.trans ih₂

/-- With one anchored cell per connected component, the kernel is trivial:
the two-point operator is definite exactly under this anchoring hypothesis. -/
theorem laplacian_definite_with_anchor
    (topology : OwnerNeighbourTopology Cell Face) (α : Face → ℝ)
    (hα : ∀ f ∈ topology.internalFaces, 0 < α f)
    (hconn : ∀ c₁ c₂ : Cell, EqvGen (cellAdjacent topology) c₁ c₂)
    (φ : Cell → ℝ) (hzero : ∀ c, twoPointLaplacian topology α φ c = 0)
    (c₀ : Cell) (hanchor : φ c₀ = 0) (c : Cell) : φ c = 0 := by
  have heq := laplacian_kernel_eq_const topology α hα hconn φ hzero c c₀
  rw [heq]
  exact hanchor

end CFDX.FiniteVolume
