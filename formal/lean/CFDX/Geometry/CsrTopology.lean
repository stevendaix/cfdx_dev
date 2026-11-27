import Mathlib
import CFDX.FiniteVolume.Conservation

/-!
# CSR cell-face storage conformance invariants

CFDX stores each cell's incident faces in CSR-like per-cell lists
(`cells()` offsets plus face references) and its mesh validation checks that
every face is referenced by its owner and, for internal faces, by its
neighbour. This module formalizes those storage invariants on top of the
`OwnerNeighbourTopology` of `CFDX.FiniteVolume.Conservation`, and derives:

- an internal face is stored by exactly two cells, its owner and its
  neighbour (`internal_face_storing_count`);
- a boundary face is stored by exactly one cell, its owner
  (`boundary_face_storing_count`);
- the CSR-stored signed incidence of an internal face cancels over all
  cells (`csr_internal_face_incidence_cancels`), connecting the storage
  contract to the conservation theorem;
- a boundary face contributes exactly its owner incidence
  (`csr_boundary_face_incidence_eq_one`).

Scope: this is an abstract storage contract over finite cell/face types. It
is not yet a proof that the C++ `Mesh` CSR construction satisfies it for
every runtime mesh; deterministic small-mesh C++ conformance tests remain the
bridge to the production structures.
-/

namespace CFDX.FiniteVolume

variable {Cell Face : Type*} [Fintype Cell] [Fintype Face]

/-- CSR-style per-cell face storage attached to an owner/neighbour topology.

`cellFaces c` is the face list stored for cell `c`, in storage order. The
fields encode the invariants validated by the C++ mesh layer: a cell never
stores the same face twice, every face is stored by its owner, every internal
face is also stored by its neighbour, and no third cell ever stores a face.
-/
structure CsrCellFaceStorage (topology : OwnerNeighbourTopology Cell Face) where
  /-- Per-cell stored face lists, as in the CSR cells() structure. -/
  cellFaces : Cell → List Face
  /-- A cell never stores the same face twice (exactly-once CSR reference). -/
  nodup : ∀ c, (cellFaces c).NoDup
  /-- Every face is stored by its owner cell. -/
  owner_stores : ∀ f, f ∈ cellFaces (topology.owner f)
  /-- Every internal face is also stored by its neighbour cell. -/
  neighbour_stores :
    ∀ f ∈ topology.internalFaces,
      ∃ n, topology.neighbour f = some n ∧ f ∈ cellFaces n
  /-- No cell besides the owner, or the neighbour of an internal face,
      stores the face. -/
  only_owner_neighbour_stores :
    ∀ (c : Cell) (f : Face), f ∈ cellFaces c →
      c = topology.owner f ∨
        (f ∈ topology.internalFaces ∧
          ∃ n, topology.neighbour f = some n ∧ c = n)

/-- Number of cells whose stored list contains face `f`. -/
def storingCellCount (topology : OwnerNeighbourTopology Cell Face)
    (storage : CsrCellFaceStorage topology) (f : Face) : ℕ :=
  ∑ c : Cell, if f ∈ storage.cellFaces c then 1 else 0

/-- Characterization of the cells storing an internal face: exactly the
owner and the given neighbour. -/
theorem storing_cells_internal_iff
    (topology : OwnerNeighbourTopology Cell Face)
    (storage : CsrCellFaceStorage topology)
    {f : Face} {n : Cell} (hf : f ∈ topology.internalFaces)
    (hn : topology.neighbour f = some n)
    {c : Cell} :
    f ∈ storage.cellFaces c ↔ (c = topology.owner f ∨ c = n) := by
  classical
  constructor
  · intro hmem
    rcases storage.only_owner_neighbour_stores c f hmem with hco | ⟨_, hn', hc⟩
    · exact Or.inl hco
    · obtain ⟨n', hn'', hc'⟩ := hn'
      have hnn : n = n' := Option.some.inj (hn.symm.trans hn'')
      subst hnn
      exact Or.inr hc'
  · rintro (rfl | rfl)
    · exact storage.owner_stores f
    · obtain ⟨n', hn'', hc'⟩ := storage.neighbour_stores f hf
      have hnn : n' = n := Option.some.inj (hn.symm.trans hn'')
      subst hnn
      exact hc'

/-- Characterization of the cells storing a boundary face: only the owner. -/
theorem storing_cells_boundary_iff
    (topology : OwnerNeighbourTopology Cell Face)
    (storage : CsrCellFaceStorage topology)
    {f : Face} (hf : f ∉ topology.internalFaces)
    {c : Cell} :
    f ∈ storage.cellFaces c ↔ c = topology.owner f := by
  classical
  constructor
  · intro hmem
    rcases storage.only_owner_neighbour_stores c f hmem with h | ⟨hint, _, _, _⟩
    · exact h
    · exact absurd hint hf
  · intro h
    rw [h]
    exact storage.owner_stores f

/-- An internal face is stored by exactly two cells: its owner and its
neighbour. -/
theorem internal_face_storing_count
    (topology : OwnerNeighbourTopology Cell Face)
    (storage : CsrCellFaceStorage topology)
    {f : Face} (hf : f ∈ topology.internalFaces) :
    storingCellCount topology storage f = 2 := by
  classical
  obtain ⟨n, hn⟩ := topology.internal_face_has_neighbour f hf
  have hdistinct : topology.owner f ≠ n :=
    fun h => topology.internal_face_owner_neighbour_distinct f hf n hn h
  have key : ∀ c : Cell,
      (if f ∈ storage.cellFaces c then (1 : ℕ) else 0)
        = (if c = topology.owner f then 1 else 0) + (if c = n then 1 else 0) := by
    intro c
    by_cases h : f ∈ storage.cellFaces c
    · rcases (storing_cells_internal_iff topology storage hf hn).mp h with hc | hc
      · have hcn : c ≠ n := fun h2 => hdistinct (hc.symm.trans h2)
        rw [if_pos h, if_pos hc, if_neg hcn, if_pos h]
        norm_num
      · have hco : c ≠ topology.owner f := fun h2 => hdistinct (h2.trans hc.symm)
        rw [if_pos h, if_neg hco, if_pos hc, if_pos h]
        norm_num
    · have hco : c ≠ topology.owner f :=
        fun h2 => h (by rw [← h2]; exact storage.owner_stores f)
      have hcn : c ≠ n := by
        intro h2
        obtain ⟨n', hn', hc'⟩ := storage.neighbour_stores f hf
        have hnn : n' = n := Option.some.inj (hn.symm.trans hn')
        subst hnn
        refine h ?_
        rw [← h2]
        exact hc'
      rw [if_neg h, if_neg hco, if_neg hcn]
      norm_num
  calc storingCellCount topology storage f
      = ∑ c : Cell, (if f ∈ storage.cellFaces c then 1 else 0) := rfl
    _ = ∑ c : Cell,
          ((if c = topology.owner f then 1 else 0) + (if c = n then 1 else 0)) :=
        Finset.sum_congr rfl (fun c _ => key c)
    _ = (∑ c : Cell, (if c = topology.owner f then 1 else 0)) +
        (∑ c : Cell, (if c = n then 1 else 0)) := Finset.sum_add_distrib
    _ = 2 := by
        rw [Finset.sum_ite_eq' Finset.univ (topology.owner f) (fun _ => (1 : ℕ)),
          Finset.sum_ite_eq' Finset.univ n (fun _ => (1 : ℕ))]
        simp [Finset.mem_univ]

/-- A boundary face is stored by exactly one cell: its owner. -/
theorem boundary_face_storing_count
    (topology : OwnerNeighbourTopology Cell Face)
    (storage : CsrCellFaceStorage topology)
    {f : Face} (hf : f ∉ topology.internalFaces) :
    storingCellCount topology storage f = 1 := by
  classical
  have key : ∀ c : Cell,
      (if f ∈ storage.cellFaces c then (1 : ℕ) else 0)
        = (if c = topology.owner f then 1 else 0) := by
    intro c
    by_cases h : f ∈ storage.cellFaces c
    · rw [if_pos h,
        if_pos ((storing_cells_boundary_iff topology storage hf).mp h)]
    · rw [if_neg h,
        if_neg (fun hc => h ((storing_cells_boundary_iff topology storage hf).mpr hc))]
  calc storingCellCount topology storage f
      = ∑ c : Cell, (if f ∈ storage.cellFaces c then 1 else 0) := rfl
    _ = ∑ c : Cell, (if c = topology.owner f then 1 else 0) :=
        Finset.sum_congr rfl (fun c _ => key c)
    _ = 1 := by
        rw [Finset.sum_ite_eq' Finset.univ (topology.owner f) (fun _ => (1 : ℕ))]
        simp [Finset.mem_univ]

/-- For every cell, the CSR-stored incidence coincides with the abstract
owner/neighbour incidence, so the storage contract carries the conservation
framework unchanged. -/
theorem csr_stored_incidence_eq_topology
    (topology : OwnerNeighbourTopology Cell Face)
    (storage : CsrCellFaceStorage topology)
    (f : Face) (c : Cell) :
    (if f ∈ storage.cellFaces c then
      ownerNeighbourIncidence topology c f else (0 : ℝ))
      = ownerNeighbourIncidence topology c f := by
  classical
  by_cases h : f ∈ storage.cellFaces c
  · rw [if_pos h]
  · rw [if_neg h]
    have hco : topology.owner f ≠ c :=
      fun h2 => h (by rw [← h2]; exact storage.owner_stores f)
    have hnc : topology.neighbour f ≠ some c := by
      intro h2
      obtain ⟨n', hn', hc'⟩ := storage.neighbour_stores f hf
      have hnn : c = n' := Option.some.inj (h2.symm.trans hn')
      refine h ?_
      rw [hnn]
      exact hc'
    simp only [ownerNeighbourIncidence]
    rw [if_neg hco, if_neg hnc]
    norm_num

/-- The CSR-stored signed incidence of an internal face cancels over all
cells. This derives the conservation balance from the storage contract
itself. -/
theorem csr_internal_face_incidence_cancels
    (topology : OwnerNeighbourTopology Cell Face)
    (storage : CsrCellFaceStorage topology)
    {f : Face} (hf : f ∈ topology.internalFaces) :
    ∑ c : Cell,
      (if f ∈ storage.cellFaces c then
        ownerNeighbourIncidence topology c f else (0 : ℝ)) = 0 := by
  classical
  calc ∑ c : Cell,
          (if f ∈ storage.cellFaces c then
            ownerNeighbourIncidence topology c f else (0 : ℝ))
      = ∑ c : Cell, ownerNeighbourIncidence topology c f :=
        Finset.sum_congr rfl (fun c _ => csr_stored_incidence_eq_topology topology storage f c)
    _ = 0 := owner_neighbour_incidence_cancels topology hf

/-- A boundary face contributes exactly its owner incidence (+1) through the
CSR storage. -/
theorem csr_boundary_face_incidence_eq_one
    (topology : OwnerNeighbourTopology Cell Face)
    (storage : CsrCellFaceStorage topology)
    {f : Face} (hf : f ∉ topology.internalFaces) :
    ∑ c : Cell,
      (if f ∈ storage.cellFaces c then
        ownerNeighbourIncidence topology c f else (0 : ℝ)) = 1 := by
  classical
  have hnone : topology.neighbour f = none := topology.boundary_face_has_no_neighbour f hf
  have key : ∀ c : Cell,
      (if f ∈ storage.cellFaces c then
        ownerNeighbourIncidence topology c f else (0 : ℝ))
        = (if c = topology.owner f then (1 : ℝ) else 0) := by
    intro c
    by_cases h : f ∈ storage.cellFaces c
    · have hco : c = topology.owner f :=
        (storing_cells_boundary_iff topology storage hf).mp h
      rw [if_pos h, if_pos hco.symm]
      simp only [ownerNeighbourIncidence]
      rw [if_pos rfl, if_neg (by rw [hnone]; simp)]
      norm_num
    · rw [if_neg h,
        if_neg (fun hc => h ((storing_cells_boundary_iff topology storage hf).mpr hc))]
  calc ∑ c : Cell,
          (if f ∈ storage.cellFaces c then
            ownerNeighbourIncidence topology c f else (0 : ℝ))
      = ∑ c : Cell, (if c = topology.owner f then (1 : ℝ) else 0) :=
        Finset.sum_congr rfl (fun c _ => key c)
    _ = 1 := by
        rw [Finset.sum_ite_eq' Finset.univ (topology.owner f) (fun _ => (1 : ℝ))]
        simp [Finset.mem_univ]

end CFDX.FiniteVolume
