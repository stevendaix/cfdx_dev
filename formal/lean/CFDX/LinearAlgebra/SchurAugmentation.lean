import Mathlib

/-!
# Conditional pressure-block row elimination

This module formalizes the exact algebraic row operation used by the
C-aware Schur construction. It assumes an explicit factorization
C = L ∘ G. It does not prove that such an L exists for a given CFDX
Rhie-Chow operator, that L is sparse, or that the transformed momentum
block is symmetric, positive definite, or well conditioned.
-/

namespace CFDX.LinearAlgebra

variable {U V P : Type*}
  [AddCommGroup U] [Module ℝ U]
  [AddCommGroup V] [Module ℝ V]
  [AddCommGroup P] [Module ℝ P]

/-- The continuity operator after eliminating the pressure block by
    subtracting L times the momentum row. -/
def schurEliminatedD
    (A : U →ₗ[ℝ] V) (D : U →ₗ[ℝ] P)
    (L : V →ₗ[ℝ] P) : U →ₗ[ℝ] P :=
  D - L.comp A

/-- If C factors as L ∘ G, subtracting L times the momentum equation
    from the continuity equation removes C exactly and replaces D by
    D - L ∘ A. This is an identity of row operators, not an approximation. -/
theorem pressure_block_row_elimination
    (A : U →ₗ[ℝ] V) (G : P →ₗ[ℝ] V)
    (D : U →ₗ[ℝ] P) (C : P →ₗ[ℝ] P)
    (L : V →ₗ[ℝ] P)
    (hC : C = L.comp G) (u : U) (p : P) :
    D u + C p - L (A u + G p) =
      schurEliminatedD A D L u := by
  rw [hC]
  simp only [LinearMap.comp_apply, map_add, schurEliminatedD,
    LinearMap.sub_apply]
  abel

/-- The complete block residual after the lower-row operation. The
    right-hand side must undergo the same operation: (f,g) becomes
    (f,g-L f). -/
def pressureBlockRowElimination
    (A : U →ₗ[ℝ] V) (G : P →ₗ[ℝ] V)
    (D : U →ₗ[ℝ] P) (C : P →ₗ[ℝ] P)
    (L : V →ₗ[ℝ] P) (u : U) (p : P) : V × P :=
  (A u + G p, D u + C p - L (A u + G p))

/-- The transformed block residual has zero pressure-pressure block
    and the exact eliminated continuity operator. -/
theorem pressure_block_row_elimination_eq
    (A : U →ₗ[ℝ] V) (G : P →ₗ[ℝ] V)
    (D : U →ₗ[ℝ] P) (C : P →ₗ[ℝ] P)
    (L : V →ₗ[ℝ] P)
    (hC : C = L.comp G) (u : U) (p : P) :
    pressureBlockRowElimination A G D C L u p =
      (A u + G p, schurEliminatedD A D L u) := by
  simp [pressureBlockRowElimination,
    pressure_block_row_elimination A G D C L hC u p]


/-- The lower-row elimination is an equivalence of block systems, provided
    the right-hand side is transformed by the same row operation. In
    particular, this theorem does not merely identify the transformed
    operator: it proves preservation of the solution set for arbitrary
    right-hand sides f and g, under the explicit factorization C = L ∘ G. -/
theorem pressure_block_row_elimination_equiv
    (A : U →ₗ[ℝ] V) (G : P →ₗ[ℝ] V)
    (D : U →ₗ[ℝ] P) (C : P →ₗ[ℝ] P)
    (L : V →ₗ[ℝ] P)
    (hC : C = L.comp G) (u : U) (p : P)
    (f : V) (g : P) :
    (A u + G p = f ∧ D u + C p = g) ↔
      (A u + G p = f ∧
        schurEliminatedD A D L u = g - L f) := by
  constructor
  · rintro ⟨h_momentum, h_continuity⟩
    refine ⟨h_momentum, ?_⟩
    have h_elimination :=
      pressure_block_row_elimination A G D C L hC u p
    rw [h_momentum, h_continuity] at h_elimination
    exact h_elimination.symm
  · rintro ⟨h_momentum, h_transformed_continuity⟩
    refine ⟨h_momentum, ?_⟩
    have h_elimination :=
      pressure_block_row_elimination A G D C L hC u p
    rw [h_momentum] at h_elimination
    calc
      D u + C p = (D u + C p - L f) + L f := by abel
      _ = (g - L f) + L f := by
        rw [h_elimination, h_transformed_continuity]
      _ = g := by abel

end CFDX.LinearAlgebra
