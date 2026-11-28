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

end CFDX.LinearAlgebra
