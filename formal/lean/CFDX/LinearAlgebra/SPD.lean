import Mathlib

/-!
# Symmetric positive-definite operator contract

This module defines the finite-dimensional mathematical contract that later
CFDX operator specifications can instantiate.

It is intentionally an abstract contract at this stage: it does not claim
that any production C++ matrix assembly has been shown to satisfy it.
-/

namespace CFDX.LinearAlgebra

/-- A finite-dimensional real matrix with symmetry and a positive
    quadratic form on every nonzero vector. -/
structure SPDMatrix (ι : Type*) [Fintype ι] where
  coeff : ι → ι → ℝ
  symmetric : ∀ i j, coeff i j = coeff j i
  positive : ∀ x : ι → ℝ, x ≠ 0 →
    0 < ∑ i, ∑ j, x i * coeff i j * x j

/-- The quadratic-form positivity guaranteed by the SPD contract. -/
theorem SPDMatrix.quadratic_form_pos
    {ι : Type*} [Fintype ι]
    (A : SPDMatrix ι) (x : ι → ℝ) (hx : x ≠ 0) :
    0 < ∑ i, ∑ j, x i * A.coeff i j * x j :=
  A.positive x hx

end CFDX.LinearAlgebra
