import Mathlib
import CFDX.FiniteVolume.Conservation

/-!
# Energy-based diffusion from a fixed linear reconstruction

This module isolates the algebraic contract for an energy-based diffusion
operator. It does not claim that a production CFDX least-squares path is
already represented by this fixed map: stencil selection, rank decisions,
boundary policies and fallback behavior must be matched separately.

For a fixed reconstruction G : Cell -> Grad and a symmetric positive
metric M on Grad, the operator is Gᵀ M G. The finite sums below make the
transpose convention explicit and avoid hiding volume/sign conventions.
-/

namespace CFDX.FiniteVolume

variable {Cell Grad : Type*} [Fintype Cell] [Fintype Grad]

/-- A fixed linear cell-to-gradient reconstruction, represented by its
matrix coefficients. The first index is a reconstructed-gradient degree
of freedom; the second is a cell degree of freedom. -/
structure FixedReconstruction (Cell Grad : Type*)
    [Fintype Cell] [Fintype Grad] where
  coeff : Grad → Cell → ℝ
  constant_annihilated : ∀ g, ∑ c : Cell, coeff g c = 0

/-- Apply the fixed reconstruction to a cell field. -/
def reconstruct
    (G : FixedReconstruction Cell Grad) (x : Cell → ℝ) : Grad → ℝ :=
  fun g => ∑ c : Cell, G.coeff g c * x c

/-- A symmetric positive-definite metric on reconstructed-gradient space.
For tensor diffusion this can represent the volume-weighted cell blocks
V_c K_c, after flattening the vector components into Grad. -/
structure SPDGradientMetric (Grad : Type*) [Fintype Grad] where
  coeff : Grad → Grad → ℝ
  symmetric : ∀ i j, coeff i j = coeff j i
  positive : ∀ y : Grad → ℝ, y ≠ 0 →
    0 < ∑ i : Grad, y i * ∑ j : Grad, coeff i j * y j

/-- Apply the gradient-space metric. -/
def metricApply
    (M : SPDGradientMetric Grad) (y : Grad → ℝ) : Grad → ℝ :=
  fun i => ∑ j : Grad, M.coeff i j * y j

/-- The energy associated with the fixed reconstruction and metric. -/
def reconstructionEnergy
    (G : FixedReconstruction Cell Grad)
    (M : SPDGradientMetric Grad) (x : Cell → ℝ) : ℝ :=
  ∑ i : Grad, reconstruct G x i * metricApply M (reconstruct G x) i

/-- The cell operator Gᵀ M G, with the transpose convention written
explicitly as a sum over reconstructed-gradient degrees of freedom. -/
def energyDiffusionOperator
    (G : FixedReconstruction Cell Grad)
    (M : SPDGradientMetric Grad) (x : Cell → ℝ) : Cell → ℝ :=
  fun c => ∑ i : Grad,
    G.coeff i c * metricApply M (reconstruct G x) i

/-- Exact energy identity for the matrix product Gᵀ M G. -/
theorem energyDiffusion_energy_identity
    (G : FixedReconstruction Cell Grad)
    (M : SPDGradientMetric Grad) (x : Cell → ℝ) :
    ∑ c : Cell, x c * energyDiffusionOperator G M x c
      = reconstructionEnergy G M x := by
  classical
  unfold energyDiffusionOperator reconstructionEnergy
    metricApply reconstruct
  simp only [Finset.mul_sum]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro i hi
  rw [← Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro c hc
  ring

/-- The energy operator is positive semidefinite. -/
theorem energyDiffusion_psd
    (G : FixedReconstruction Cell Grad)
    (M : SPDGradientMetric Grad) (x : Cell → ℝ) :
    0 ≤ ∑ c : Cell, x c * energyDiffusionOperator G M x c := by
  rw [energyDiffusion_energy_identity]
  unfold reconstructionEnergy
  by_cases hy : reconstruct G x = 0
  · simp [hy, metricApply]
  · exact le_of_lt (M.positive (reconstruct G x) hy)

/-- Constant fields are annihilated by the fixed reconstruction. -/
theorem reconstruct_constant
    (G : FixedReconstruction Cell Grad) (μ : ℝ) :
    reconstruct G (fun _ : Cell => μ) = 0 := by
  funext i
  simp [reconstruct, G.constant_annihilated i]

/-- Therefore constants are in the kernel of Gᵀ M G. -/
theorem energyDiffusion_constant
    (G : FixedReconstruction Cell Grad)
    (M : SPDGradientMetric Grad) (μ : ℝ) :
    energyDiffusionOperator G M (fun _ : Cell => μ) = 0 := by
  funext c
  simp [energyDiffusionOperator, metricApply,
    reconstruct_constant G μ]

/-- For an SPD metric, the kernel of Gᵀ M G is exactly the kernel of G.
This is an algebraic statement about the fixed map above; it does not
establish that a field-dependent production reconstruction is linear. -/
theorem energyDiffusion_kernel_iff
    (G : FixedReconstruction Cell Grad)
    (M : SPDGradientMetric Grad) (x : Cell → ℝ) :
    energyDiffusionOperator G M x = 0 ↔ reconstruct G x = 0 := by
  constructor
  · intro hA
    have henergy :
        reconstructionEnergy G M x = 0 := by
      rw [← energyDiffusion_energy_identity]
      simp [hA]
    by_contra hG
    have hpos := M.positive (reconstruct G x) hG
    unfold reconstructionEnergy metricApply at henergy
    linarith
  · intro hG
    funext c
    simp [energyDiffusionOperator, metricApply, hG]


/-- Bilinear form associated with the energy diffusion operator. -/
theorem energyDiffusion_bilinear_identity
    (G : FixedReconstruction Cell Grad)
    (M : SPDGradientMetric Grad) (x y : Cell → ℝ) :
    ∑ c : Cell, x c * energyDiffusionOperator G M y c
      = ∑ i : Grad,
          reconstruct G x i * metricApply M (reconstruct G y) i := by
  classical
  unfold energyDiffusionOperator metricApply reconstruct
  -- Expand both products of finite sums to the same four-index sum.
  simp only [Finset.mul_sum, Finset.sum_mul]
  -- The only difference is the order of the outer cell/gradient sums.
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro i hi
  apply Finset.sum_congr rfl
  intro c hc
  apply Finset.sum_congr rfl
  intro j hj
  apply Finset.sum_congr rfl
  intro d hd
  ring

/-- Symmetry of the metric pairing follows from symmetry of its coefficient
matrix. -/
theorem metric_pairing_symmetric
    (M : SPDGradientMetric Grad) (x y : Grad → ℝ) :
    ∑ i : Grad, x i * metricApply M y i
      = ∑ i : Grad, y i * metricApply M x i := by
  classical
  unfold metricApply
  calc
    (∑ i : Grad, x i * ∑ j : Grad, M.coeff i j * y j)
        = ∑ i : Grad, ∑ j : Grad, x i * (M.coeff i j * y j) := by
          simp only [Finset.mul_sum]
    _ = ∑ j : Grad, ∑ i : Grad, y j * (M.coeff j i * x i) := by
          rw [Finset.sum_comm]
          apply Finset.sum_congr rfl
          intro j hj
          apply Finset.sum_congr rfl
          intro i hi
          rw [M.symmetric i j]
          ring
    _ = ∑ i : Grad, y i * ∑ j : Grad, M.coeff i j * x j := by
          simp only [Finset.mul_sum]
          rw [Finset.sum_comm]

/-- The induced cell-space bilinear form is symmetric. This is the
operator-level symmetry contract for the fixed reconstruction model. -/
theorem energyDiffusion_bilinear_symmetric
    (G : FixedReconstruction Cell Grad)
    (M : SPDGradientMetric Grad) (x y : Cell → ℝ) :
    ∑ c : Cell, x c * energyDiffusionOperator G M y c
      = ∑ c : Cell, y c * energyDiffusionOperator G M x c := by
  rw [energyDiffusion_bilinear_identity, energyDiffusion_bilinear_identity]
  exact metric_pairing_symmetric M (reconstruct G x) (reconstruct G y)


/-! ## Exact three-cell counterexample: PSD does not imply monotonicity

The following symmetric matrix is the Gram matrix of the two-row
reconstruction
  G = [[1, 0, -1], [1/5, 1, -6/5]].
It has zero row sums and is positive semidefinite, but its (0,1)
off-diagonal coefficient is positive. Thus the energy structure alone
does not imply the non-positive off-diagonal condition of an M-matrix.
-/

/-- Coefficients of the exact two-by-three reconstruction used in the
three-cell counterexample. -/
def threeCellReconstructionCoeff (g : Fin 2) (c : Fin 3) : ℝ :=
  if g.val = 0 then
    if c.val = 0 then 1 else if c.val = 1 then 0 else -1
  else
    if c.val = 0 then 1 / 5 else if c.val = 1 then 1 else -6 / 5

/-- Explicit entries of GᵀG for the three-cell counterexample. -/
def threeCellGramEntry (i j : Fin 3) : ℝ :=
  if i.val = 0 then
    if j.val = 0 then 26 / 25 else if j.val = 1 then 1 / 5 else -31 / 25
  else if i.val = 1 then
    if j.val = 0 then 1 / 5 else if j.val = 1 then 1 else -6 / 5
  else
    if j.val = 0 then -31 / 25 else if j.val = 1 then -6 / 5 else 61 / 25

/-- The displayed matrix is exactly the Gram matrix of the stated
reconstruction, not a separately chosen example. -/
theorem threeCellGramEntry_eq_reconstruction_gram (i j : Fin 3) :
    threeCellGramEntry i j =
      ∑ g : Fin 2,
        threeCellReconstructionCoeff g i *
          threeCellReconstructionCoeff g j := by
  fin_cases i <;> fin_cases j <;>
    norm_num [threeCellGramEntry, threeCellReconstructionCoeff,
      Fin.sum_univ_succ]

/-- The counterexample matrix is symmetric. -/
theorem threeCellGramEntry_symmetric (i j : Fin 3) :
    threeCellGramEntry i j = threeCellGramEntry j i := by
  fin_cases i <;> fin_cases j <;> norm_num [threeCellGramEntry]

/-- Every row sums to zero, so constants are in the nullspace. -/
theorem threeCellGramEntry_row_sum_zero (i : Fin 3) :
    ∑ j : Fin 3, threeCellGramEntry i j = 0 := by
  fin_cases i <;>
    norm_num [threeCellGramEntry, Fin.sum_univ_succ]

/-- The quadratic form is a sum of the two reconstruction components
squared, hence nonnegative for every cell field. -/
theorem threeCellGramEntry_psd (x : Fin 3 → ℝ) :
    0 ≤ ∑ i : Fin 3,
      x i * ∑ j : Fin 3, threeCellGramEntry i j * x j := by
  simp [threeCellGramEntry, Fin.sum_univ_succ]
  nlinarith [
    sq_nonneg (x 0 - x 2),
    sq_nonneg (x 0 / 5 + x 1 - 6 * x 2 / 5)
  ]

/-- A positive off-diagonal entry violates the Z-matrix sign condition,
despite symmetry, conservation and positive semidefiniteness. -/
theorem threeCellGramEntry_positive_off_diagonal :
    0 < threeCellGramEntry (0 : Fin 3) (1 : Fin 3) := by
  norm_num [threeCellGramEntry]

end CFDX.FiniteVolume
