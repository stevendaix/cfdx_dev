import Mathlib

/-!
# Discrete energy transfer: conservation is not stability

This module isolates a minimal finite-dimensional counterexample. A two-cell
internal face contributes equal and opposite cell forces, so the global force
sum is zero. Nevertheless, its work against a chosen velocity can be positive,
negative, or zero.

The point is logical: a topological conservation theorem alone cannot imply
an energy-stability theorem. This example is not a CFDX discretization and
does not assert that physical turbulence closures must always be dissipative.
Backscatter and boundary work must be represented explicitly in a physical
energy budget.
-/

namespace CFDX.FiniteVolume

/-- Two-cell velocity used by the counterexample: (1, -1). -/
def twoCellVelocity (i : Fin 2) : ℝ :=
  if i.val = 0 then 1 else -1

/-- Equal-and-opposite internal force with positive amplitude: (1, -1). -/
def positiveWorkForce (i : Fin 2) : ℝ :=
  if i.val = 0 then 1 else -1

/-- Equal-and-opposite internal force with negative amplitude: (-1, 1). -/
def negativeWorkForce (i : Fin 2) : ℝ :=
  if i.val = 0 then -1 else 1

/-- Uniform two-cell velocity: (1, 1). -/
def uniformVelocity (_i : Fin 2) : ℝ := 1

/-- Discrete global force balance over the two cells. -/
def globalForce (force : Fin 2 → ℝ) : ℝ :=
  ∑ i : Fin 2, force i

/-- Discrete power/work in the Euclidean cell inner product. -/
def discretePower (velocity force : Fin 2 → ℝ) : ℝ :=
  ∑ i : Fin 2, velocity i * force i

theorem positive_work_force_is_conservative :
    globalForce positiveWorkForce = 0 := by
  norm_num [globalForce, positiveWorkForce]

theorem positive_work_force_injects_energy :
    discretePower twoCellVelocity positiveWorkForce = 2 := by
  norm_num [discretePower, twoCellVelocity, positiveWorkForce]

theorem negative_work_force_is_conservative :
    globalForce negativeWorkForce = 0 := by
  norm_num [globalForce, negativeWorkForce]

theorem negative_work_force_removes_energy :
    discretePower twoCellVelocity negativeWorkForce = -2 := by
  norm_num [discretePower, twoCellVelocity, negativeWorkForce]

theorem conservative_force_is_conservative_for_uniform_velocity :
    globalForce positiveWorkForce = 0 := by
  exact positive_work_force_is_conservative

theorem conservative_force_can_have_zero_power :
    discretePower uniformVelocity positiveWorkForce = 0 := by
  norm_num [discretePower, uniformVelocity, positiveWorkForce]

/-- A single counterexample disproves the implication
"zero global force balance implies non-positive discrete power". -/
theorem conservation_does_not_imply_nonpositive_power :
    ∃ force : Fin 2 → ℝ,
      globalForce force = 0 ∧ 0 < discretePower twoCellVelocity force := by
  refine ⟨positiveWorkForce, positive_work_force_is_conservative, ?_⟩
  norm_num [positive_work_force_injects_energy]

/-- Conservation is compatible with either sign of energy transfer. -/
theorem conservative_power_can_have_either_sign :
    (∃ force : Fin 2 → ℝ,
      globalForce force = 0 ∧ 0 < discretePower twoCellVelocity force) ∧
    (∃ force : Fin 2 → ℝ,
      globalForce force = 0 ∧ discretePower twoCellVelocity force < 0) := by
  constructor
  · exact conservation_does_not_imply_nonpositive_power
  · refine ⟨negativeWorkForce, negative_work_force_is_conservative, ?_⟩
    norm_num [negative_work_force_removes_energy]

/-- Finite-volume kinetic-energy work with positive cell-volume weights.
The density can be included in the weight when it is positive and cellwise
constant over the velocity degree of freedom. -/
def volumeWeightedPower
    (volume velocity force : Fin 2 → ℝ) : ℝ :=
  ∑ i : Fin 2, volume i * velocity i * force i

/-- Positive cell volumes do not rescue the implication from conservation
to non-positive power: the same conservative internal force injects energy
for every strictly positive pair of cell-volume weights. -/
theorem conservative_force_injects_volume_weighted_energy
    (volume : Fin 2 → ℝ)
    (hvolume : ∀ i, 0 < volume i) :
    globalForce positiveWorkForce = 0 ∧
      0 < volumeWeightedPower volume twoCellVelocity positiveWorkForce := by
  constructor
  · exact positive_work_force_is_conservative
  · simpa [volumeWeightedPower, twoCellVelocity, positiveWorkForce] using
      add_pos (hvolume (0 : Fin 2)) (hvolume (1 : Fin 2))

end CFDX.FiniteVolume
