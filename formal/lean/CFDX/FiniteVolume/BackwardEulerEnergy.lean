import Mathlib
import CFDX.FiniteVolume.VectorFaceWork

/-!
# Backward-Euler kinetic-energy work identity

For a positive density and time step, the backward-Euler momentum increment
paired with the new velocity splits exactly into the kinetic-energy increment
plus a non-negative quadratic time-discretization defect.

This is a time-term contract only. It does not establish stability of the
coupled CFDX pressure-velocity algorithm, nor does it identify pseudo-time
or under-relaxation terms with physical time derivatives.
-/

namespace CFDX.FiniteVolume

/-- Kinetic energy per unit mass for a three-component velocity. -/
def kineticEnergy3 (velocity : Fin 3 → ℝ) : ℝ :=
  (1 / 2 : ℝ) * dot3 velocity velocity

/-- Backward-Euler momentum increment paired with the new velocity,
per unit volume, for density rho and time step dt. -/
def backwardEulerKineticPower
    (rho dt : ℝ) (newVelocity oldVelocity : Fin 3 → ℝ) : ℝ :=
  rho / dt * dot3 newVelocity (fun i => newVelocity i - oldVelocity i)

/-- The backward-Euler kinetic-energy identity. The final term is the
quadratic numerical dissipation associated with the implicit time step. -/
theorem backward_euler_kinetic_energy_identity
    (rho dt : ℝ) (newVelocity oldVelocity : Fin 3 → ℝ) :
    backwardEulerKineticPower rho dt newVelocity oldVelocity =
      rho / dt *
        (kineticEnergy3 newVelocity - kineticEnergy3 oldVelocity +
          (1 / 2 : ℝ) *
            dot3 (fun i => newVelocity i - oldVelocity i)
              (fun i => newVelocity i - oldVelocity i)) := by
  unfold backwardEulerKineticPower kineticEnergy3 dot3
  ring

/-- Squared Euclidean norm of a three-component real vector is non-negative. -/
theorem dot3_self_nonnegative (velocity : Fin 3 → ℝ) :
    0 ≤ dot3 velocity velocity := by
  unfold dot3
  positivity

/-- With positive density and time step, the backward-Euler defect is
non-negative. This is a local algebraic property, not a full solver-stability
theorem. -/
theorem backward_euler_defect_nonnegative
    (rho dt : ℝ) (newVelocity oldVelocity : Fin 3 → ℝ)
    (hrho : 0 < rho) (hdt : 0 < dt) :
    0 ≤ rho / dt * ((1 / 2 : ℝ) *
      dot3 (fun i => newVelocity i - oldVelocity i)
        (fun i => newVelocity i - oldVelocity i)) := by
  have hnorm := dot3_self_nonnegative
    (fun i => newVelocity i - oldVelocity i)
  have hscale : 0 ≤ rho / dt := div_nonneg hrho.le hdt.le
  positivity

end CFDX.FiniteVolume


/-- The BDF2 G-energy for two consecutive three-component velocity states.
It is a sum of squared Euclidean norms and is therefore non-negative. -/
def bdf2GEnergy3
    (newVelocity currentVelocity : Fin 3 → ℝ) : ℝ :=
  dot3 newVelocity newVelocity +
    dot3 (fun i => 2 * newVelocity i - currentVelocity i)
      (fun i => 2 * newVelocity i - currentVelocity i)

/-- The second temporal difference in the constant-step BDF2 formula. -/
def bdf2SecondDifference3
    (newVelocity currentVelocity previousVelocity : Fin 3 → ℝ) :
    Fin 3 → ℝ :=
  fun i => newVelocity i - 2 * currentVelocity i + previousVelocity i

/-- Constant-step BDF2 momentum increment paired with the new velocity,
per unit volume. -/
def constantStepBdf2KineticPower
    (rho dt : ℝ)
    (newVelocity currentVelocity previousVelocity : Fin 3 → ℝ) : ℝ :=
  rho / (2 * dt) *
    dot3 newVelocity
      (fun i => 3 * newVelocity i - 4 * currentVelocity i + previousVelocity i)

/-- G-stability identity for constant-step BDF2. The time derivative paired
with the new velocity is the increment of a two-level G-energy plus a
non-negative second-difference defect. This theorem is deliberately limited
to constant time step; variable-step BDF2 requires a step-ratio-dependent
energy estimate and is not covered here. -/
theorem constant_step_bdf2_g_stability_identity
    (rho dt : ℝ)
    (newVelocity currentVelocity previousVelocity : Fin 3 → ℝ) :
    constantStepBdf2KineticPower rho dt
      newVelocity currentVelocity previousVelocity =
      rho / (4 * dt) *
        (bdf2GEnergy3 newVelocity currentVelocity -
          bdf2GEnergy3 currentVelocity previousVelocity +
          dot3 (bdf2SecondDifference3
            newVelocity currentVelocity previousVelocity)
            (bdf2SecondDifference3
              newVelocity currentVelocity previousVelocity)) := by
  unfold constantStepBdf2KineticPower bdf2GEnergy3
    bdf2SecondDifference3 dot3
  ring

/-- The constant-step BDF2 second-difference contribution is non-negative
for positive density and time step. -/
theorem constant_step_bdf2_defect_nonnegative
    (rho dt : ℝ)
    (newVelocity currentVelocity previousVelocity : Fin 3 → ℝ)
    (hrho : 0 < rho) (hdt : 0 < dt) :
    0 ≤ rho / (4 * dt) *
      dot3 (bdf2SecondDifference3
        newVelocity currentVelocity previousVelocity)
        (bdf2SecondDifference3
          newVelocity currentVelocity previousVelocity) := by
  have hnorm := dot3_self_nonnegative
    (bdf2SecondDifference3 newVelocity currentVelocity previousVelocity)
  have hscale : 0 ≤ rho / (4 * dt) := by positivity
  positivity
