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
