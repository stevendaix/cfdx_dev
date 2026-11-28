import Mathlib
import CFDX.FiniteVolume.MeshVectorWork

/-!
# Vector momentum residual power

This module lifts the owner/neighbour face-work identity to an integrated
vector momentum residual with arbitrary cell source densities. A shared
oriented face flux contributes with owner-positive/neighbour-negative
incidence; cell sources (for example pressure-gradient, body-force and
temporal terms) are multiplied by cell volume.

The theorem is an accounting identity over a finite abstract topology. It
does not prove a sign for pressure work, viscous work, closure work, or the
complete kinetic-energy balance. It also does not establish that the C++
momentum assembly or its runtime mesh conforms to this abstraction.
-/

namespace CFDX.FiniteVolume

/-- Integrated vector momentum residual in one cell. faceFlux is the
shared owner-oriented integrated vector flux; sourceDensity contains
per-volume vector sources with their physical signs already applied. -/
def cellVectorIntegratedMomentumResidual
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (faceFlux : Face → Fin 3 → ℝ)
    (sourceDensity : Cell → Fin 3 → ℝ)
    (volume : Cell → ℝ) (cell : Cell) : Fin 3 → ℝ :=
  fun i =>
    cellVectorFaceForce topology faceFlux cell i +
      volume cell * sourceDensity cell i

/-- Global velocity pairing with the integrated cell momentum residual.
Its physical dimension is power when the inputs carry consistent SI units. -/
def meshIntegratedMomentumResidualPower
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (faceFlux : Face → Fin 3 → ℝ)
    (sourceDensity : Cell → Fin 3 → ℝ)
    (volume : Cell → ℝ) : ℝ :=
  ∑ cell : Cell,
    dot3 (velocity cell)
      (cellVectorIntegratedMomentumResidual
        topology faceFlux sourceDensity volume cell)

/-- Dot product distributes over vector addition. -/
theorem dot3_add_right
    (a b c : Fin 3 → ℝ) :
    dot3 a (fun i => b i + c i) = dot3 a b + dot3 a c := by
  simp [dot3]
  ring

/-- A scalar cell-volume factor can be taken outside the vector dot product. -/
theorem dot3_smul_right
    (a b : Fin 3 → ℝ) (scale : ℝ) :
    dot3 a (fun i => scale * b i) = scale * dot3 a b := by
  simp [dot3]
  ring

/-- The integrated momentum residual power is exactly the sum of
owner/neighbour face work and cell-source work.

Internal faces reduce to the owner-minus-neighbour velocity jump dotted with
faceFlux; boundary faces reduce to owner velocity dotted with faceFlux.
Pressure, body-force, temporal and other cell sources remain explicit in the
second sum. No sign restriction is imposed on any contribution. -/
theorem mesh_integrated_momentum_power_eq_face_and_source_work
    {Cell Face : Type*} [Fintype Cell] [Fintype Face]
    (topology : OwnerNeighbourTopology Cell Face)
    (velocity : Cell → Fin 3 → ℝ)
    (faceFlux : Face → Fin 3 → ℝ)
    (sourceDensity : Cell → Fin 3 → ℝ)
    (volume : Cell → ℝ) :
    meshIntegratedMomentumResidualPower
        topology velocity faceFlux sourceDensity volume =
      (∑ face : Face,
        ownerNeighbourFaceWork topology velocity faceFlux face) +
      (∑ cell : Cell,
        volume cell * dot3 (velocity cell) (sourceDensity cell)) := by
  classical
  calc
    meshIntegratedMomentumResidualPower
        topology velocity faceFlux sourceDensity volume =
      ∑ cell : Cell,
        (dot3 (velocity cell)
            (cellVectorFaceForce topology faceFlux cell) +
          volume cell * dot3 (velocity cell) (sourceDensity cell)) := by
        unfold meshIntegratedMomentumResidualPower
          cellVectorIntegratedMomentumResidual
        apply Finset.sum_congr rfl
        intro cell hcell
        rw [dot3_add_right, dot3_smul_right]
    _ =
      meshCellVectorWork topology velocity faceFlux +
        (∑ cell : Cell,
          volume cell * dot3 (velocity cell) (sourceDensity cell)) := by
        rw [Finset.sum_add_distrib]
        rfl
    _ =
      (∑ face : Face,
        ownerNeighbourFaceWork topology velocity faceFlux face) +
      (∑ cell : Cell,
        volume cell * dot3 (velocity cell) (sourceDensity cell)) := by
        rw [mesh_vector_work_eq_owner_neighbour_face_work]

end CFDX.FiniteVolume
