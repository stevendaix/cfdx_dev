#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

/**
 * Discrete pressure-velocity block roles.
 *
 * N9 owns the pressure-velocity system contract, while N2-N4 own the
 * discretisation details and N8 owns the linear solution machinery.
 *
 * The canonical incompressible block system is
 *
 *     [ M  G ] [ U ] = [ b_U ]
 *     [ D  C ] [ p ]   [ b_p ]
 *
 * M : momentum-to-velocity operator
 * G : pressure-to-momentum operator
 * D : velocity-to-continuity operator
 * C : pressure-to-continuity contribution (including only explicitly
 *     documented pressure coupling/stabilisation terms).
 */
enum class PressureVelocityBlock {
    MOMENTUM,
    PRESSURE_GRADIENT,
    CONTINUITY,
    PRESSURE_CONTINUITY
};

enum class PressureGaugePolicy {
    NONE,
    REFERENCE_CELL,
    ZERO_MEAN
};

enum class ContinuityFluxPolicy {
    CELL_VELOCITY,
    PRESSURE_CORRECTED,
    RHIE_CHOW
};

enum class PressureStabilizationPolicy {
    NONE,
    RHIE_CHOW
};

/**
 * Immutable contract describing the discrete pressure-velocity system.
 *
 * This type deliberately contains metadata only. It does not assemble or
 * solve a matrix. That separation is essential: SIMPLE/PISO/PIMPLE/
 * fractional-step/COUPLED must operate on the same discrete contract, while
 * N8 supplies the linear algebra used to solve its blocks.
 */
struct PressureVelocitySystemContract {
    std::size_t n_velocity_unknowns = 0;
    std::size_t n_pressure_unknowns = 0;

    PressureGaugePolicy pressure_gauge = PressureGaugePolicy::REFERENCE_CELL;
    ContinuityFluxPolicy continuity_flux =
        ContinuityFluxPolicy::PRESSURE_CORRECTED;
    PressureStabilizationPolicy pressure_stabilization =
        PressureStabilizationPolicy::RHIE_CHOW;

    // Sign convention for the block equations is fixed by the matrix form
    // documented above. These flags are retained in the contract so future
    // matrix-free implementations cannot silently change conventions.
    bool momentum_block_positive_diagonal = true;
    bool pressure_gradient_is_subtracted_from_momentum = true;
    bool continuity_is_outward_flux_balance = true;

    void validate() const
    {
        if (n_velocity_unknowns == 0 || n_pressure_unknowns == 0) {
            throw std::invalid_argument(
                "pressure-velocity contract requires non-zero unknown counts");
        }

        if (pressure_gauge == PressureGaugePolicy::NONE &&
            pressure_stabilization == PressureStabilizationPolicy::NONE) {
            // A pure incompressible pressure system is expected to retain its
            // null-space. This combination is valid, but it must be explicit
            // so callers cannot accidentally infer uniqueness.
            return;
        }

        if (pressure_gauge == PressureGaugePolicy::REFERENCE_CELL &&
            n_pressure_unknowns < 1) {
            throw std::invalid_argument(
                "reference-cell pressure gauge requires pressure unknowns");
        }
    }
};

inline const char* pressure_velocity_block_name(PressureVelocityBlock block)
{
    switch (block) {
    case PressureVelocityBlock::MOMENTUM:
        return "M";
    case PressureVelocityBlock::PRESSURE_GRADIENT:
        return "G";
    case PressureVelocityBlock::CONTINUITY:
        return "D";
    case PressureVelocityBlock::PRESSURE_CONTINUITY:
        return "C";
    }
    return "unknown";
}

inline const char* pressure_gauge_policy_name(PressureGaugePolicy policy)
{
    switch (policy) {
    case PressureGaugePolicy::NONE:
        return "none";
    case PressureGaugePolicy::REFERENCE_CELL:
        return "reference_cell";
    case PressureGaugePolicy::ZERO_MEAN:
        return "zero_mean";
    }
    return "unknown";
}

inline const char* continuity_flux_policy_name(ContinuityFluxPolicy policy)
{
    switch (policy) {
    case ContinuityFluxPolicy::CELL_VELOCITY:
        return "cell_velocity";
    case ContinuityFluxPolicy::PRESSURE_CORRECTED:
        return "pressure_corrected";
    case ContinuityFluxPolicy::RHIE_CHOW:
        return "rhie_chow";
    }
    return "unknown";
}

inline const char* pressure_stabilization_policy_name(
    PressureStabilizationPolicy policy)
{
    switch (policy) {
    case PressureStabilizationPolicy::NONE:
        return "none";
    case PressureStabilizationPolicy::RHIE_CHOW:
        return "rhie_chow";
    }
    return "unknown";
}

} // namespace cfdx::physics
