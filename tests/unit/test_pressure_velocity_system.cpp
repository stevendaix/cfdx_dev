#include "cfdx/physics/pressure_velocity_system.h"

#include <cassert>
#include <string>

int main()
{
    using namespace cfdx::physics;

    PressureVelocitySystemContract contract;
    contract.n_velocity_unknowns = 12;
    contract.n_pressure_unknowns = 4;
    contract.validate();

    assert(pressure_velocity_block_name(PressureVelocityBlock::MOMENTUM) ==
           std::string("M"));
    assert(pressure_velocity_block_name(PressureVelocityBlock::PRESSURE_GRADIENT) ==
           std::string("G"));
    assert(pressure_velocity_block_name(PressureVelocityBlock::CONTINUITY) ==
           std::string("D"));
    assert(pressure_velocity_block_name(PressureVelocityBlock::PRESSURE_CONTINUITY) ==
           std::string("C"));

    assert(std::string(pressure_gauge_policy_name(
               PressureGaugePolicy::REFERENCE_CELL)) == "reference_cell");
    assert(std::string(continuity_flux_policy_name(
               ContinuityFluxPolicy::PRESSURE_CORRECTED)) == "pressure_corrected");
    assert(std::string(pressure_stabilization_policy_name(
               PressureStabilizationPolicy::RHIE_CHOW)) == "rhie_chow");

    contract.n_velocity_unknowns = 0;
    bool rejected = false;
    try {
        contract.validate();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);

    return 0;
}
