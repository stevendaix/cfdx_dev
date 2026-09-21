#pragma once
#include "cfdx/thermodynamics/thermo_state.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfdx::physics {

struct LowMachControls {
    double mach_threshold = 0.3;
    double pressure_floor = 1.0;
    bool thermal_density = true;
};

inline double low_mach_precondition_factor(double mach, const LowMachControls& c = {}) {
    if (mach < 0.0 || c.mach_threshold <= 0.0) throw std::invalid_argument("invalid low-Mach controls");
    return std::max(mach, 1e-3) / std::max(c.mach_threshold, mach);
}

inline double low_mach_density(double pressure, double temperature,
                               const cfdx::thermodynamics::IdealGasThermoModel& eos,
                               const LowMachControls& c = {}) {
    if (pressure < c.pressure_floor || temperature <= 0.0)
        throw std::invalid_argument("invalid low-Mach thermodynamic state");
    return eos.state(pressure, temperature).rho;
}

} // namespace cfdx::physics
