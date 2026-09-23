#pragma once

#include <cmath>
#include <stdexcept>

namespace cfdx::physics {

enum class ThermalConductivityModel {
    CONSTANT,
    LINEAR_TEMPERATURE,
    POWER_LAW
};

struct ThermalConductivityControls {
    ThermalConductivityModel model = ThermalConductivityModel::CONSTANT;
    double k0 = 1.0;
    double reference_temperature = 300.0;
    double slope = 0.0;
    double exponent = 0.0;
};

inline double evaluate_thermal_conductivity(
    const ThermalConductivityControls& c, double temperature)
{
    if (!std::isfinite(temperature) || !std::isfinite(c.k0) ||
        !std::isfinite(c.reference_temperature) || !std::isfinite(c.slope) ||
        !std::isfinite(c.exponent) || c.k0 < 0.0)
        throw std::invalid_argument("invalid thermal conductivity controls");

    double k = c.k0;
    switch (c.model) {
    case ThermalConductivityModel::CONSTANT:
        break;
    case ThermalConductivityModel::LINEAR_TEMPERATURE:
        k = c.k0 + c.slope * (temperature - c.reference_temperature);
        break;
    case ThermalConductivityModel::POWER_LAW:
        if (c.reference_temperature <= 0.0 || temperature <= 0.0)
            throw std::invalid_argument("power-law conductivity requires positive temperatures");
        k = c.k0 * std::pow(temperature / c.reference_temperature, c.exponent);
        break;
    }
    if (!(k >= 0.0) || !std::isfinite(k))
        throw std::invalid_argument("thermal conductivity model produced invalid conductivity");
    return k;
}

} // namespace cfdx::physics
