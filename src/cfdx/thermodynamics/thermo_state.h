#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfdx::thermodynamics {

struct ThermoState {
    double rho = 1.0;
    double mu = 1e-5;
    double cp = 1000.0;
    double conductivity = 0.02;
    double enthalpy = 0.0;
    double speed_of_sound = 0.0;

    bool is_valid() const
    {
        return std::isfinite(rho) && std::isfinite(mu) &&
               std::isfinite(cp) && std::isfinite(conductivity) &&
               std::isfinite(enthalpy) && std::isfinite(speed_of_sound) &&
               rho > 0.0 && mu >= 0.0 && cp > 0.0 &&
               conductivity >= 0.0 && speed_of_sound >= 0.0;
    }
};

struct IdealGasThermoModel {
    double R = 287.05;
    double gamma = 1.4;
    double cp = 1.4 * 287.05 / 0.4;
    double mu_ref = 1.716e-5;
    double T_ref = 273.15;
    double S = 110.4;
    double k_ref = 0.0241;

    void validate_parameters() const
    {
        if (!std::isfinite(R) || R <= 0.0 ||
            !std::isfinite(gamma) || gamma <= 1.0 ||
            !std::isfinite(cp) || cp <= 0.0 ||
            !std::isfinite(mu_ref) || mu_ref <= 0.0 ||
            !std::isfinite(T_ref) || T_ref <= 0.0 ||
            !std::isfinite(S) || S < 0.0 ||
            !std::isfinite(k_ref) || k_ref <= 0.0) {
            throw std::invalid_argument("IdealGasThermoModel: invalid parameters");
        }

        const double cp_expected = gamma * R / (gamma - 1.0);
        const double scale = std::max({1.0, std::abs(cp), std::abs(cp_expected)});
        if (std::abs(cp - cp_expected) > 1e-12 * scale)
            throw std::invalid_argument(
                "IdealGasThermoModel: cp must equal gamma*R/(gamma-1)");
    }

    ThermoState state(double p, double T) const
    {
        validate_parameters();
        if (!std::isfinite(p) || !std::isfinite(T) || p <= 0.0 || T <= 0.0)
            throw std::invalid_argument("IdealGasThermoModel: invalid thermodynamic state");

        ThermoState s;
        s.rho = p / (R * T);
        s.cp = cp;
        s.mu = mu_ref * std::pow(T / T_ref, 1.5) * (T_ref + S) / (T + S);
        s.conductivity = k_ref * std::pow(T / T_ref, 0.76);
        s.enthalpy = cp * T;
        s.speed_of_sound = std::sqrt(gamma * R * T);
        if (!s.is_valid())
            throw std::runtime_error("IdealGasThermoModel: non-finite state result");
        return s;
    }
};

} // namespace cfdx::thermodynamics
