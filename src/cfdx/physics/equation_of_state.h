// M0.14-T01 — Equations of State
//
// Spécification CFDX v0.7 §36 :
//   Modèles thermodynamiques pour ρ, p, T, h, c
//
//   Modèles supportés :
//     - Incompressible (densité constante)
//     - Ideal Gas (loi des gaz parfaits)
//     - Peng-Robinson (pour haute pression) - futur
//     - IAPWS-IF97 (vapeur d'eau) - futur
//
//   Interface unifiée :
//     ρ = f(p, T, composition)
//     h = f(p, T, composition)
//     c = sqrt(∂p/∂ρ)_s

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/index_types.h"
#include <algorithm>
#include <cstddef>
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <cmath>
#include <limits>

namespace cfdx {
namespace physics {

using cfdx::core::Location;
using cfdx::core::Field;

enum class EquationOfStateType : std::uint8_t {
    INCOMPRESSIBLE = 0,
    IDEAL_GAS,
    PENG_ROBINSON,
    IAPWS_IF97
};

inline const char* to_string(EquationOfStateType t) {
    switch (t) {
        case EquationOfStateType::INCOMPRESSIBLE: return "incompressible";
        case EquationOfStateType::IDEAL_GAS:      return "ideal_gas";
        case EquationOfStateType::PENG_ROBINSON:  return "peng_robinson";
        case EquationOfStateType::IAPWS_IF97:     return "iapws_if97";
        default:                                  return "unknown";
    }
}

struct IdealGasParams {
    double M = 0.02896546;     // Mass molaire [kg/mol] (air)
    double gamma = 1.4;        // Ratio des chaleurs spécifiques Cp/Cv
    double R_univ = 8.314462618;  // Constante des gaz universelle [J/mol/K]
    double Cp = 1004.5;        // Chaleur spécifique à pression constante [J/kg/K]
    double T_ref = 300.0;      // Température de référence [K]
    double p_ref = 101325.0;   // Pression de référence [Pa]
};

struct IncompressibleParams {
    double rho = 1.0;          // Densité constante [kg/m³]
    double Cp = 4180.0;        // Chaleur spécifique [J/kg/K] (eau)
    double beta = 0.0;         // Coefficient de dilatation thermique [1/K]
    double T_ref = 300.0;      // Température de référence [K]
};

class EquationOfState {
public:
    virtual ~EquationOfState() = default;
    virtual EquationOfStateType type() const = 0;

    virtual double density(double p, double T) const = 0;
    virtual double enthalpy(double p, double T) const = 0;
    virtual double entropy(double p, double T) const = 0;
    virtual double speed_of_sound(double p, double T) const = 0;
    virtual double temperature_from_enthalpy(double p, double h) const = 0;
    virtual double pressure_from_density_temp(double rho, double T) const = 0;
    virtual double dp_drho_s(double p, double T) const = 0;
    virtual double dp_dT_rho(double p, double T) const = 0;
    virtual double drho_dp_T(double p, double T) const = 0;
    virtual double drho_dT_p(double p, double T) const = 0;

    virtual double cp(double p, double T) const = 0;
    virtual double cv(double p, double T) const = 0;

    virtual double internal_energy(double p, double T) const = 0;
    virtual double total_energy(double p, double T, double u_mag2) const = 0;
};

class IncompressibleEOS : public EquationOfState {
    IncompressibleParams params_;

public:
    IncompressibleEOS(const IncompressibleParams& params = {}) : params_(params) { validate_params(); }

    EquationOfStateType type() const override { return EquationOfStateType::INCOMPRESSIBLE; }

    void validate_params() const {
        if (!std::isfinite(params_.rho) || !(params_.rho > 0.0) ||
            !std::isfinite(params_.Cp) || !(params_.Cp > 0.0) ||
            !std::isfinite(params_.T_ref) || !(params_.T_ref > 0.0) ||
            !std::isfinite(params_.beta) || params_.beta < 0.0) {
            throw std::invalid_argument("invalid incompressible EOS parameters");
        }
    }

    double density(double /*p*/, double /*T*/) const override {
        return params_.rho;
    }

    double enthalpy(double /*p*/, double T) const override {
        return params_.Cp * (T - params_.T_ref);
    }

    double entropy(double /*p*/, double T) const override {
        return params_.Cp * std::log(T / params_.T_ref);
    }

    double speed_of_sound(double /*p*/, double /*T*/) const override {
        if (params_.beta <= 0) {
            return std::numeric_limits<double>::infinity();
        }
        return std::sqrt(params_.rho / params_.beta);
    }

    double temperature_from_enthalpy(double /*p*/, double h) const override {
        return params_.T_ref + h / params_.Cp;
    }

    double pressure_from_density_temp(double /*rho*/, double /*T*/) const override {
        throw std::runtime_error("IncompressibleEOS: pressure not defined from rho,T");
    }

    double dp_drho_s(double /*p*/, double /*T*/) const override {
        return 1.0 / (params_.beta * params_.rho);
    }

    double dp_dT_rho(double /*p*/, double /*T*/) const override {
        return params_.rho * params_.beta;
    }

    double drho_dp_T(double /*p*/, double /*T*/) const override {
        return 0.0;
    }

    double drho_dT_p(double /*p*/, double /*T*/) const override {
        return -params_.rho * params_.beta;
    }

    double internal_energy(double /*p*/, double T) const override {
        return params_.Cp * (T - params_.T_ref);
    }

    double cp(double /*p*/, double /*T*/) const override {
        return params_.Cp;
    }

    double cv(double /*p*/, double /*T*/) const override {
        return params_.Cp;
    }

    double total_energy(double /*p*/, double T, double u_mag2) const override {
        return internal_energy(0, T) + 0.5 * u_mag2;
    }

    void set_params(const IncompressibleParams& params) { params_ = params; validate_params(); }
    const IncompressibleParams& params() const { return params_; }
};

class IdealGasEOS : public EquationOfState {
    IdealGasParams params_;
    double R = 0.0;

    void update_R() {
        R = params_.R_univ / params_.M;
    }

    void validate_thermodynamic_consistency() const {
        const double cp_expected = params_.gamma * R / (params_.gamma - 1.0);
        const double scale = std::max({1.0, std::abs(params_.Cp), std::abs(cp_expected)});
        if (std::abs(params_.Cp - cp_expected) > 1.0e-10 * scale)
            throw std::invalid_argument("IdealGasEOS: Cp is inconsistent with gamma and R");
    }

public:
    IdealGasEOS(const IdealGasParams& params = {}) : params_(params) {
        validate_params();
        update_R();
        validate_thermodynamic_consistency();
    }

    EquationOfStateType type() const override { return EquationOfStateType::IDEAL_GAS; }

    void validate_params() const {
        if (!std::isfinite(params_.M) || !(params_.M > 0.0) ||
            !std::isfinite(params_.gamma) || !(params_.gamma > 1.0) ||
            !std::isfinite(params_.R_univ) || !(params_.R_univ > 0.0) ||
            !std::isfinite(params_.Cp) || !(params_.Cp > 0.0) ||
            !std::isfinite(params_.T_ref) || !(params_.T_ref > 0.0) ||
            !std::isfinite(params_.p_ref) || !(params_.p_ref > 0.0)) {
            throw std::invalid_argument("invalid ideal-gas EOS parameters");
        }
    }

    void set_params(const IdealGasParams& params) {
        params_ = params;
        validate_params();
        update_R();
        validate_thermodynamic_consistency();
    }

    const IdealGasParams& params() const { return params_; }

    double density(double p, double T) const override {
        return p / (R * T);
    }

    double enthalpy(double /*p*/, double T) const override {
        double Cp_derived = params_.gamma * R / (params_.gamma - 1.0);
        return Cp_derived * (T - params_.T_ref);
    }

    double entropy(double p, double T) const override {
        double Cp_derived = params_.gamma * R / (params_.gamma - 1.0);
        return Cp_derived * std::log(T / params_.T_ref) - R * std::log(p / params_.p_ref);
    }

    double speed_of_sound(double /*p*/, double T) const override {
        double gamma = params_.gamma;
        return std::sqrt(gamma * R * T);
    }

    double temperature_from_enthalpy(double /*p*/, double h) const override {
        double Cp_derived = params_.gamma * R / (params_.gamma - 1.0);
        return params_.T_ref + h / Cp_derived;
    }

    double pressure_from_density_temp(double rho, double T) const override {
        return rho * R * T;
    }

    double dp_drho_s(double p, double T) const override {
        return params_.gamma * p / density(p, T);
    }

    double dp_dT_rho(double p, double T) const override {
        return density(p, T) * R;
    }

    double drho_dp_T(double p, double T) const override {
        return 1.0 / (R * T);
    }

    double drho_dT_p(double p, double T) const override {
        return -p / (R * T * T);
    }

    double internal_energy(double /*p*/, double T) const override {
        double Cv_derived = R / (params_.gamma - 1.0);
        return Cv_derived * (T - params_.T_ref);
    }

    double cp(double /*p*/, double /*T*/) const override {
        return params_.gamma * R / (params_.gamma - 1.0);
    }

    double cv(double /*p*/, double /*T*/) const override {
        return R / (params_.gamma - 1.0);
    }

    double total_energy(double p, double T, double u_mag2) const override {
        return internal_energy(p, T) + 0.5 * u_mag2;
    }
};

inline std::unique_ptr<EquationOfState> create_eos(EquationOfStateType type,
                                                       const IdealGasParams& ideal_params = {}) {
    switch (type) {
        case EquationOfStateType::INCOMPRESSIBLE:
            return std::make_unique<IncompressibleEOS>();
        case EquationOfStateType::IDEAL_GAS:
            return std::make_unique<IdealGasEOS>(ideal_params);
        default:
            throw std::runtime_error("create_eos: unsupported type");
    }
}

inline std::unique_ptr<EquationOfState> create_eos(EquationOfStateType type,
                                                       const IncompressibleParams& incompressible_params) {
    switch (type) {
        case EquationOfStateType::INCOMPRESSIBLE:
            return std::make_unique<IncompressibleEOS>(incompressible_params);
        case EquationOfStateType::IDEAL_GAS:
            return std::make_unique<IdealGasEOS>();
        default:
            throw std::runtime_error("create_eos: unsupported type for parameters");
    }
}

template <typename EOS>
void compute_thermo_fields(const EOS& eos,
                           const Field<double, Location::CELL>& pressure,
                           const Field<double, Location::CELL>& temperature,
                           Field<double, Location::CELL>& density,
                           Field<double, Location::CELL>& enthalpy,
                           Field<double, Location::CELL>& entropy,
                           Field<double, Location::CELL>& speed_of_sound)
{
    const std::size_t n = pressure.size();
    for (std::size_t c = 0; c < n; ++c) {
        const double p = pressure(c);
        const double T = temperature(c);
        density(c) = eos.density(p, T);
        enthalpy(c) = eos.enthalpy(p, T);
        entropy(c) = eos.entropy(p, T);
        speed_of_sound(c) = eos.speed_of_sound(p, T);
    }
}

}  // namespace physics
}  // namespace cfdx