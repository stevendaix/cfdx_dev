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
#include <cstddef>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>

namespace cfdx {
namespace physics {

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

// Paramètres pour gaz parfait
struct IdealGasParams {
    double R = 287.058;    // Constante des gaz [J/kg/K] (air)
    double gamma = 1.4;    // Ratio des chaleurs spécifiques Cp/Cv
    double Cp = 1004.5;    // Chaleur spécifique à pression constante [J/kg/K]
    double Cv = 717.5;     // Chaleur spécifique à volume constant [J/kg/K]
    double T_ref = 300.0;  // Température de référence [K]
    double p_ref = 101325.0; // Pression de référence [Pa]
};

// Paramètres pour fluide incompressible
struct IncompressibleParams {
    double rho = 1.0;      // Densité constante [kg/m³]
    double Cp = 4180.0;    // Chaleur spécifique [J/kg/K] (eau)
    double beta = 0.0;     // Coefficient de dilatation thermique [1/K]
    double T_ref = 300.0;  // Température de référence [K]
};

// Base class pour les EOS
class EquationOfState {
public:
    virtual ~EquationOfState() = default;
    virtual EquationOfStateType type() const = 0;

    // Propriétés thermodynamiques
    virtual double density(double p, double T) const = 0;
    virtual double enthalpy(double p, double T) const = 0;
    virtual double entropy(double p, double T) const = 0;
    virtual double speed_of_sound(double p, double T) const = 0;
    virtual double temperature_from_enthalpy(double p, double h) const = 0;
    virtual double pressure_from_density_temp(double rho, double T) const = 0;
    virtual double dp_drho_s(double p, double T) const = 0;  // Dérivée pour vitesse du son
    virtual double dp_dT_rho(double p, double T) const = 0;  // Dérivée partielle
    virtual double drho_dp_T(double p, double T) const = 0;  // Dérivée partielle
    virtual double drho_dT_p(double p, double T) const = 0;  // Dérivée partielle

    // Pour flux compressible : flux = ρ*u, energie = ρ*(e + 0.5*u²)
    virtual double internal_energy(double p, double T) const = 0;
    virtual double total_energy(double p, double T, double u_mag2) const = 0;
};

// EOS Incompressible : ρ = constant
class IncompressibleEOS : public EquationOfState {
    IncompressibleParams params_;

public:
    IncompressibleEOS(const IncompressibleParams& params = {}) : params_(params) {}

    EquationOfStateType type() const override { return EquationOfStateType::INCOMPRESSIBLE; }

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
        return 1e6;  // "Infini" pour incompressible (limite rigide)
    }

    double temperature_from_enthalpy(double /*p*/, double h) const override {
        return params_.T_ref + h / params_.Cp;
    }

    double pressure_from_density_temp(double /*rho*/, double /*T*/) const override {
        throw std::runtime_error("IncompressibleEOS: pressure not defined from rho,T");
    }

    double dp_drho_s(double /*p*/, double /*T*/) const override {
        return 1e12;  // Presque incompressible
    }

    double dp_dT_rho(double /*p*/, double /*T*/) const override {
        return 0.0;  // Incompressible: p indépendant de T à rho constant
    }

    double drho_dp_T(double /*p*/, double /*T*/) const override {
        return 0.0;  // rho constant
    }

    double drho_dT_p(double /*p*/, double /*T*/) const override {
        return -params_.rho * params_.beta;  // dilatation thermique
    }

    double internal_energy(double /*p*/, double T) const override {
        return params_.Cv * (T - params_.T_ref);
    }

    double total_energy(double /*p*/, double T, double u_mag2) const override {
        return internal_energy(0, T) + 0.5 * u_mag2;
    }

    void set_params(const IncompressibleParams& params) { params_ = params; }
    const IncompressibleParams& params() const { return params_; }
};

// EOS Gaz Parfait : p = ρ*R*T
class IdealGasEOS : public EquationOfState {
    IdealGasParams params_;

public:
    IdealGasEOS(const IdealGasParams& params = {}) : params_(params) {}

    EquationOfStateType type() const override { return EquationOfStateType::IDEAL_GAS; }

    double density(double p, double T) const override {
        return p / (params_.R * T);
    }

    double enthalpy(double /*p*/, double T) const override {
        return params_.Cp * (T - params_.T_ref);
    }

    double entropy(double p, double T) const override {
        return params_.Cp * std::log(T / params_.T_ref) - params_.R * std::log(p / params_.p_ref);
    }

    double speed_of_sound(double /*p*/, double T) const override {
        return std::sqrt(params_.gamma * params_.R * T);
    }

    double temperature_from_enthalpy(double /*p*/, double h) const override {
        return params_.T_ref + h / params_.Cp;
    }

    double pressure_from_density_temp(double rho, double T) const override {
        return rho * params_.R * T;
    }

    double dp_drho_s(double p, double T) const override {
        return params_.gamma * p / density(p, T);  // c² = γ*p/ρ
    }

    double dp_dT_rho(double /*p*/, double T) const override {
        return params_.R * T;  // p = ρ*R*T
    }

    double drho_dp_T(double p, double T) const override {
        return 1.0 / (params_.R * T);
    }

    double drho_dT_p(double p, double T) const override {
        return -p / (params_.R * T * T);
    }

    double internal_energy(double /*p*/, double T) const override {
        return params_.Cv * (T - params_.T_ref);
    }

    double total_energy(double p, double T, double u_mag2) const override {
        return internal_energy(p, T) + 0.5 * u_mag2;
    }

    void set_params(const IdealGasParams& params) { params_ = params; }
    const IdealGasParams& params() const { return params_; }
};

// Factory pour créer des EOS
inline std::unique_ptr<EquationOfState> create_eos(EquationOfStateType type) {
    switch (type) {
        case EquationOfStateType::INCOMPRESSIBLE:
            return std::make_unique<IncompressibleEOS>();
        case EquationOfStateType::IDEAL_GAS:
            return std::make_unique<IdealGasEOS>();
        default:
            throw std::runtime_error("create_eos: unsupported type");
    }
}

// Helper pour calculer les propriétés sur un champ
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