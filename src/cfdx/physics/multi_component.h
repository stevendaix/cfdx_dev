// M0.14-T03 — Multi-Component Mixture
//
// Transport d'espèces, mélange multi-composants
// M[kg/mol] = 0.02896546 pour l'air, Ru = 8.314... J/mol/K

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include <cstddef>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <stdexcept>

namespace cfdx {
namespace physics {

using cfdx::core::Location;
using cfdx::core::Field;

constexpr double R_univ = 8.314462618;

struct Species {
    std::string name;
    double M = 0.02896546;    // Masse molaire [kg/mol] (air)
    double h0 = 0.0;          // Enthalpie de formation [J/kg]
};

class Mixture {
    std::vector<Species> species_;
    std::map<std::string, std::size_t> name_to_idx_;

public:
    Mixture() = default;

    void add_species(const Species& sp) {
        name_to_idx_[sp.name] = species_.size();
        species_.push_back(sp);
    }

    std::size_t n_species() const { return species_.size(); }
    const Species& species(std::size_t i) const { return species_[i]; }
    std::size_t index(const std::string& name) const { return name_to_idx_.at(name); }

    static constexpr double air_molar_mass() { return 0.02896546; }
    static constexpr double air_gas_constant() { return R_univ / air_molar_mass(); }

    double molar_mass(const std::vector<double>& Y) const {
        double M_mix_inv = 0.0;
        for (std::size_t i = 0; i < species_.size(); ++i) {
            M_mix_inv += Y[i] / species_[i].M;
        }
        return 1.0 / M_mix_inv;
    }

    double gas_constant(const std::vector<double>& Y) const {
        return R_univ / molar_mass(Y);
    }

    double species_cp(const Species& sp, double T) const {
        (void)T;
        (void)sp;
        return 1004.5;
    }

    double mixture_cp(const std::vector<double>& Y, double T) const {
        double Cp_mix = 0.0;
        for (std::size_t i = 0; i < species_.size(); ++i) {
            Cp_mix += Y[i] * species_cp(species_[i], T);
        }
        return Cp_mix;
    }

    void validate_mass_fractions(const std::vector<double>& Y, double tol = 1e-10) const {
        double sum = 0.0;
        for (size_t i = 0; i < Y.size(); ++i) {
            if (Y[i] < -tol) {
                throw std::runtime_error("Mass fraction Y[" + std::to_string(i) + "] must be >= 0, got " + std::to_string(Y[i]));
            }
            sum += Y[i];
        }
        if (std::abs(sum - 1.0) > tol) {
            throw std::runtime_error("Sum of mass fractions must equal 1, got " + std::to_string(sum));
        }
    }
};

inline void validate_mass_fractions(const std::vector<double>& Y, double tol = 1e-10) {
    double sum = 0.0;
    for (size_t i = 0; i < Y.size(); ++i) {
        if (Y[i] < -tol) {
            throw std::runtime_error("Mass fraction Y[" + std::to_string(i) + "] must be >= 0");
        }
        sum += Y[i];
    }
    if (std::abs(sum - 1.0) > tol) {
        throw std::runtime_error("Sum of mass fractions must equal 1, got " + std::to_string(sum));
    }
}

inline void compute_mixture_properties(const Mixture& mix,
                                      const std::vector<Field<double, Location::CELL>>& Y,
                                      Field<double, Location::CELL>& M_mix,
                                      Field<double, Location::CELL>& R_mix,
                                      Field<double, Location::CELL>& Cp_mix,
                                      double T = 300.0)
{
    const std::size_t n_cells = Y[0].size();
    std::vector<double> Y_cell(mix.n_species());
    
    for (std::size_t c = 0; c < n_cells; ++c) {
        for (std::size_t i = 0; i < mix.n_species(); ++i) {
            Y_cell[i] = Y[i](c);
        }
        mix.validate_mass_fractions(Y_cell);
        M_mix(c) = mix.molar_mass(Y_cell);
        R_mix(c) = mix.gas_constant(Y_cell);
        Cp_mix(c) = mix.mixture_cp(Y_cell, T);
    }
}

}  // namespace physics
}  // namespace cfdx