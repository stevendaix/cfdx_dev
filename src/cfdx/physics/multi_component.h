// M0.14-T03 — Multi-Component Mixture
//
// Transport d'espèces, mélange multi-composants

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include <cstddef>
#include <vector>
#include <string>
#include <map>

namespace cfdx {
namespace physics {

struct Species {
    std::string name;
    double M = 28.97;      // Masse molaire [g/mol]
    double h0 = 0.0;       // Enthalpie de formation [J/kg]
};

// Mélange de gaz parfaits
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

    // Masse molaire du mélange [g/mol]
    double molar_mass(const std::vector<double>& Y) const {
        double M_mix = 0.0;
        for (std::size_t i = 0; i < species_.size(); ++i) {
            M_mix += Y[i] / species_[i].M;
        }
        return 1.0 / M_mix * 1000.0;  // [kg/mol]
    }

    // Constante des gaz du mélange [J/kg/K]
    double gas_constant(const std::vector<double>& Y) const {
        const double R_univ = 8.314462618;
        return R_univ / molar_mass(Y);
    }

    // Cp du mélange [J/kg/K]
    double mixture_cp(const std::vector<double>& Y, double T) const {
        double Cp_mix = 0.0;
        for (std::size_t i = 0; i < species_.size(); ++i) {
            Cp_mix += Y[i] * species_cp(species_[i], T);
        }
        return Cp_mix;
    }

    // Cp d'une espèce (approximation polynomiale NASA)
    double species_cp(const Species& sp, double T) const {
        // Simplifié : Cp constant par espèce
        (void)sp; (void)T;
        return 1004.5;
    }
};

// Compute champs de fractions massiques et propriétés
inline void compute_mixture_properties(const Mixture& mix,
                                        const std::vector<Field<double, Location::CELL>>& Y,
                                        Field<double, Location::CELL>& M_mix,
                                        Field<double, Location::CELL>& R_mix,
                                        Field<double, Location::CELL>& Cp_mix,
                                        double T = 300.0)
{
    const std::size_t n_cells = Y[0].size();
    for (std::size_t c = 0; c < n_cells; ++c) {
        std::vector<double> Y_cell(mix.n_species());
        for (std::size_t i = 0; i < mix.n_species(); ++i) {
            Y_cell[i] = Y[i](c);
        }
        M_mix(c) = mix.molar_mass(Y_cell);
        R_mix(c) = mix.gas_constant(Y_cell);
        Cp_mix(c) = mix.mixture_cp(Y_cell, T);
    }
}

}  // namespace physics
}  // namespace cfdx