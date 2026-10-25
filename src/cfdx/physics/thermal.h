#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/convection.h"
#include "cfdx/core/numerics/laplacian.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {

struct EnergyTerms {
    cfdx::core::Field<double, cfdx::core::Location::CELL> convection;
    cfdx::core::Field<double, cfdx::core::Location::CELL> conduction;
};

inline EnergyTerms evaluate_energy_terms(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Mesh& mesh, double rho, double cp, double conductivity,
    cfdx::core::InterpScheme scheme = cfdx::core::InterpScheme::UPWIND)
{
    const std::size_t n = mesh.n_cells();
    if (temperature.size() != n || temperature.dimension() != 1)
        throw std::runtime_error("evaluate_energy_terms: invalid temperature");
    if (mass_flux.size() != mesh.n_faces() || mass_flux.dimension() != 1)
        throw std::runtime_error("evaluate_energy_terms: invalid face mass flux");
    if (!(rho > 0.0) || !(cp > 0.0) || !(conductivity >= 0.0))
        throw std::invalid_argument("evaluate_energy_terms: invalid material properties");

    // F13 (audit 2026-10-24): both returned terms must share one unit
    // contract. The helper returns rate-of-change terms in K/s for
    //   dT/dt = -convection + conduction + Q/(rho*cp)
    //   convection(c) = sum_f(mdot_f * T_f) / (rho * V_c)   [K/s]
    //   conduction(c) = k/(rho*cp) * laplacian(T)            [K/s]
    // The previous contract multiplied the convective power by cp (a
    // flux-sum in W) while conduction was a specific rate in K/s, so the
    // two fields could not be added into a coherent temperature equation.
    auto convection = cfdx::core::compute_convection(temperature, mass_flux, mesh, scheme);
    const auto f13_geometry = cfdx::core::make_geometry_cache(mesh);
    for (std::size_t c = 0; c < n; ++c) {
        const double volume = f13_geometry.cell_volumes[c];
        if (!(volume > 0.0) || !std::isfinite(volume))
            throw std::runtime_error("evaluate_energy_terms: invalid cell volume");
        convection(c) /= rho * volume;
    }

    auto conduction = cfdx::core::compute_laplacian(temperature, mesh);
    for (std::size_t c = 0; c < n; ++c)
        conduction(c) *= conductivity / (rho * cp);

    return {std::move(convection), std::move(conduction)};
}

inline double heat_flux_conduction(double conductivity, double temperature_owner,
                                   double temperature_neighbour, double distance)
{
    if (conductivity < 0.0 || distance <= 0.0)
        throw std::invalid_argument("heat_flux_conduction: invalid conductivity or distance");
    return -conductivity * (temperature_neighbour - temperature_owner) / distance;
}

inline double cht_interface_conductance(double k_fluid, double k_solid,
                                        double fluid_distance, double solid_distance,
                                        double area)
{
    if (k_fluid <= 0.0 || k_solid <= 0.0 ||
        fluid_distance <= 0.0 || solid_distance <= 0.0 || area < 0.0)
        throw std::invalid_argument("cht_interface_conductance: invalid input");
    const double resistance = fluid_distance / (k_fluid * area) +
                              solid_distance / (k_solid * area);
    return 1.0 / resistance;
}

inline double cht_interface_heat_flux(double conductance,
                                      double temperature_fluid,
                                      double temperature_solid)
{
    if (conductance < 0.0)
        throw std::invalid_argument("cht_interface_heat_flux: conductance must be non-negative");
    return conductance * (temperature_fluid - temperature_solid);
}

}  // namespace cfdx::physics
