#pragma once

#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::validation {

// Opt-in snapshot of the exact fields consumed immediately before momentum
// assembly. This records evidence only and never modifies solver state.
class NacaMomentumTrace {
public:
    explicit NacaMomentumTrace(const cfdx::core::Mesh& mesh,
                               std::ostream& output = std::cerr)
        : mesh_(mesh), output_(output) {}

    void operator()(
        std::size_t iteration,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>& grad_pressure,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
        const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
        const std::vector<double>& effective_viscosity)
    {
        struct Extremum {
            double value = 0.0;
            std::size_t index = 0;
            std::size_t nonfinite = 0;
            std::size_t first_nonfinite = static_cast<std::size_t>(-1);
        };
        const auto max_abs = [](const auto& values) {
            Extremum result;
            for (std::size_t i = 0; i < values.size(); ++i) {
                const double value = values[i];
                if (std::isfinite(value)) {
                    if (std::abs(value) > result.value) {
                        result.value = std::abs(value);
                        result.index = i;
                    }
                } else {
                    ++result.nonfinite;
                    if (result.first_nonfinite == static_cast<std::size_t>(-1))
                        result.first_nonfinite = i;
                }
            }
            return result;
        };

        std::vector<double> pressure_values(pressure.size());
        for (std::size_t i = 0; i < pressure.size(); ++i)
            pressure_values[i] = pressure(i);
        const auto p_stats = max_abs(pressure_values);
        const auto mu_stats = max_abs(effective_viscosity);

        std::vector<double> grad_values;
        grad_values.reserve(3 * grad_pressure.size());
        std::vector<double> velocity_values;
        velocity_values.reserve(3 * velocity.size());
        for (std::size_t d = 0; d < 3; ++d) {
            for (std::size_t i = 0; i < grad_pressure.size(); ++i)
                grad_values.push_back(grad_pressure.component_data(d)[i]);
            for (std::size_t i = 0; i < velocity.size(); ++i)
                velocity_values.push_back(velocity.component_data(d)[i]);
        }
        const auto grad_stats = max_abs(grad_values);
        const auto u_stats = max_abs(velocity_values);

        std::vector<double> phi_values(mass_flux.size());
        for (std::size_t i = 0; i < mass_flux.size(); ++i)
            phi_values[i] = mass_flux(i);
        const auto phi_stats = max_abs(phi_values);

        std::vector<double> divergence(mesh_.n_cells(), 0.0);
        double boundary_flux = 0.0;
        double boundary_flux_abs = 0.0;
        std::size_t nonfinite_flux = 0;
        std::size_t first_nonfinite_flux = static_cast<std::size_t>(-1);
        for (std::size_t face = 0; face < mesh_.n_faces(); ++face) {
            const double flux = mass_flux(face);
            if (!std::isfinite(flux)) {
                ++nonfinite_flux;
                if (first_nonfinite_flux == static_cast<std::size_t>(-1))
                    first_nonfinite_flux = face;
                continue;
            }
            const std::size_t owner = mesh_.ownership().owner(face);
            divergence[owner] += flux;
            const auto neighbour = mesh_.ownership().neighbour(face);
            if (neighbour >= 0)
                divergence[static_cast<std::size_t>(neighbour)] -= flux;
            else {
                boundary_flux += flux;
                boundary_flux_abs += std::abs(flux);
            }
        }
        const auto div_stats = max_abs(divergence);
        const std::size_t nonfinite_total =
            p_stats.nonfinite + mu_stats.nonfinite + grad_stats.nonfinite +
            u_stats.nonfinite + phi_stats.nonfinite + div_stats.nonfinite +
            nonfinite_flux;

        const auto face_description = [&](std::size_t face) {
            if (face >= mesh_.n_faces())
                return std::string("none");
            std::ostringstream value;
            const auto owner = mesh_.ownership().owner(face);
            const auto neighbour = mesh_.ownership().neighbour(face);
            value << face << ":" << owner << ":" << neighbour;
            return value.str();
        };

        std::ostringstream line;
        line << std::setprecision(17)
             << "NACA_MOMENTUM_TRACE"
             << " seq=" << sequence_++
             << " iter=" << iteration
             << " max_p=" << p_stats.value
             << " max_p_cell=" << p_stats.index
             << " max_gradp=" << grad_stats.value
             << " max_gradp_component=" << grad_stats.index / std::max<std::size_t>(grad_pressure.size(), 1)
             << " max_gradp_cell=" << grad_stats.index % std::max<std::size_t>(grad_pressure.size(), 1)
             << " max_U=" << u_stats.value
             << " max_U_component=" << u_stats.index / std::max<std::size_t>(velocity.size(), 1)
             << " max_U_cell=" << u_stats.index % std::max<std::size_t>(velocity.size(), 1)
             << " max_phi=" << phi_stats.value
             << " max_phi_face_owner_neighbour=" << face_description(phi_stats.index)
             << " max_divphi=" << div_stats.value
             << " max_divphi_cell=" << div_stats.index
             << " boundary_flux_sum=" << boundary_flux
             << " boundary_flux_abs_sum=" << boundary_flux_abs
             << " max_mu_eff=" << mu_stats.value
             << " nonfinite_pressure=" << p_stats.nonfinite
             << " first_nonfinite_pressure_cell=" << p_stats.first_nonfinite
             << " nonfinite_gradp=" << grad_stats.nonfinite
             << " first_nonfinite_gradp_component_cell="
             << grad_stats.first_nonfinite
             << " nonfinite_U=" << u_stats.nonfinite
             << " first_nonfinite_U_component_cell=" << u_stats.first_nonfinite
             << " nonfinite_phi=" << (phi_stats.nonfinite + nonfinite_flux)
             << " first_nonfinite_phi_face="
             << (phi_stats.first_nonfinite != static_cast<std::size_t>(-1)
                     ? phi_stats.first_nonfinite : first_nonfinite_flux)
             << " nonfinite_divphi=" << div_stats.nonfinite
             << " first_nonfinite_divphi_cell=" << div_stats.first_nonfinite
             << " nonfinite_mu_eff=" << mu_stats.nonfinite
             << " nonfinite_total=" << nonfinite_total
             << '\n';
        output_ << line.str();
        output_.flush();
    }

private:
    const cfdx::core::Mesh& mesh_;
    std::ostream& output_;
    std::size_t sequence_ = 0;
};

} // namespace cfdx::validation
