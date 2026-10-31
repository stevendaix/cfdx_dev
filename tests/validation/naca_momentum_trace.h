#pragma once

#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
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
        using cfdx::core::Location;
        const auto max_abs = [](const auto& values) {
            double result = 0.0;
            std::size_t nonfinite = 0;
            for (const double value : values) {
                if (std::isfinite(value))
                    result = std::max(result, std::abs(value));
                else
                    ++nonfinite;
            }
            return std::pair<double, std::size_t>{result, nonfinite};
        };

        const auto p_stats = max_abs(pressure.data());
        const auto mu_stats = max_abs(effective_viscosity);
        std::vector<double> grad_values;
        grad_values.reserve(3 * grad_pressure.size());
        for (std::size_t d = 0; d < 3; ++d)
            grad_values.insert(grad_values.end(),
                               grad_pressure.component_data(d).begin(),
                               grad_pressure.component_data(d).end());
        const auto grad_stats = max_abs(grad_values);

        std::vector<double> velocity_values;
        velocity_values.reserve(3 * velocity.size());
        for (std::size_t d = 0; d < 3; ++d)
            velocity_values.insert(velocity_values.end(),
                                   velocity.component_data(d).begin(),
                                   velocity.component_data(d).end());
        const auto u_stats = max_abs(velocity_values);
        const auto phi_stats = max_abs(mass_flux.data());

        std::vector<double> divergence(mesh_.n_cells(), 0.0);
        double boundary_flux = 0.0;
        double boundary_flux_abs = 0.0;
        std::size_t nonfinite_flux = 0;
        for (std::size_t face = 0; face < mesh_.n_faces(); ++face) {
            const double flux = mass_flux(face);
            if (!std::isfinite(flux)) {
                ++nonfinite_flux;
                continue;
            }
            const std::size_t owner = mesh_.ownership().owner(face);
            divergence[owner] += flux;
            const auto neighbour = mesh_.ownership().neighbour(face);
            if (neighbour >= 0) {
                divergence[static_cast<std::size_t>(neighbour)] -= flux;
            } else {
                boundary_flux += flux;
                boundary_flux_abs += std::abs(flux);
            }
        }
        const auto div_stats = max_abs(divergence);
        const std::size_t nonfinite_total =
            p_stats.second + mu_stats.second + grad_stats.second +
            u_stats.second + phi_stats.second + div_stats.second +
            nonfinite_flux;

        std::ostringstream line;
        line << std::setprecision(17)
             << "NACA_MOMENTUM_TRACE"
             << " seq=" << sequence_++
             << " iter=" << iteration
             << " max_p=" << p_stats.first
             << " max_gradp=" << grad_stats.first
             << " max_U=" << u_stats.first
             << " max_phi=" << phi_stats.first
             << " max_divphi=" << div_stats.first
             << " boundary_flux_sum=" << boundary_flux
             << " boundary_flux_abs_sum=" << boundary_flux_abs
             << " max_mu_eff=" << mu_stats.first
             << " nonfinite_pressure=" << p_stats.second
             << " nonfinite_gradp=" << grad_stats.second
             << " nonfinite_U=" << u_stats.second
             << " nonfinite_phi=" << (phi_stats.second + nonfinite_flux)
             << " nonfinite_divphi=" << div_stats.second
             << " nonfinite_mu_eff=" << mu_stats.second
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
