#pragma once
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct LocalTimeStepControls {
    double cfl = 1.0;
    double dt_min = 1e-12;
    double dt_max = 1e12;
    double growth_limit = 2.0;
};

inline void compute_local_time_step(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& rho,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& dt,
    const LocalTimeStepControls& c = {}) {
    const std::size_t n = mesh.n_cells();
    if (mass_flux.size() != mesh.n_faces() || rho.size() != n ||
        !std::isfinite(c.cfl) || c.cfl <= 0.0 || !std::isfinite(c.dt_min) || c.dt_min <= 0.0 ||
        !std::isfinite(c.dt_max) || c.dt_max < c.dt_min) {
        throw std::invalid_argument("invalid local time-step inputs");
    }
    dt.resize(n);
    dt.fill(c.dt_max);
    const auto& cells = mesh.cells();
    const auto* faces = cells.faces_data();
    const auto* offsets = cells.offsets_data();

    for (std::size_t cell = 0; cell < n; ++cell) {
        double sum_flux = 0.0;
        const auto off = offsets[cell];
        const auto count = offsets[cell + 1] - off;
        for (std::size_t k = 0; k < count; ++k) {
            const std::size_t f = faces[off + k];
            const double flux = mass_flux(f);
            if (!std::isfinite(flux)) throw std::invalid_argument("non-finite face mass flux");
            sum_flux += std::abs(flux);
        }
        const double density = rho(cell);
        if (!std::isfinite(density) || !(density > 0.0)) {
            throw std::invalid_argument("non-positive or non-finite cell density");
        }
        if (!std::isfinite(sum_flux)) {
            throw std::invalid_argument("non-finite cell flux sum");
        }
        if (!(sum_flux > 0.0)) {
            dt(cell) = c.dt_max;
        } else {
            dt(cell) = std::clamp(c.cfl * density / sum_flux, c.dt_min, c.dt_max);
        }
    }
}

inline double global_pseudo_time_step(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& local_dt) {
    if (local_dt.empty()) throw std::invalid_argument("empty local time-step field");
    double value = local_dt(0);
    for (std::size_t i = 1; i < local_dt.size(); ++i) value = std::min(value, local_dt(i));
    return value;
}

} // namespace cfdx::physics
