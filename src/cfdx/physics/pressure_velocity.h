#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/mesh/mesh.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct PressureCorrectionSystem {
    cfdx::core::SparseMatrix matrix;
    cfdx::core::Vector rhs;
};

inline PressureCorrectionSystem assemble_pressure_correction(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& momentum_diagonal,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& continuity)
{
    using namespace cfdx::core;
    const std::size_t n = mesh.n_cells();
    if (momentum_diagonal.dimension() != 1 || momentum_diagonal.size() != n ||
        continuity.dimension() != 1 || continuity.size() != n) {
        throw std::runtime_error("assemble_pressure_correction: invalid fields");
    }

    std::vector<Vec3> centres(mesh.n_faces());
    std::vector<Vec3> sf(mesh.n_faces());
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const auto off = mesh.faces().offsets_data()[f];
        const auto count = mesh.faces().offsets_data()[f + 1] - off;
        const auto fg = compute_face_geometry(
            mesh.points().x_data(), mesh.points().y_data(), mesh.points().z_data(),
            mesh.faces().vertices_data(), off, count);
        centres[f] = fg.centre;
        sf[f] = fg.Sf;
    }

    std::vector<double> volume(n, 0.0);
    std::vector<Vec3> cell_centre(n);
    for (std::size_t c = 0; c < n; ++c) {
        const auto off = mesh.cells().offsets_data()[c];
        const auto count = mesh.cells().offsets_data()[c + 1] - off;
        const auto cg = compute_cell_geometry(
            centres.data(), sf.data(), mesh.cells().faces_data() + off, count);
        volume[c] = cg.volume;
        cell_centre[c] = cg.centre;
    }

    std::vector<std::vector<std::pair<std::size_t, double>>> rows(n);
    std::vector<double> diagonal(n, 0.0);
    const auto& own = mesh.ownership();
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const auto nb_raw = own.neighbour(f);
        if (nb_raw < 0) continue;
        const std::size_t o = own.owner(f);
        const std::size_t nb = static_cast<std::size_t>(nb_raw);
        if (o >= n || nb >= n) throw std::runtime_error("assemble_pressure_correction: invalid ownership");

        const double ao = momentum_diagonal.component_data(0)[o];
        const double an = momentum_diagonal.component_data(0)[nb];
        if (!(ao > 0.0) || !(an > 0.0)) {
            throw std::runtime_error("assemble_pressure_correction: momentum diagonal must be positive");
        }

        const Vec3 d = cell_centre[nb] - cell_centre[o];
        const double distance = d.mag();
        if (!(distance > 1e-14)) {
            throw std::runtime_error("assemble_pressure_correction: degenerate cell-centre distance");
        }
        const double face_area_sq = sf[f].mag2();
        const double coefficient =
            0.5 * (1.0 / ao + 1.0 / an) * face_area_sq / distance;
        diagonal[o] += coefficient;
        diagonal[nb] += coefficient;
        rows[o].push_back({nb, -coefficient});
        rows[nb].push_back({o, -coefficient});
    }

    SparseMatrix matrix(n, n);
    for (std::size_t row = 0; row < n; ++row) {
        rows[row].push_back({row, diagonal[row]});
        std::sort(rows[row].begin(), rows[row].end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        for (const auto& [col, value] : rows[row]) matrix.push_back(row, col, value);
    }
    matrix.finalize();

    Vector rhs(n, 0.0);
    for (std::size_t c = 0; c < n; ++c)
        rhs(c) = -continuity.component_data(0)[c] * volume[c];

    return {std::move(matrix), std::move(rhs)};
}

}  // namespace cfdx::physics
