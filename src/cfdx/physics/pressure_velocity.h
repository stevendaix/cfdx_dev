#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/physics/finite_volume_transport.h"
#include <cmath>
#include <cstddef>
#include <algorithm>
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
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& continuity,
    double density = 1.0)
{
    using namespace cfdx::core;
    if (!(density > 0.0) || !std::isfinite(density))
        throw std::invalid_argument(
            "assemble_pressure_correction: density must be finite and positive");
    const std::size_t n = mesh.n_cells();
    if (momentum_diagonal.dimension() != 1 || momentum_diagonal.size() != n ||
        continuity.dimension() != 1 || continuity.size() != n) {
        throw std::runtime_error("assemble_pressure_correction: invalid fields");
    }

    // Use the canonical FV geometry builder so imported face winding is
    // normalized before cell volumes and centres are computed. Reconstructing
    // this geometry locally used the raw face winding and produced a different
    // pressure coefficient on meshes whose face orientation was not canonical.
    const FvGeometry geometry = build_fv_geometry(mesh);

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

        const Vec3 d = geometry.cell_centres[nb] - geometry.cell_centres[o];
        const double distance = d.mag();
        if (!(distance > 1e-14)) {
            throw std::runtime_error("assemble_pressure_correction: degenerate cell-centre distance");
        }
        const double face_area = geometry.face_area_vectors[f].mag();
        if (!(face_area > 0.0) || !std::isfinite(face_area))
            throw std::runtime_error(
                "assemble_pressure_correction: degenerate face area");
        const Vec3 e{d.x / distance, d.y / distance, d.z / distance};
        const Vec3 Sf = geometry.face_area_vectors[f];
        const double orthogonal_area = Sf.x * e.x +
                                       Sf.y * e.y +
                                       Sf.z * e.z;
        // This legacy scalar helper follows the production convention
        // dAU = V/A_P. The pressure-flux derivative is dAU_f*(Sf·e)/d;
        // using (1/A_P)_f*|Sf|^2/d omits the cell volume and assumes
        // orthogonality, so it disagrees with the production operator on
        // non-unit or skewed cells.
        const double dAU_face = 0.5 * (
            geometry.cell_volumes[o] / ao + geometry.cell_volumes[nb] / an);
        const double coefficient = density * dAU_face * orthogonal_area / distance;
        if (!(coefficient > 0.0) || !std::isfinite(coefficient))
            throw std::runtime_error(
                "assemble_pressure_correction: invalid face pressure coefficient");
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
    // Continuity is already a cell-integrated mass/volume flux imbalance.
    // Do not multiply it by cell volume a second time: that would destroy
    // the global conservation constraint on non-uniform meshes.
    for (std::size_t c = 0; c < n; ++c)
        rhs(c) = -continuity.component_data(0)[c];

    return {std::move(matrix), std::move(rhs)};
}

}  // namespace cfdx::physics
