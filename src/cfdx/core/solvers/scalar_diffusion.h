#pragma once

#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/physics/poisson/poisson.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::core {

/// Backward-compatible name for the canonical Poisson boundary contract.
using DirichletBoundary = PoissonBoundaryCondition;

struct ScalarDiffusionConfig {
    double diffusivity = 1.0;
    std::size_t max_iterations = 2000;
    double tolerance = 1e-10;
    PoissonGauge gauge{};
    double compatibility_tolerance = 1e-12;
};

struct ScalarDiffusionResult {
    SolverResult linear_result;
    Vector solution;
    SparseMatrix matrix;
    Vector rhs;
};

inline SparseMatrix assemble_cell_diffusion_matrix(
    const Mesh& mesh,
    const DirichletBoundary& boundary,
    double diffusivity,
    const std::vector<double>& source,
    Vector& rhs,
    const GeometryCache& geometry)
{
    if (!(diffusivity > 0.0) || !std::isfinite(diffusivity))
        throw std::invalid_argument(
            "assemble_cell_diffusion_matrix: diffusivity must be positive");
    if (source.size() != mesh.n_cells())
        throw std::invalid_argument(
            "assemble_cell_diffusion_matrix: source size mismatch");
    boundary.validate(mesh.n_faces());
    for (double value : source) {
        if (!std::isfinite(value))
            throw std::invalid_argument(
                "assemble_cell_diffusion_matrix: source must be finite");
    }

    const std::size_t nc = mesh.n_cells();
    const std::size_t nf = mesh.n_faces();
    if (!is_valid(geometry, mesh))
        throw std::invalid_argument(
            "assemble_cell_diffusion_matrix: invalid geometry cache");

    // Canonical PDE convention:
    //
    //     -div(Gamma grad(phi)) = S
    //
    // with S a volumetric source density. The finite-volume balance therefore
    // contributes S_c * V_c to the integrated right-hand side.
    rhs = Vector(nc);
    for (std::size_t c = 0; c < nc; ++c)
        rhs(c) = source[c] * geometry.cell_volumes[c];

    // Assemble row-wise in deterministic face order. For each internal face,
    // use the two-point finite-volume conductance Gamma A / d. For Dirichlet
    // boundaries, eliminate the boundary unknown into the RHS.
    std::vector<std::vector<std::pair<std::size_t, double>>> rows(nc);
    for (std::size_t c = 0; c < nc; ++c)
        rows[c].push_back({c, 0.0});

    const auto& own = mesh.ownership();
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t o = own.owner(f);
        if (o >= nc)
            throw std::runtime_error(
                "assemble_cell_diffusion_matrix: invalid owner");

        const auto n = own.neighbour(f);
        const double area = geometry.face_Sf[f].mag();
        if (!(area > 0.0) || !std::isfinite(area))
            throw std::runtime_error(
                "assemble_cell_diffusion_matrix: invalid face area");

        if (n >= 0) {
            const std::size_t j = static_cast<std::size_t>(n);
            if (j >= nc)
                throw std::runtime_error(
                    "assemble_cell_diffusion_matrix: invalid neighbour");
            const double d =
                (geometry.cell_centres[j] - geometry.cell_centres[o]).mag();
            if (!(d > std::numeric_limits<double>::epsilon()) ||
                !std::isfinite(d))
                throw std::runtime_error(
                    "assemble_cell_diffusion_matrix: invalid cell-centre "
                    "distance");

            const double g = diffusivity * area / d;
            rows[o][0].second += g;
            rows[j][0].second += g;
            rows[o].push_back({j, -g});
            rows[j].push_back({o, -g});
        } else if (std::isfinite(boundary.face_values[f])) {
            const double d =
                (geometry.face_centres[f] - geometry.cell_centres[o]).mag();
            if (!(d > std::numeric_limits<double>::epsilon()) ||
                !std::isfinite(d))
                throw std::runtime_error(
                    "assemble_cell_diffusion_matrix: invalid boundary "
                    "distance");
            const double g = diffusivity * area / d;
            if (boundary.type_for_face(f) == PoissonBoundaryType::DIRICHLET) {
                rows[o][0].second += g;
                rhs(o) += g * boundary.face_values[f];
            } else {
                // q_n is the prescribed outward flux
                // Gamma * grad(phi) . n_out. Since
                // -div(Gamma grad(phi)) = S, integration gives
                // -sum_boundary q_n A = S V. Moving a prescribed
                // boundary flux to the algebraic RHS therefore gives
                // A * phi = S V + q_n A.
                rhs(o) += boundary.face_values[f] * area;
            }
        }
        // A NaN boundary value is deliberately an unprescribed face.
    }

    // Consolidate duplicates and sort columns for stable CSR.
    SparseMatrix A(nc, nc);
    for (std::size_t i = 0; i < nc; ++i) {
        std::sort(rows[i].begin(), rows[i].end(),
                  [](const auto& a, const auto& b) {
                      return a.first < b.first;
                  });
        for (std::size_t k = 0; k < rows[i].size();) {
            const std::size_t col = rows[i][k].first;
            double value = 0.0;
            while (k < rows[i].size() && rows[i][k].first == col) {
                value += rows[i][k].second;
                ++k;
            }
            if (value != 0.0 || col == i)
                A.push_back(i, col, value);
        }
    }

    A.finalize();
    if (!A.is_consistent())
        throw std::runtime_error(
            "assemble_cell_diffusion_matrix: invalid CSR");
    return A;
}

inline SparseMatrix assemble_cell_diffusion_matrix(
    const Mesh& mesh,
    const DirichletBoundary& boundary,
    double diffusivity,
    const std::vector<double>& source,
    Vector& rhs)
{
    const GeometryCache geometry = make_geometry_cache(mesh);
    return assemble_cell_diffusion_matrix(
        mesh, boundary, diffusivity, source, rhs, geometry);
}

inline SparseMatrix apply_reference_cell_gauge(
    const SparseMatrix& matrix,
    Vector& rhs,
    std::size_t reference_cell,
    double reference_value)
{
    const std::size_t n = matrix.n_rows();
    if (matrix.n_cols() != n)
        throw std::invalid_argument("apply_reference_cell_gauge: matrix must be square");
    if (rhs.size() != n)
        throw std::invalid_argument("apply_reference_cell_gauge: RHS size mismatch");
    if (reference_cell >= n)
        throw std::invalid_argument("apply_reference_cell_gauge: reference cell out of range");
    if (!std::isfinite(reference_value))
        throw std::invalid_argument("apply_reference_cell_gauge: reference value must be finite");

    SparseMatrix gauged(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        if (i != reference_cell)
            rhs(i) -= matrix(i, reference_cell) * reference_value;
        const auto begin = matrix.row_offsets_data()[i];
        const auto end = matrix.row_offsets_data()[i + 1];
        for (std::uint32_t k = begin; k < end; ++k) {
            const std::size_t j = matrix.columns_data()[k];
            if (i == reference_cell || j == reference_cell)
                continue;
            gauged.push_back(i, j, matrix.values_data()[k]);
        }
    }
    gauged.push_back(reference_cell, reference_cell, 1.0);
    rhs(reference_cell) = reference_value;
    gauged.finalize();
    return gauged;
}

inline bool has_dirichlet_boundary(
    const Mesh& mesh, const PoissonBoundaryCondition& boundary)
{
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (std::isfinite(boundary.face_values[f]) &&
            boundary.type_for_face(f) == PoissonBoundaryType::DIRICHLET)
            return true;
    }
    return false;
}

inline double pure_neumann_compatibility_residual(
    const Mesh& mesh,
    const PoissonBoundaryCondition& boundary,
    const std::vector<double>& source,
    const GeometryCache& geometry)
{
    double balance = 0.0;
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        balance += source[c] * geometry.cell_volumes[c];
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (std::isfinite(boundary.face_values[f]) &&
            boundary.type_for_face(f) == PoissonBoundaryType::NEUMANN)
            balance += boundary.face_values[f] * geometry.face_Sf[f].mag();
    }
    return balance;
}

inline ScalarDiffusionResult solve_poisson_mixed(
    const Mesh& mesh,
    const PoissonBoundaryCondition& boundary,
    const std::vector<double>& source,
    const ScalarDiffusionConfig& config = {})
{
    const GeometryCache geometry = make_geometry_cache(mesh);
    ScalarDiffusionResult out;
    out.solution = Vector(mesh.n_cells(), 0.0);
    out.matrix = assemble_cell_diffusion_matrix(
        mesh, boundary, config.diffusivity, source, out.rhs, geometry);

    if (!has_dirichlet_boundary(mesh, boundary)) {
        const double balance =
            pure_neumann_compatibility_residual(mesh, boundary, source, geometry);
        double scale = 1.0;
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
            scale += std::abs(source[cell]) * geometry.cell_volumes[cell];
        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
            if (std::isfinite(boundary.face_values[f]) &&
                boundary.type_for_face(f) == PoissonBoundaryType::NEUMANN)
                scale += std::abs(boundary.face_values[f]) *
                         geometry.face_Sf[f].mag();
        }
        if (!std::isfinite(balance) ||
            !std::isfinite(config.compatibility_tolerance) ||
            config.compatibility_tolerance < 0.0 ||
            std::abs(balance) > config.compatibility_tolerance * scale) {
            throw std::invalid_argument(
                "solve_poisson_mixed: incompatible pure-Neumann data; "
                "integral(S dV) + integral(q_n dA) must be zero");
        }
        if (config.gauge.type != PoissonGaugeType::REFERENCE_CELL)
            throw std::invalid_argument(
                "solve_poisson_mixed: unsupported Poisson gauge");
        out.matrix = apply_reference_cell_gauge(
            out.matrix, out.rhs, config.gauge.reference_cell,
            config.gauge.reference_value);
        out.solution(config.gauge.reference_cell) =
            config.gauge.reference_value;
    }

    out.linear_result = solve_cg(
        out.matrix, out.rhs, out.solution,
        config.max_iterations, config.tolerance);
    return out;
}

inline ScalarDiffusionResult solve_poisson_dirichlet(
    const Mesh& mesh,
    const DirichletBoundary& boundary,
    const std::vector<double>& source,
    const ScalarDiffusionConfig& config = {})
{
    boundary.validate(mesh.n_faces());
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (std::isfinite(boundary.face_values[f]) &&
            boundary.type_for_face(f) != PoissonBoundaryType::DIRICHLET)
            throw std::invalid_argument(
                "solve_poisson_dirichlet: all prescribed boundary faces "
                "must be Dirichlet");
    }
    return solve_poisson_mixed(mesh, boundary, source, config);
}

inline ScalarDiffusionResult solve_laplace_dirichlet(
    const Mesh& mesh,
    const DirichletBoundary& boundary,
    const ScalarDiffusionConfig& config = {})
{
    return solve_poisson_dirichlet(
        mesh, boundary,
        std::vector<double>(mesh.n_cells(), 0.0), config);
}

} // namespace cfdx::core
