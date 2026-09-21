#pragma once

#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/mesh/mesh.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
#include <utility>

namespace cfdx::core {

struct DirichletBoundary {
    // NaN means that the face is not prescribed. Boundary faces without a
    // value are treated as zero-normal-gradient (no contribution).
    std::vector<double> face_values;
};

struct ScalarDiffusionConfig {
    double diffusivity = 1.0;
    std::size_t max_iterations = 2000;
    double tolerance = 1e-10;
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
    Vector& rhs)
{
    if (!(diffusivity > 0.0) || !std::isfinite(diffusivity))
        throw std::invalid_argument("assemble_cell_diffusion_matrix: diffusivity must be positive");
    if (source.size() != mesh.n_cells())
        throw std::invalid_argument("assemble_cell_diffusion_matrix: source size mismatch");
    if (!boundary.face_values.empty() && boundary.face_values.size() != mesh.n_faces())
        throw std::invalid_argument("assemble_cell_diffusion_matrix: boundary size mismatch");

    const std::size_t nc = mesh.n_cells();
    const std::size_t nf = mesh.n_faces();
    rhs = Vector(nc);
    for (std::size_t c = 0; c < nc; ++c) rhs(c) = source[c];

    std::vector<Vec3> fc(nf);
    std::vector<Vec3> sf(nf);
    for (std::size_t f = 0; f < nf; ++f) {
        const auto off = mesh.faces().offsets_data()[f];
        const auto n = mesh.faces().offsets_data()[f + 1] - off;
        auto g = compute_face_geometry(
            mesh.points().x_data(), mesh.points().y_data(), mesh.points().z_data(),
            mesh.faces().vertices_data(), off, n);
        fc[f] = g.centre;
        sf[f] = g.Sf;
    }

    std::vector<Vec3> cc(nc);
    for (std::size_t c = 0; c < nc; ++c) {
        const auto off = mesh.cells().offsets_data()[c];
        const auto n = mesh.cells().offsets_data()[c + 1] - off;
        cc[c] = compute_cell_geometry(fc.data(), sf.data(),
                                      mesh.cells().faces_data() + off, n).centre;
    }

    // Assemble row-wise in deterministic face order. For each internal face,
    // use the two-point finite-volume conductance k A / d. For Dirichlet
    // boundaries, eliminate the boundary unknown into the RHS.
    std::vector<std::vector<std::pair<std::size_t, double>>> rows(nc);
    for (std::size_t c = 0; c < nc; ++c) rows[c].push_back({c, 0.0});

    const auto& own = mesh.ownership();
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t o = own.owner(f);
        if (o >= nc) throw std::runtime_error("assemble_cell_diffusion_matrix: invalid owner");
        const auto n = own.neighbour(f);
        const double area = sf[f].mag();
        if (!(area > 0.0)) throw std::runtime_error("assemble_cell_diffusion_matrix: zero face area");

        if (n >= 0) {
            const std::size_t j = static_cast<std::size_t>(n);
            if (j >= nc) throw std::runtime_error("assemble_cell_diffusion_matrix: invalid neighbour");
            const double d = (cc[j] - cc[o]).mag();
            if (!(d > std::numeric_limits<double>::epsilon()))
                throw std::runtime_error("assemble_cell_diffusion_matrix: zero cell-centre distance");
            const double g = diffusivity * area / d;
            rows[o][0].second += g;
            rows[j][0].second += g;
            rows[o].push_back({j, -g});
            rows[j].push_back({o, -g});
        } else if (!boundary.face_values.empty() &&
                   std::isfinite(boundary.face_values[f])) {
            const double d = (fc[f] - cc[o]).mag();
            if (!(d > std::numeric_limits<double>::epsilon()))
                throw std::runtime_error("assemble_cell_diffusion_matrix: zero boundary distance");
            const double g = diffusivity * area / d;
            rows[o][0].second += g;
            rhs(o) += g * boundary.face_values[f];
        }
    }

    // Consolidate duplicates (a polyhedral cell can be connected to another
    // cell through more than one face) and sort columns for stable CSR.
    SparseMatrix A(nc, nc);
    for (std::size_t i = 0; i < nc; ++i) {
        std::sort(rows[i].begin(), rows[i].end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        for (std::size_t k = 0; k < rows[i].size();) {
            const std::size_t col = rows[i][k].first;
            double value = 0.0;
            while (k < rows[i].size() && rows[i][k].first == col) {
                value += rows[i][k].second;
                ++k;
            }
            if (value != 0.0 || col == i) A.push_back(i, col, value);
        }
    }
    A.finalize();
    if (!A.is_consistent()) throw std::runtime_error("assemble_cell_diffusion_matrix: invalid CSR");
    return A;
}

inline ScalarDiffusionResult solve_poisson_dirichlet(
    const Mesh& mesh,
    const DirichletBoundary& boundary,
    const std::vector<double>& source,
    const ScalarDiffusionConfig& config = {})
{
    ScalarDiffusionResult out;
    out.solution = Vector(mesh.n_cells(), 0.0);
    out.matrix = assemble_cell_diffusion_matrix(
        mesh, boundary, config.diffusivity, source, out.rhs);
    out.linear_result = solve_cg(
        out.matrix, out.rhs, out.solution,
        config.max_iterations, config.tolerance);
    return out;
}

inline ScalarDiffusionResult solve_laplace_dirichlet(
    const Mesh& mesh,
    const DirichletBoundary& boundary,
    const ScalarDiffusionConfig& config = {})
{
    return solve_poisson_dirichlet(mesh, boundary,
                                   std::vector<double>(mesh.n_cells(), 0.0),
                                   config);
}

} // namespace cfdx::core
