#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/gradient.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::physics {

enum class ScalarBoundaryType {
    FIXED_VALUE,
    ZERO_GRADIENT,
    FIXED_GRADIENT
};

struct ScalarBoundaryCondition {
    ScalarBoundaryType type = ScalarBoundaryType::ZERO_GRADIENT;
    double value = 0.0;
    double gradient = 0.0;
};

using ScalarBoundaryConditions = std::map<std::string, ScalarBoundaryCondition>;

enum class ConvectionScheme {
    UPWIND,
    SECOND_ORDER_UPWIND
};

struct ScalarBoundaryFaceValues {
    // Values are indexed by global face id. Only boundary faces need entries.
    std::map<std::string, std::vector<double>> values;
};

struct FvGeometry {
    std::vector<cfdx::core::Vec3> face_centres;
    std::vector<cfdx::core::Vec3> face_area_vectors;
    std::vector<cfdx::core::Vec3> cell_centres;
    std::vector<double> cell_volumes;
    std::vector<std::size_t> face_patch;
};

inline FvGeometry build_fv_geometry(const cfdx::core::Mesh& mesh)
{
    using namespace cfdx::core;
    const std::size_t nf = mesh.n_faces();
    const std::size_t nc = mesh.n_cells();

    FvGeometry g;
    g.face_centres.resize(nf);
    g.face_area_vectors.resize(nf);
    g.cell_centres.resize(nc);
    g.cell_volumes.resize(nc);
    g.face_patch.assign(nf, mesh.boundary().n_patches());

    for (std::size_t f = 0; f < nf; ++f) {
        const VertexIndex off = mesh.faces().offsets_data()[f];
        const VertexIndex count = mesh.faces().offsets_data()[f + 1] - off;
        const auto fg = compute_face_geometry(
            mesh.points().x_data(), mesh.points().y_data(), mesh.points().z_data(),
            mesh.faces().vertices_data(), off, count);
        g.face_centres[f] = fg.centre;
        g.face_area_vectors[f] = fg.Sf;
    }

    for (std::size_t c = 0; c < nc; ++c) {
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        const auto cg = compute_cell_geometry(
            g.face_centres.data(), g.face_area_vectors.data(),
            mesh.cells().faces_data() + off, count);
        g.cell_centres[c] = cg.centre;

        // Face area vectors are stored once, oriented from owner to neighbour.
        // For a neighbour cell the same face therefore has the opposite
        // outward normal. Cell geometry must account for that orientation;
        // otherwise internal faces corrupt neighbour-cell volumes.
        double signed_volume = 0.0;
        for (std::size_t k = 0; k < count; ++k) {
            const std::size_t f = mesh.cells().faces_data()[off + k];
            const auto owner = mesh.ownership().owner(f);
            const auto neighbour = mesh.ownership().neighbour(f);
            if (owner != c && neighbour != static_cast<int>(c))
                throw std::runtime_error("build_fv_geometry: cell-face ownership mismatch");
            const double orientation = owner == c ? 1.0 : -1.0;
            signed_volume +=
                (g.face_centres[f] - g.cell_centres[c]).dot(
                    g.face_area_vectors[f] * orientation);
        }
        signed_volume /= 3.0;
        g.cell_volumes[c] = std::abs(signed_volume);
        if (!(g.cell_volumes[c] > 0.0) || !std::isfinite(g.cell_volumes[c]))
            throw std::runtime_error("build_fv_geometry: non-positive cell volume");
    }

    for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p) {
        for (const auto f : mesh.boundary().patch(p).face_ids) {
            if (f >= nf || mesh.ownership().neighbour(f) >= 0)
                throw std::runtime_error("build_fv_geometry: invalid boundary face");
            g.face_patch[f] = p;
        }
    }

    return g;
}

struct ScalarEquation {
    cfdx::core::SparseMatrix matrix;
    cfdx::core::Vector rhs;
    std::vector<double> diagonal;
    double max_imbalance = 0.0;
};

struct ScalarSolveControls {
    std::size_t max_iterations = 1000;
    double tolerance = 1e-10;
    double relaxation = 1.0;
};

inline void validate_scalar_controls(const ScalarSolveControls& c)
{
    if (c.max_iterations == 0 || !(c.tolerance > 0.0) ||
        !(c.relaxation > 0.0 && c.relaxation <= 1.0))
        throw std::invalid_argument("invalid scalar solver controls");
}

// Assemble:
//     div(phi*psi) - div(gamma grad(psi)) = Su + Sp*psi
// with first-order upwind convection and two-point orthogonal diffusion.
// For steady problems, bounded convection subtracts div(phi)*psi, matching
// the standard bounded finite-volume treatment.
inline ScalarEquation assemble_scalar_equation(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& face_flux,
    double diffusion_coefficient,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& source_explicit,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& source_implicit,
    const ScalarBoundaryConditions& boundary_conditions = {},
    bool bounded_convection = true,
    const ScalarBoundaryFaceValues* face_values = nullptr,
    const std::vector<double>* extra_diagonal = nullptr,
    const std::vector<double>* extra_rhs = nullptr,
    const std::vector<double>* cell_diffusion = nullptr,
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* convected_field = nullptr)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    const std::size_t nf = mesh.n_faces();

    if (face_flux.size() != nf || face_flux.dimension() != 1 ||
        source_explicit.size() != nc || source_explicit.dimension() != 1 ||
        source_implicit.size() != nc || source_implicit.dimension() != 1 ||
        (extra_diagonal && extra_diagonal->size() != nc) ||
        (extra_rhs && extra_rhs->size() != nc) ||
        (cell_diffusion && cell_diffusion->size() != nc))
        throw std::invalid_argument("assemble_scalar_equation: field dimensions do not match mesh");
    if (!(diffusion_coefficient >= 0.0) || !std::isfinite(diffusion_coefficient))
        throw std::invalid_argument("assemble_scalar_equation: invalid diffusion coefficient");

    std::vector<std::map<std::size_t, double>> rows(nc);
    std::vector<double> rhs(nc, 0.0);
    std::vector<double> diag(nc, 0.0);
    std::vector<double> div_phi(nc, 0.0);
    std::vector<double> deferred_rhs(nc, 0.0);
    cfdx::core::Field<double, cfdx::core::Location::CELL> reconstructed_gradient;
    if (convection_scheme == ConvectionScheme::SECOND_ORDER_UPWIND) {
        if (convected_field == nullptr || convected_field->size() != nc || convected_field->dimension() != 1)
            throw std::invalid_argument("assemble_scalar_equation: second-order upwind requires a scalar convected field");
        reconstructed_gradient = cfdx::core::compute_gradient_gauss(*convected_field, mesh);
    }

    const auto& own = mesh.ownership();
    const double* phi = face_flux.component_data(0);

    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t o = own.owner(f);
        const auto nraw = own.neighbour(f);
        const double F = phi[f];
        if (!std::isfinite(F)) throw std::invalid_argument("assemble_scalar_equation: non-finite face flux");

        if (nraw >= 0) {
            const std::size_t n = static_cast<std::size_t>(nraw);
            const double distance = (geometry.cell_centres[n] - geometry.cell_centres[o]).mag();
            const double area = geometry.face_area_vectors[f].mag();
            if (!(distance > 0.0) || !(area > 0.0))
                throw std::runtime_error("assemble_scalar_equation: degenerate internal face");
            double gamma_face = diffusion_coefficient;
            if (cell_diffusion) {
                const double go = (*cell_diffusion)[o];
                const double gn = (*cell_diffusion)[n];
                if (go < 0.0 || gn < 0.0)
                    throw std::invalid_argument("assemble_scalar_equation: negative cell diffusion");
                gamma_face = (go > 0.0 && gn > 0.0)
                    ? 2.0 * go * gn / (go + gn)
                    : 0.0;
            }
            const double D = gamma_face * area / distance;

            const double a_on = D + std::max(F, 0.0);
            const double a_no = D + std::max(-F, 0.0);
            diag[o] += a_on;
            diag[n] += a_no;
            rows[o][n] -= a_no;
            rows[n][o] -= a_on;

            div_phi[o] += F;
            div_phi[n] -= F;

            if (convection_scheme == ConvectionScheme::SECOND_ORDER_UPWIND) {
                const std::size_t upwind = F >= 0.0 ? o : n;
                const double phi_up = (*convected_field)(upwind);
                const auto& C_up = geometry.cell_centres[upwind];
                const auto& Cf = geometry.face_centres[f];
                const double* gx = reconstructed_gradient.component_data(0);
                const double* gy = reconstructed_gradient.component_data(1);
                const double* gz = reconstructed_gradient.component_data(2);
                double phi_high = phi_up
                    + gx[upwind] * (Cf.x - C_up.x)
                    + gy[upwind] * (Cf.y - C_up.y)
                    + gz[upwind] * (Cf.z - C_up.z);
                const double phi_other = (*convected_field)(F >= 0.0 ? n : o);
                const double lo = std::min(phi_up, phi_other);
                const double hi = std::max(phi_up, phi_other);
                phi_high = std::clamp(phi_high, lo, hi);
                const double correction = F * (phi_high - phi_up);
                deferred_rhs[o] -= correction;
                deferred_rhs[n] += correction;
            }
        } else {
            const std::size_t patch = geometry.face_patch[f];
            ScalarBoundaryCondition bc;
            if (patch < mesh.boundary().n_patches()) {
                const auto& patch_name = mesh.boundary().patch(patch).name;
                const auto it = boundary_conditions.find(patch_name);
                if (it != boundary_conditions.end()) bc = it->second;
                if (face_values) {
                    const auto fv = face_values->values.find(patch_name);
                    if (fv != face_values->values.end() && f < fv->second.size()) {
                        bc.type = ScalarBoundaryType::FIXED_VALUE;
                        bc.value = fv->second[f];
                    }
                }
            }

            const double area = geometry.face_area_vectors[f].mag();
            const double distance = (geometry.face_centres[f] - geometry.cell_centres[o]).mag();
            if (!(distance > 0.0) || !(area > 0.0))
                throw std::runtime_error("assemble_scalar_equation: degenerate boundary face");

            if (bc.type == ScalarBoundaryType::FIXED_VALUE) {
                const double gamma_owner = cell_diffusion
                    ? (*cell_diffusion)[o] : diffusion_coefficient;
                if (gamma_owner < 0.0)
                    throw std::invalid_argument("assemble_scalar_equation: negative boundary diffusion");
                const double D = gamma_owner * area / distance;
                const double a_owner = D + std::max(F, 0.0);
                const double a_boundary = D + std::max(-F, 0.0);
                diag[o] += a_owner;
                rhs[o] += a_boundary * bc.value;
                div_phi[o] += F;
            } else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT) {
                diag[o] += std::max(F, 0.0);
                rhs[o] += (F < 0.0 ? -F * bc.value : 0.0);
                rhs[o] += diffusion_coefficient * area * bc.gradient;
                div_phi[o] += F;
            } else {
                // Zero-gradient upwind boundary: outflow contributes to the
                // owner diagonal. For inflow, the external state is supplied
                // by the boundary condition and must not make the diagonal
                // negative.
                diag[o] += std::max(F, 0.0);
                div_phi[o] += F;
            }
        }
    }

    ScalarEquation eq{SparseMatrix(nc, nc), Vector(nc, 0.0), diag, 0.0};
    for (std::size_t c = 0; c < nc; ++c) {
        // Bounded steady convection: subtract div(phi)*psi from the
        // discretized convection operator. The term vanishes when continuity
        // is satisfied, while preserving constants during intermediate SIMPLE
        // iterations.
        if (bounded_convection) diag[c] -= div_phi[c];

        diag[c] -= source_implicit(c) * geometry.cell_volumes[c];
        if (extra_diagonal) diag[c] += (*extra_diagonal)[c];
        if (!(diag[c] > 0.0) || !std::isfinite(diag[c])) {
            throw std::runtime_error(
                "assemble_scalar_equation: non-positive diagonal at cell " +
                std::to_string(c) + " diag=" + std::to_string(diag[c]) +
                " div_phi=" + std::to_string(div_phi[c]) +
                " Sp=" + std::to_string(source_implicit(c)) +
                " V=" + std::to_string(geometry.cell_volumes[c]));
        }

        rhs[c] += source_explicit(c) * geometry.cell_volumes[c];
        rhs[c] += deferred_rhs[c];
        if (extra_rhs) rhs[c] += (*extra_rhs)[c];

        std::vector<std::pair<std::size_t, double>> entries;
        entries.reserve(rows[c].size() + 1);
        for (const auto& [col, value] : rows[c]) entries.push_back({col, value});
        entries.push_back({c, diag[c]});
        std::sort(entries.begin(), entries.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        std::vector<std::pair<std::size_t, double>> merged;
        for (const auto& e : entries) {
            if (!merged.empty() && merged.back().first == e.first)
                merged.back().second += e.second;
            else
                merged.push_back(e);
        }
        for (const auto& [col, value] : merged) eq.matrix.push_back(c, col, value);
        eq.rhs(c) = rhs[c];
        eq.diagonal[c] = diag[c];
    }
    eq.matrix.finalize();

    double total_flux = 0.0;
    for (std::size_t c = 0; c < nc; ++c) total_flux += div_phi[c];
    eq.max_imbalance = std::abs(total_flux);
    return eq;
}

inline cfdx::core::SolverResult solve_scalar_equation(
    const ScalarEquation& equation,
    cfdx::core::Vector& solution,
    const ScalarSolveControls& controls = {})
{
    validate_scalar_controls(controls);
    if (solution.size() != equation.rhs.size())
        throw std::invalid_argument("solve_scalar_equation: solution size mismatch");

    if (solution.size() == 1) {
        const auto begin = equation.matrix.row_offsets_data()[0];
        const auto end = equation.matrix.row_offsets_data()[1];
        double diagonal = 0.0;
        for (std::uint32_t k = begin; k < end; ++k) {
            if (equation.matrix.columns_data()[k] == 0)
                diagonal += equation.matrix.values_data()[k];
        }
        if (!(std::abs(diagonal) > 0.0) || !std::isfinite(diagonal))
            throw std::runtime_error("solve_scalar_equation: singular 1x1 system");
        const double candidate_value = equation.rhs(0) / diagonal;
        const double residual = std::abs(diagonal * candidate_value - equation.rhs(0));
        solution(0) += controls.relaxation * (candidate_value - solution(0));
        return {
            residual <= controls.tolerance
                ? cfdx::core::SolverStatus::CONVERGED
                : cfdx::core::SolverStatus::MAX_ITER_REACHED,
            1, residual, residual
        };
    }

    cfdx::core::Vector candidate = solution;
    auto result = cfdx::core::solve_bicgstab(
        equation.matrix, equation.rhs, candidate,
        controls.max_iterations, controls.tolerance);

    if (result.status == cfdx::core::SolverStatus::CONVERGED ||
        result.status == cfdx::core::SolverStatus::MAX_ITER_REACHED) {
        for (std::size_t i = 0; i < candidate.size(); ++i) {
            if (!std::isfinite(candidate(i))) {
                result.status = cfdx::core::SolverStatus::DIVERGED;
                return result;
            }
        }
        if (!std::isfinite(result.residual) || !std::isfinite(result.residual_relative)) {
            result.status = cfdx::core::SolverStatus::DIVERGED;
            return result;
        }
        for (std::size_t i = 0; i < solution.size(); ++i)
            solution(i) += controls.relaxation * (candidate(i) - solution(i));
    }
    return result;
}

inline double scalar_equation_residual_inf(
    const ScalarEquation& equation,
    const cfdx::core::Vector& solution)
{
    if (solution.size() != equation.rhs.size())
        throw std::invalid_argument("scalar_equation_residual_inf: dimension mismatch");
    double max_r = 0.0;
    for (std::size_t i = 0; i < solution.size(); ++i) {
        double r = -equation.rhs(i);
        const auto row_begin = equation.matrix.row_offsets_data()[i];
        const auto row_end = equation.matrix.row_offsets_data()[i + 1];
        for (std::uint32_t k = row_begin; k < row_end; ++k)
            r += equation.matrix.values_data()[k] * solution(equation.matrix.columns_data()[k]);
        max_r = std::max(max_r, std::abs(r));
    }
    return max_r;
}

} // namespace cfdx::physics
