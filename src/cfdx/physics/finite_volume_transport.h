#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/boundary_constraint_fvm.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/gauss_seidel_solver.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/numerics/interpolation.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>
#include <limits>

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
    SECOND_ORDER_UPWIND,
    // Monotone TVD reconstruction using a minmod limiter on the MUSCL slope.
    TVD,
    CENTRAL,
    BLENDED,
    QUICK,
    QUICK_BOUNDED
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
        const auto cg = compute_cell_geometry_oriented(
            g.face_centres.data(), g.face_area_vectors.data(),
            mesh.cells().faces_data() + off, count, static_cast<CellIndex>(c),
            mesh.ownership());
        g.cell_centres[c] = cg.centre;
        g.cell_volumes[c] = cg.volume;
        if (!(cg.volume > 0.0) || !std::isfinite(cg.volume))
            throw std::runtime_error("build_fv_geometry: non-positive cell volume");
    }

    // Physics kernels use face_area_vectors with the invariant convention
    // Sf · (C_neighbour - C_owner) > 0 for internal faces and
    // Sf · (Cf - C_owner) > 0 for boundary faces. The raw face vertex order
    // is not guaranteed to encode this owner orientation, so normalize the
    // global face vectors after cell centres are available.
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t o = mesh.ownership().owner(f);
        if (o >= nc)
            throw std::runtime_error("build_fv_geometry: invalid face owner");
        const auto nraw = mesh.ownership().neighbour(f);
        if (nraw >= 0) {
            const std::size_t n = static_cast<std::size_t>(nraw);
            if (n >= nc)
                throw std::runtime_error("build_fv_geometry: invalid face neighbour");
            ensure_face_orientation(
                g.face_area_vectors[f], g.face_centres[f],
                g.cell_centres[o], &g.cell_centres[n]);
        } else {
            ensure_face_orientation(
                g.face_area_vectors[f], g.face_centres[f],
                g.cell_centres[o], nullptr);
        }
    }

    // Recompute cell geometry after the face vectors have been normalized. This
    // second pass is required because the first pass may have used an imported
    // face winding that was inward with respect to the owner.
    for (std::size_t c = 0; c < nc; ++c) {
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        const auto cg = compute_cell_geometry_oriented(
            g.face_centres.data(), g.face_area_vectors.data(),
            mesh.cells().faces_data() + off, count, static_cast<CellIndex>(c),
            mesh.ownership());
        g.cell_centres[c] = cg.centre;
        g.cell_volumes[c] = cg.volume;
        if (!(cg.volume > 0.0) || !std::isfinite(cg.volume))
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

// Method contract: docs/audits/01-diffusion-laplacien.md
struct FluxImbalanceMetrics {
    double signed_global = 0.0;
    double global_l1 = 0.0;
    double local_linf = 0.0;
};

inline FluxImbalanceMetrics flux_imbalance_metrics(
    const std::vector<double>& cell_flux_residuals)
{
    FluxImbalanceMetrics metrics;
    for (const double residual : cell_flux_residuals) {
        if (!std::isfinite(residual))
            throw std::invalid_argument("flux_imbalance_metrics: non-finite residual");
        metrics.signed_global += residual;
        metrics.global_l1 += std::abs(residual);
        metrics.local_linf = std::max(metrics.local_linf, std::abs(residual));
    }
    return metrics;
}

struct ScalarEquation {
    cfdx::core::SparseMatrix matrix;
    cfdx::core::Vector rhs;
    std::vector<double> diagonal;
    // Compatibility alias for the maximum cell-local face-flux imbalance.
    double max_imbalance = 0.0;
    double signed_global_imbalance = 0.0;
    double global_l1_imbalance = 0.0;
    double max_local_imbalance = 0.0;
    // Assembly diagnostics retained temporarily for Issue #377 root-cause analysis.
    double max_abs_face_flux = 0.0;
    double min_internal_distance = std::numeric_limits<double>::infinity();
    double max_internal_distance = 0.0;
    double min_face_area = std::numeric_limits<double>::infinity();
    double max_face_area = 0.0;
};

struct ScalarSolveControls {
    std::size_t max_iterations = 1000;
    double tolerance = 1e-10;
    double relaxation = 1.0;
    cfdx::core::LinearProblemKind problem =
        cfdx::core::LinearProblemKind::ScalarTransport;
    cfdx::core::LinearSolverRequest linear_solver;

    ScalarSolveControls() = default;
    ScalarSolveControls(
        std::size_t iterations,
        double solve_tolerance,
        double solve_relaxation,
        cfdx::core::LinearProblemKind problem_kind =
            cfdx::core::LinearProblemKind::ScalarTransport,
        cfdx::core::LinearSolverRequest solver = {})
        : max_iterations(iterations),
          tolerance(solve_tolerance),
          relaxation(solve_relaxation),
          problem(problem_kind),
          linear_solver(solver) {}
};

inline void validate_scalar_controls(const ScalarSolveControls& c)
{
    if (c.max_iterations == 0 || !(c.tolerance > 0.0) ||
        !(c.relaxation > 0.0 && c.relaxation <= 1.0))
        throw std::invalid_argument("invalid scalar solver controls");
}

inline double boundary_normal_distance(
    const FvGeometry& geometry, std::size_t face, std::size_t cell)
{
    const auto& Sf = geometry.face_area_vectors[face];
    const double area = Sf.mag();
    if (!(area > 0.0) || !std::isfinite(area))
        throw std::runtime_error("boundary_normal_distance: degenerate face");
    const auto d = geometry.face_centres[face] - geometry.cell_centres[cell];
    const double dn = std::abs(d.dot(Sf)) / area;
    if (!(dn > 0.0) || !std::isfinite(dn)) {
        // 2D fallback: face centre coincides with cell centre (degenerate topology
        // where 2D cells are represented as single faces). Approximate the normal
        // distance using the cell's characteristic length sqrt(area).
        const double cell_vol = geometry.cell_volumes[cell];
        if (cell_vol > 0.0 && std::isfinite(cell_vol)) {
            const double dn2 = std::sqrt(cell_vol);
            if (dn2 > 0.0 && std::isfinite(dn2))
                return dn2;
        }
        throw std::runtime_error("boundary_normal_distance: invalid normal distance");
    }
    return dn;
}

inline ScalarBoundaryCondition resolve_scalar_boundary_condition(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    std::size_t face,
    const ScalarBoundaryConditions& boundary_conditions,
    const ScalarBoundaryFaceValues* face_values = nullptr,
    const ScalarBoundaryFaceConditions* face_conditions = nullptr)
{
    if (face >= mesh.n_faces() || face >= geometry.face_patch.size())
        throw std::out_of_range("resolve_scalar_boundary_condition: invalid face");
    ScalarBoundaryCondition bc;
    if (face_conditions && face_conditions->has(face)) {
        const auto& resolved = face_conditions->conditions[face];
        if (!std::isfinite(resolved.value))
            throw std::invalid_argument(
                "resolve_scalar_boundary_condition: non-finite face condition");
        if (resolved.type == ScalarBoundaryFaceCondition::Type::FIXED_VALUE) {
            bc.type = ScalarBoundaryType::FIXED_VALUE;
            bc.value = resolved.value;
        } else {
            bc.type = ScalarBoundaryType::FIXED_GRADIENT;
            bc.gradient = resolved.value;
        }
        return bc;
    }

    const std::size_t patch = geometry.face_patch[face];
    if (patch < mesh.boundary().n_patches()) {
        const auto& patch_name = mesh.boundary().patch(patch).name;
        const auto it = boundary_conditions.find(patch_name);
        if (it != boundary_conditions.end()) bc = it->second;
        if (face_values) {
            const auto fv = face_values->values.find(patch_name);
            if (fv != face_values->values.end() &&
                face < fv->second.size() && std::isfinite(fv->second[face])) {
                bc.type = ScalarBoundaryType::FIXED_VALUE;
                bc.value = fv->second[face];
            }
        }
    }
    if (bc.type == ScalarBoundaryType::FIXED_VALUE && !std::isfinite(bc.value))
        throw std::invalid_argument(
            "resolve_scalar_boundary_condition: non-finite fixed value");
    if (bc.type == ScalarBoundaryType::FIXED_GRADIENT && !std::isfinite(bc.gradient))
        throw std::invalid_argument(
            "resolve_scalar_boundary_condition: non-finite fixed gradient");
    return bc;
}

// Boundary-aware scalar gradient shared by transport and diffusion-flux
// reconstruction. Face-specific values/conditions take precedence over patch
// conditions so the gradient uses the same boundary state as the assembly.
inline cfdx::core::Field<double, cfdx::core::Location::CELL>
compute_scalar_gradient_with_boundary(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& field,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarBoundaryConditions& boundary_conditions,
    const ScalarBoundaryFaceValues* face_values = nullptr,
    const ScalarBoundaryFaceConditions* face_conditions = nullptr)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    const std::size_t nf = mesh.n_faces();
    if (field.dimension() != 1 || field.size() != nc ||
        geometry.cell_centres.size() != nc ||
        geometry.cell_volumes.size() != nc ||
        geometry.face_centres.size() != nf ||
        geometry.face_area_vectors.size() != nf ||
        geometry.face_patch.size() != nf ||
        (face_conditions &&
         (face_conditions->conditions.size() != nf ||
          face_conditions->specified.size() != nf)))
        throw std::invalid_argument(
            "compute_scalar_gradient_with_boundary: field/geometry shape mismatch");

    Field<double, Location::CELL> gradient(
        nc, field.name() + "_gradient_bc", "", 3);
    const auto& ownership = mesh.ownership();
    for (std::size_t cell = 0; cell < nc; ++cell) {
        Vec3 sum{0.0, 0.0, 0.0};
        const Offset off = mesh.cells().offsets_data()[cell];
        const Offset count = mesh.cells().offsets_data()[cell + 1] - off;
        for (Offset k = 0; k < count; ++k) {
            const std::size_t face = mesh.cells().faces_data()[off + k];
            const bool is_owner = ownership.owner(face) == cell;
            const Vec3 sf = is_owner
                ? geometry.face_area_vectors[face]
                : geometry.face_area_vectors[face] * -1.0;
            double face_value = field(cell);
            const auto neighbour_raw = ownership.neighbour(face);
            if (neighbour_raw >= 0) {
                const std::size_t other = is_owner
                    ? static_cast<std::size_t>(neighbour_raw)
                    : ownership.owner(face);
                if (other >= nc)
                    throw std::runtime_error(
                        "compute_scalar_gradient_with_boundary: invalid neighbour");
                const double dc =
                    (geometry.face_centres[face] - geometry.cell_centres[cell]).mag();
                const double dn =
                    (geometry.cell_centres[other] - geometry.face_centres[face]).mag();
                if (!(dc + dn > 0.0) || !std::isfinite(dc + dn))
                    throw std::runtime_error(
                        "compute_scalar_gradient_with_boundary: degenerate face distance");
                const double w = dn / (dc + dn);
                face_value = w * field(cell) + (1.0 - w) * field(other);
            } else {
                const std::size_t patch = geometry.face_patch[face];
                if (patch < mesh.boundary().n_patches() &&
                    mesh.boundary().patch(patch).type == PatchType::EMPTY)
                    continue;

                const auto bc = resolve_scalar_boundary_condition(
                    mesh, geometry, face, boundary_conditions,
                    face_values, face_conditions);
                if (bc.type == ScalarBoundaryType::FIXED_VALUE) {
                    face_value = bc.value;
                } else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT) {
                    face_value = field(cell) + bc.gradient *
                        boundary_normal_distance(geometry, face, cell);
                }
            }
            if (!std::isfinite(face_value))
                throw std::runtime_error(
                    "compute_scalar_gradient_with_boundary: non-finite face value");
            sum += sf * face_value;
        }
        const double volume = geometry.cell_volumes[cell];
        if (!(volume > 0.0) || !std::isfinite(volume))
            throw std::runtime_error(
                "compute_scalar_gradient_with_boundary: invalid cell volume");
        gradient.component_data(0)[cell] = sum.x / volume;
        gradient.component_data(1)[cell] = sum.y / volume;
        gradient.component_data(2)[cell] = sum.z / volume;
    }
    return gradient;
}

// Integrated tangential diffusion flux per global face, shared by matrix
// assembly and post-solve flux reconstruction. Positive flux is outward from
// the owner cell. The implicit orthogonal coefficient and this correction use
// the same over-relaxed decomposition and variable-diffusion face conductance.
inline std::vector<double> compute_nonorthogonal_diffusion_face_fluxes(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& field,
    double diffusion_coefficient,
    const ScalarBoundaryConditions& boundary_conditions,
    const std::vector<double>* cell_diffusion = nullptr,
    const ScalarBoundaryFaceValues* face_values = nullptr,
    const ScalarBoundaryFaceConditions* face_conditions = nullptr)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells(), nf = mesh.n_faces();
    if (field.dimension() != 1 || field.size() != nc ||
        geometry.cell_centres.size() != nc || geometry.cell_volumes.size() != nc ||
        geometry.face_centres.size() != nf || geometry.face_area_vectors.size() != nf ||
        geometry.face_patch.size() != nf ||
        (cell_diffusion && cell_diffusion->size() != nc))
        throw std::invalid_argument(
            "compute_nonorthogonal_diffusion_face_fluxes: shape mismatch");
    if (!(diffusion_coefficient >= 0.0) || !std::isfinite(diffusion_coefficient))
        throw std::invalid_argument(
            "compute_nonorthogonal_diffusion_face_fluxes: invalid diffusion coefficient");
    bool active = diffusion_coefficient > 0.0;
    if (cell_diffusion) {
        active = false;
        for (double gamma : *cell_diffusion) {
            if (!(gamma >= 0.0) || !std::isfinite(gamma))
                throw std::invalid_argument(
                    "compute_nonorthogonal_diffusion_face_fluxes: invalid cell diffusion");
            active = active || gamma > 0.0;
        }
    }
    std::vector<double> flux(nf, 0.0);
    if (!active) return flux;

    // Green-Gauss with two-point interpolation is biased when face centroids
    // do not lie on owner-neighbour centre lines. Use least squares where the
    // cell has at least three internal neighbours (the LS kernel is rank-aware);
    // retain boundary-aware Green-Gauss on smaller stencils where prescribed
    // boundary data supply information absent from the cell-centred stencil.
    auto boundary_gradient = compute_scalar_gradient_with_boundary(
        field, mesh, geometry, boundary_conditions, face_values, face_conditions);
    auto gradient = boundary_gradient;
    const auto least_squares_gradient =
        cfdx::core::compute_gradient_least_squares(field, mesh);
    const auto& stencil_ownership = mesh.ownership();
    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    for (std::size_t cell = 0; cell < nc; ++cell) {
        std::size_t internal_neighbours = 0;
        for (Offset k = cell_offsets[cell]; k < cell_offsets[cell + 1]; ++k) {
            const std::size_t face = cell_faces[k];
            const auto neighbour = stencil_ownership.neighbour(face);
            if (neighbour >= 0 || stencil_ownership.owner(face) != cell)
                ++internal_neighbours;
        }
        if (internal_neighbours >= 3) {
            for (std::size_t direction = 0; direction < 3; ++direction)
                gradient.component_data(direction)[cell] =
                    least_squares_gradient.component_data(direction)[cell];
        }
    }
    const auto& ownership = mesh.ownership();
    for (std::size_t face = 0; face < nf; ++face) {
        const std::size_t owner = ownership.owner(face);
        if (owner >= nc)
            throw std::runtime_error(
                "compute_nonorthogonal_diffusion_face_fluxes: invalid owner");
        const auto neighbour_raw = ownership.neighbour(face);
        if (neighbour_raw < 0) {
            const std::size_t patch = geometry.face_patch[face];
            if (patch < mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch).type == PatchType::EMPTY)
                continue;
            const auto bc = resolve_scalar_boundary_condition(
                mesh, geometry, face, boundary_conditions, face_values, face_conditions);
            if (bc.type == ScalarBoundaryType::ZERO_GRADIENT)
                continue;
        }

        const Vec3 sf = geometry.face_area_vectors[face];
        const double area = sf.mag();
        const Vec3 d = neighbour_raw >= 0
            ? geometry.cell_centres[static_cast<std::size_t>(neighbour_raw)] -
                  geometry.cell_centres[owner]
            : geometry.face_centres[face] - geometry.cell_centres[owner];
        const double sd = sf.dot(d);
        const double d2 = d.dot(d);
        if (!(area > 0.0) || !(sd > 0.0) || !(d2 > 0.0) ||
            !std::isfinite(area) || !std::isfinite(sd) || !std::isfinite(d2))
            throw std::runtime_error(
                "compute_nonorthogonal_diffusion_face_fluxes: invalid face geometry");
        const double alpha = area * area / sd;
        const Vec3 tangential = sf - d * alpha;
        if (tangential.mag() <= 1e-12 * area) continue;

        double gamma_face = cell_diffusion
            ? (*cell_diffusion)[owner] : diffusion_coefficient;
        Vec3 grad_face{
            gradient.component_data(0)[owner],
            gradient.component_data(1)[owner],
            gradient.component_data(2)[owner]};
        if (neighbour_raw >= 0) {
            const std::size_t neighbour = static_cast<std::size_t>(neighbour_raw);
            if (neighbour >= nc)
                throw std::runtime_error(
                    "compute_nonorthogonal_diffusion_face_fluxes: invalid neighbour");
            const double area_projection = std::abs(sf.dot(d)) / area;
            const double d_owner =
                std::abs((geometry.face_centres[face] -
                          geometry.cell_centres[owner]).dot(sf)) / area;
            const double d_neighbour =
                std::abs((geometry.cell_centres[neighbour] -
                          geometry.face_centres[face]).dot(sf)) / area;
            if (!(d_owner > 0.0) || !(d_neighbour > 0.0) ||
                !std::isfinite(d_owner) || !std::isfinite(d_neighbour) ||
                !(area_projection > 0.0))
                throw std::runtime_error(
                    "compute_nonorthogonal_diffusion_face_fluxes: invalid normal distances");
            double conductance = 0.0;
            if (cell_diffusion) {
                const double go = (*cell_diffusion)[owner];
                const double gn = (*cell_diffusion)[neighbour];
                if (go > 0.0 && gn > 0.0)
                    conductance = area / (d_owner / go + d_neighbour / gn);
            } else {
                conductance = diffusion_coefficient * area / (d_owner + d_neighbour);
            }
            gamma_face = conductance / alpha;
            const double dc = (geometry.face_centres[face] -
                               geometry.cell_centres[owner]).mag();
            const double dn = (geometry.cell_centres[neighbour] -
                               geometry.face_centres[face]).mag();
            if (!(dc + dn > 0.0) || !std::isfinite(dc + dn))
                throw std::runtime_error(
                    "compute_nonorthogonal_diffusion_face_fluxes: invalid gradient interpolation distance");
            const double w = dn / (dc + dn);
            const Vec3 grad_neighbour{
                gradient.component_data(0)[neighbour],
                gradient.component_data(1)[neighbour],
                gradient.component_data(2)[neighbour]};
            grad_face = grad_face * w + grad_neighbour * (1.0 - w);
        }

        if (neighbour_raw < 0) {
            // Boundary correction must honour resolved face/patch data even
            // when the cell-centred internal-face gradient uses least squares.
            grad_face = Vec3{
                boundary_gradient.component_data(0)[owner],
                boundary_gradient.component_data(1)[owner],
                boundary_gradient.component_data(2)[owner]};
        }
        const double correction = gamma_face * grad_face.dot(tangential);
        if (!std::isfinite(correction))
            throw std::runtime_error(
                "compute_nonorthogonal_diffusion_face_fluxes: non-finite correction");
        flux[face] = correction;
    }
    return flux;
}

// Assemble:
//     div(phi*psi) - div(gamma grad(psi)) = Su + Sp*psi
// with first-order upwind convection and a normal implicit diffusion
// conductance plus conservative non-orthogonal tangential correction.
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
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* convected_field = nullptr,
    const ScalarBoundaryFaceConditions* face_conditions = nullptr,
    bool apply_nonorthogonal_correction = true)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    const std::size_t nf = mesh.n_faces();

    if ((face_conditions && (face_conditions->conditions.size() != nf || face_conditions->specified.size() != nf)) ||
        face_flux.size() != nf || face_flux.dimension() != 1 ||
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
    cfdx::core::Field<double, cfdx::core::Location::FACE> reconstructed_face;
    const bool needs_field_reconstruction =
        convection_scheme != ConvectionScheme::UPWIND;
    if (needs_field_reconstruction) {
        if (convected_field == nullptr || convected_field->size() != nc || convected_field->dimension() != 1)
            throw std::invalid_argument(
                "assemble_scalar_equation: selected convection scheme requires a scalar convected field");
    }
    if (convection_scheme == ConvectionScheme::SECOND_ORDER_UPWIND ||
        convection_scheme == ConvectionScheme::TVD) {
        reconstructed_gradient = cfdx::core::compute_gradient_gauss(*convected_field, mesh);
    }
    if (convection_scheme == ConvectionScheme::CENTRAL ||
        convection_scheme == ConvectionScheme::BLENDED ||
        convection_scheme == ConvectionScheme::QUICK ||
        convection_scheme == ConvectionScheme::QUICK_BOUNDED) {
        const auto geometry_cache = cfdx::core::make_geometry_cache(mesh);
        const auto core_scheme =
            convection_scheme == ConvectionScheme::CENTRAL ? cfdx::core::InterpScheme::LINEAR :
            convection_scheme == ConvectionScheme::BLENDED ? cfdx::core::InterpScheme::BLENDED :
            convection_scheme == ConvectionScheme::QUICK ? cfdx::core::InterpScheme::QUICK :
            cfdx::core::InterpScheme::QUICK_BOUNDED;
        reconstructed_face = cfdx::core::interpolate_cell_to_face(
            *convected_field, mesh, geometry_cache, core_scheme, &face_flux,
            cfdx::core::LimiterType::NONE, nullptr, 0.75);
    }

    const auto& own = mesh.ownership();
    const double* phi = face_flux.component_data(0);

    bool diffusion_active = diffusion_coefficient > 0.0;
    if (cell_diffusion) {
        diffusion_active = false;
        for (const double gamma : *cell_diffusion) {
            if (!(gamma >= 0.0) || !std::isfinite(gamma))
                throw std::invalid_argument(
                    "assemble_scalar_equation: invalid cell diffusion");
            diffusion_active = diffusion_active || gamma > 0.0;
        }
    }
    bool needs_nonorthogonal_gradient = false;
    if (apply_nonorthogonal_correction && diffusion_active) {
        for (std::size_t face = 0; face < nf; ++face) {
            const auto nraw = own.neighbour(face);
            if (nraw < 0) {
                const std::size_t patch = geometry.face_patch[face];
                if (patch < mesh.boundary().n_patches() &&
                    mesh.boundary().patch(patch).type == PatchType::EMPTY)
                    continue;
                const auto bc = resolve_scalar_boundary_condition(
                    mesh, geometry, face, boundary_conditions,
                    face_values, face_conditions);
                if (bc.type == ScalarBoundaryType::ZERO_GRADIENT)
                    continue;
            }
            const std::size_t owner = own.owner(face);
            if (owner >= nc)
                throw std::runtime_error(
                    "assemble_scalar_equation: invalid diffusion face owner");
            const Vec3 sf = geometry.face_area_vectors[face];
            const double area = sf.mag();
            const Vec3 d = nraw >= 0
                ? geometry.cell_centres[static_cast<std::size_t>(nraw)] -
                      geometry.cell_centres[owner]
                : geometry.face_centres[face] - geometry.cell_centres[owner];
            const double sd = sf.dot(d);
            const double d2 = d.dot(d);
            if (!(area > 0.0) || !(sd > 0.0) || !(d2 > 0.0) ||
                !std::isfinite(area) || !std::isfinite(sd) || !std::isfinite(d2))
                throw std::runtime_error(
                    "assemble_scalar_equation: invalid diffusion decomposition");
            const double alpha = area * area / sd;
            const Vec3 tangential = sf - d * alpha;
            if (tangential.mag() > 1e-12 * area)
                needs_nonorthogonal_gradient = true;
        }
    }
    if (apply_nonorthogonal_correction &&
        needs_nonorthogonal_gradient && convected_field == nullptr)
        throw std::invalid_argument(
            "assemble_scalar_equation: non-orthogonal diffusion requires the current scalar field");
    if (apply_nonorthogonal_correction && needs_nonorthogonal_gradient &&
        (convected_field->size() != nc || convected_field->dimension() != 1))
        throw std::invalid_argument(
            "assemble_scalar_equation: invalid scalar field for non-orthogonal diffusion");
    std::vector<double> nonorthogonal_diffusion_flux(nf, 0.0);
    if (apply_nonorthogonal_correction && needs_nonorthogonal_gradient) {
        nonorthogonal_diffusion_flux = compute_nonorthogonal_diffusion_face_fluxes(
            mesh, geometry, *convected_field, diffusion_coefficient,
            boundary_conditions, cell_diffusion, face_values, face_conditions);
    }

    double max_abs_face_flux = 0.0;
    double min_internal_distance = std::numeric_limits<double>::infinity();
    double max_internal_distance = 0.0;
    double min_face_area = std::numeric_limits<double>::infinity();
    double max_face_area = 0.0;

    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t o = own.owner(f);
        const auto nraw = own.neighbour(f);
        const double F = phi[f];
        max_abs_face_flux = std::max(max_abs_face_flux, std::abs(F));
        if (!std::isfinite(F)) throw std::invalid_argument("assemble_scalar_equation: non-finite face flux");

        if (nraw >= 0) {
            const std::size_t n = static_cast<std::size_t>(nraw);
            const double distance = (geometry.cell_centres[n] - geometry.cell_centres[o]).mag();
            const double area = geometry.face_area_vectors[f].mag();
            min_internal_distance = std::min(min_internal_distance, distance);
            max_internal_distance = std::max(max_internal_distance, distance);
            min_face_area = std::min(min_face_area, area);
            max_face_area = std::max(max_face_area, area);
            if (!(distance > 0.0) || !(area > 0.0))
                throw std::runtime_error("assemble_scalar_equation: degenerate internal face");
            const auto& Sf = geometry.face_area_vectors[f];
            const auto& Cf = geometry.face_centres[f];
            const auto& Co = geometry.cell_centres[o];
            const auto& Cn = geometry.cell_centres[n];
            // The implicit orthogonal conductance is a series resistance:
            // A / (d_Pf/Gamma_P + d_Nf/Gamma_N). Euclidean centre distance
            // and an unweighted harmonic mean are only equivalent on a
            // uniform orthogonal mesh.
            const double d_owner = std::abs((Cf - Co).dot(Sf)) / area;
            const double d_neighbour = std::abs((Cn - Cf).dot(Sf)) / area;
            if (!(d_owner > 0.0) || !(d_neighbour > 0.0) ||
                !std::isfinite(d_owner) || !std::isfinite(d_neighbour))
                throw std::runtime_error("assemble_scalar_equation: degenerate normal half-distance");
            double D = 0.0;
            if (cell_diffusion) {
                const double go = (*cell_diffusion)[o];
                const double gn = (*cell_diffusion)[n];
                if (!(go >= 0.0) || !(gn >= 0.0) ||
                    !std::isfinite(go) || !std::isfinite(gn))
                    throw std::invalid_argument("assemble_scalar_equation: invalid cell diffusion");
                if (go > 0.0 && gn > 0.0)
                    D = area / (d_owner / go + d_neighbour / gn);
            } else {
                D = diffusion_coefficient * area / (d_owner + d_neighbour);
            }

            const double a_on = D + std::max(F, 0.0);
            const double a_no = D + std::max(-F, 0.0);
            diag[o] += a_on;
            diag[n] += a_no;
            rows[o][n] -= a_no;
            rows[n][o] -= a_on;

            div_phi[o] += F;
            div_phi[n] -= F;

            if (needs_nonorthogonal_gradient) {
                // The implicit matrix represents -div(F_orth); moving the
                // shared tangential face flux to the RHS contributes
                // +div(F_tangent), with conservative owner/neighbour signs.
                deferred_rhs[o] += nonorthogonal_diffusion_flux[f];
                deferred_rhs[n] -= nonorthogonal_diffusion_flux[f];
            }

            if (convection_scheme == ConvectionScheme::SECOND_ORDER_UPWIND ||
                convection_scheme == ConvectionScheme::TVD ||
                convection_scheme == ConvectionScheme::CENTRAL ||
                convection_scheme == ConvectionScheme::BLENDED ||
                convection_scheme == ConvectionScheme::QUICK ||
                convection_scheme == ConvectionScheme::QUICK_BOUNDED) {
                const std::size_t upwind = F >= 0.0 ? o : n;
                const std::size_t downwind = F >= 0.0 ? n : o;
                const double phi_up = (*convected_field)(upwind);
                const double phi_other = (*convected_field)(downwind);
                double phi_high = phi_up;

                if (convection_scheme == ConvectionScheme::SECOND_ORDER_UPWIND) {
                    const auto& C_up = geometry.cell_centres[upwind];
                    const auto& Cf = geometry.face_centres[f];
                    const double* gx = reconstructed_gradient.component_data(0);
                    const double* gy = reconstructed_gradient.component_data(1);
                    const double* gz = reconstructed_gradient.component_data(2);
                    const double high_increment =
                        gx[upwind] * (Cf.x - C_up.x) +
                        gy[upwind] * (Cf.y - C_up.y) +
                        gz[upwind] * (Cf.z - C_up.z);
                    phi_high = phi_up + high_increment;
                    // Keep the historical bounded SOU reconstruction.
                    phi_high = std::clamp(
                        phi_high, std::min(phi_up, phi_other), std::max(phi_up, phi_other));
                } else if (convection_scheme == ConvectionScheme::TVD) {
                    // Darwish-Moukalled successive-slope ratio for an
                    // unstructured cell-centred face. The gradient must be
                    // projected over the full upwind-to-downwind centre
                    // distance, not the half-distance from the upwind centre
                    // to the face. For an exact linear field this gives r=1,
                    // so MINMOD retains the centred linear face value.
                    const auto& C_up = geometry.cell_centres[upwind];
                    const auto& Cf = geometry.face_centres[f];
                    const double* gx = reconstructed_gradient.component_data(0);
                    const double* gy = reconstructed_gradient.component_data(1);
                    const double* gz = reconstructed_gradient.component_data(2);
                    const double high_increment =
                        gx[upwind] * (Cf.x - C_up.x) +
                        gy[upwind] * (Cf.y - C_up.y) +
                        gz[upwind] * (Cf.z - C_up.z);
                    const double delta_down = phi_other - phi_up;
                    const auto& C_down = geometry.cell_centres[downwind];
                    const double delta_upwind_to_down =
                        gx[upwind] * (C_down.x - C_up.x) +
                        gy[upwind] * (C_down.y - C_up.y) +
                        gz[upwind] * (C_down.z - C_up.z);
                    const double ratio = std::abs(delta_down) > 1e-14
                        ? (2.0 * delta_upwind_to_down - delta_down) / delta_down
                        : 0.0;
                    const double limiter =
                        cfdx::core::limiter_psi(ratio, cfdx::core::LimiterType::MINMOD);
                    phi_high = phi_up + limiter * high_increment;
                    phi_high = std::clamp(
                        phi_high, std::min(phi_up, phi_other), std::max(phi_up, phi_other));
                } else {
                    // Central, blended and quadratic schemes are explicit
                    // deferred corrections around the conservative upwind
                    // matrix. The face value itself is produced by the core
                    // reconstruction library, so assembly and post-solve
                    // flux reconstruction share exactly the same scheme.
                    phi_high = reconstructed_face(f);
                }

                const double correction = F * (phi_high - phi_up);
                deferred_rhs[o] -= correction;
                deferred_rhs[n] += correction;
            }
        } else {
            const std::size_t patch = geometry.face_patch[f];
            if (patch < mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch).type == PatchType::EMPTY)
                continue;

            const ScalarBoundaryCondition bc = resolve_scalar_boundary_condition(
                mesh, geometry, f, boundary_conditions, face_values, face_conditions);

            const double area = geometry.face_area_vectors[f].mag();
            const double distance = boundary_normal_distance(geometry, f, o);
            min_face_area = std::min(min_face_area, area);
            max_face_area = std::max(max_face_area, area);
            if (!(distance > 0.0) || !(area > 0.0))
                throw std::runtime_error("assemble_scalar_equation: degenerate boundary face");

            if (needs_nonorthogonal_gradient)
                deferred_rhs[o] += nonorthogonal_diffusion_flux[f];

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
                const double gamma_owner = cell_diffusion
                    ? (*cell_diffusion)[o] : diffusion_coefficient;
                if (!(gamma_owner >= 0.0) || !std::isfinite(gamma_owner))
                    throw std::invalid_argument(
                        "assemble_scalar_equation: invalid boundary diffusion");
                // phi_b = phi_P + grad(phi).d: convection and diffusion use
                // the same boundary state.
                diag[o] += F;
                rhs[o] -= F * bc.gradient * distance;
                rhs[o] += gamma_owner * area * bc.gradient;
                div_phi[o] += F;
            } else {
                // Zero-gradient: phi_b = phi_P. For bounded convection the
                // later -div(phi)*phi_P correction cancels this term. For the
                // unbounded form, defer inflow rather than creating a sink or
                // a non-positive diagonal.
                if (bounded_convection) {
                    diag[o] += F;
                } else {
                    diag[o] += std::max(F, 0.0);
                    if (F < 0.0 && convected_field)
                        deferred_rhs[o] -= F * (*convected_field)(o);
                }
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
    eq.max_abs_face_flux = max_abs_face_flux;
    eq.min_internal_distance = min_internal_distance;
    eq.max_internal_distance = max_internal_distance;
    eq.min_face_area = min_face_area;
    eq.max_face_area = max_face_area;

    double signed_global = 0.0;
    double global_l1 = 0.0;
    double local_linf = 0.0;
    for (const double residual : div_phi) {
        if (!std::isfinite(residual))
            throw std::invalid_argument("assemble_scalar_equation: non-finite cell flux residual");
        signed_global += residual;
        global_l1 += std::abs(residual);
        local_linf = std::max(local_linf, std::abs(residual));
    }
    eq.signed_global_imbalance = signed_global;
    eq.global_l1_imbalance = global_l1;
    eq.max_local_imbalance = local_linf;
    eq.max_imbalance = local_linf;
    return eq;
}


// Reconstruct the conservative scalar transport flux directly from the final
// field and authoritative face mass flux. This is intentionally separate from
// the matrix assembly path: it is suitable for post-solve conservation audits.
inline cfdx::core::Field<double, cfdx::core::Location::FACE>
reconstruct_scalar_transport_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& face_flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& field,
    double diffusion_coefficient,
    const ScalarBoundaryConditions& boundary_conditions = {},
    bool bounded_convection = true,
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND,
    const ScalarBoundaryFaceConditions* face_conditions = nullptr,
    const std::vector<double>* cell_diffusion = nullptr,
    const ScalarBoundaryFaceValues* face_values = nullptr)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells(), nf = mesh.n_faces();
    if (face_flux.size() != nf || face_flux.dimension() != 1 ||
        field.size() != nc || field.dimension() != 1 ||
        geometry.face_centres.size() != nf ||
        geometry.face_area_vectors.size() != nf ||
        geometry.cell_centres.size() != nc ||
        geometry.cell_volumes.size() != nc ||
        geometry.face_patch.size() != nf)
        throw std::invalid_argument("reconstruct_scalar_transport_flux: dimensions do not match mesh");
    if (!(diffusion_coefficient >= 0.0) || !std::isfinite(diffusion_coefficient))
        throw std::invalid_argument("reconstruct_scalar_transport_flux: invalid diffusion coefficient");
    if (cell_diffusion && cell_diffusion->size() != nc)
        throw std::invalid_argument("reconstruct_scalar_transport_flux: invalid cell diffusion");

    Field<double, Location::FACE> result(nf, field.name()+"_flux", "", 1);
    const auto nonorthogonal_diffusion_flux = compute_nonorthogonal_diffusion_face_fluxes(
        mesh, geometry, field, diffusion_coefficient, boundary_conditions,
        cell_diffusion, face_values, face_conditions);
    Field<double, Location::CELL> grad;
    Field<double, Location::FACE> reconstructed_face;
    // Only SECOND_ORDER_UPWIND and TVD consume the cell gradient below.
    // CENTRAL/BLENDED/QUICK variants use reconstructed_face directly.
    if (convection_scheme == ConvectionScheme::SECOND_ORDER_UPWIND ||
        convection_scheme == ConvectionScheme::TVD)
        grad = compute_gradient_gauss(field, mesh);
    if (convection_scheme == ConvectionScheme::CENTRAL ||
        convection_scheme == ConvectionScheme::BLENDED ||
        convection_scheme == ConvectionScheme::QUICK ||
        convection_scheme == ConvectionScheme::QUICK_BOUNDED) {
        const auto geometry_cache = cfdx::core::make_geometry_cache(mesh);
        const auto core_scheme =
            convection_scheme == ConvectionScheme::CENTRAL ? cfdx::core::InterpScheme::LINEAR :
            convection_scheme == ConvectionScheme::BLENDED ? cfdx::core::InterpScheme::BLENDED :
            convection_scheme == ConvectionScheme::QUICK ? cfdx::core::InterpScheme::QUICK :
            cfdx::core::InterpScheme::QUICK_BOUNDED;
        reconstructed_face = cfdx::core::interpolate_cell_to_face(
            field, mesh, geometry_cache, core_scheme, &face_flux,
            cfdx::core::LimiterType::NONE, nullptr, 0.75);
    }

    const auto& own = mesh.ownership();
    for (std::size_t f = 0; f < nf; ++f) {
        const double F = face_flux(f);
        if (!std::isfinite(F))
            throw std::runtime_error("reconstruct_scalar_transport_flux: non-finite mass flux");
        const std::size_t o = own.owner(f);
        if (o >= nc) throw std::runtime_error("reconstruct_scalar_transport_flux: invalid owner");

        double gamma = diffusion_coefficient;
        double psi_face = field(o);
        const auto nraw = own.neighbour(f);
        if (nraw >= 0) {
            const std::size_t n = static_cast<std::size_t>(nraw);
            if (n >= nc) throw std::runtime_error("reconstruct_scalar_transport_flux: invalid neighbour");
            const double area = geometry.face_area_vectors[f].mag();
            if (!(area > 0.0) || !std::isfinite(area))
                throw std::runtime_error("reconstruct_scalar_transport_flux: degenerate internal face");
            const auto& Sf = geometry.face_area_vectors[f];
            const auto& Cf = geometry.face_centres[f];
            const double d_owner =
                std::abs((Cf - geometry.cell_centres[o]).dot(Sf)) / area;
            const double d_neighbour =
                std::abs((geometry.cell_centres[n] - Cf).dot(Sf)) / area;
            if (!(d_owner > 0.0) || !(d_neighbour > 0.0) ||
                !std::isfinite(d_owner) || !std::isfinite(d_neighbour))
                throw std::runtime_error(
                    "reconstruct_scalar_transport_flux: degenerate normal half-distance");
            double diffusion_conductance = 0.0;
            if (cell_diffusion) {
                const double go = (*cell_diffusion)[o], gn = (*cell_diffusion)[n];
                if (!(go >= 0.0) || !(gn >= 0.0) ||
                    !std::isfinite(go) || !std::isfinite(gn))
                    throw std::invalid_argument(
                        "reconstruct_scalar_transport_flux: invalid cell diffusion");
                if (go > 0.0 && gn > 0.0)
                    diffusion_conductance = area / (d_owner / go + d_neighbour / gn);
            } else {
                diffusion_conductance =
                    diffusion_coefficient * area / (d_owner + d_neighbour);
            }
            const double up = F >= 0.0 ? field(o) : field(n);
            psi_face = up;
            if (convection_scheme == ConvectionScheme::SECOND_ORDER_UPWIND ||
                convection_scheme == ConvectionScheme::TVD ||
                convection_scheme == ConvectionScheme::CENTRAL ||
                convection_scheme == ConvectionScheme::BLENDED ||
                convection_scheme == ConvectionScheme::QUICK ||
                convection_scheme == ConvectionScheme::QUICK_BOUNDED) {
                const std::size_t u = F >= 0.0 ? o : n;
                const std::size_t dcell = F >= 0.0 ? n : o;
                if (convection_scheme == ConvectionScheme::SECOND_ORDER_UPWIND ||
                    convection_scheme == ConvectionScheme::TVD) {
                    const double gx=grad.component_data(0)[u];
                    const double gy=grad.component_data(1)[u];
                    const double gz=grad.component_data(2)[u];
                    const auto& Cu=geometry.cell_centres[u];
                    const auto& Cf=geometry.face_centres[f];
                    const double increment=gx*(Cf.x-Cu.x)+gy*(Cf.y-Cu.y)+gz*(Cf.z-Cu.z);
                    double high=field(u)+increment;
                    if (convection_scheme == ConvectionScheme::TVD) {
                        const auto& Cd=geometry.cell_centres[dcell];
                        const double delta=field(dcell)-field(u);
                        const double full=gx*(Cd.x-Cu.x)+gy*(Cd.y-Cu.y)+gz*(Cd.z-Cu.z);
                        const double ratio=std::abs(delta)>1e-14 ? (2.0*full-delta)/delta : 0.0;
                        high=field(u)+limiter_psi(ratio, LimiterType::MINMOD)*increment;
                    }
                    // SECOND_ORDER_UPWIND and TVD face states are bounded by
                    // scheme contract, independently of bounded_convection.
                    // That option controls a separate bounded-divergence term.
                    psi_face = std::clamp(
                        high, std::min(field(u), field(dcell)),
                        std::max(field(u), field(dcell)));
                } else {
                    psi_face = reconstructed_face(f);
                }
            }
            result(f)=F*psi_face - diffusion_conductance*(field(n)-field(o))
                - nonorthogonal_diffusion_flux[f];
        } else {
            const std::size_t patch=geometry.face_patch[f];
            if (patch < mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch).type == PatchType::EMPTY) {
                result(f)=0.0;
                continue;
            }
            const ScalarBoundaryCondition bc = resolve_scalar_boundary_condition(
                mesh, geometry, f, boundary_conditions, face_values, face_conditions);
            const double area=geometry.face_area_vectors[f].mag();
            const double d=boundary_normal_distance(geometry,f,o);
            if (!(area > 0.0) || !(d > 0.0))
                throw std::runtime_error("reconstruct_scalar_transport_flux: degenerate boundary face");
            gamma = cell_diffusion ? (*cell_diffusion)[o] : diffusion_coefficient;
            if (!(gamma >= 0.0) || !std::isfinite(gamma))
                throw std::invalid_argument("reconstruct_scalar_transport_flux: invalid boundary diffusion");
            if (bc.type == ScalarBoundaryType::FIXED_VALUE) {
                // For a fixed-value boundary, the prescribed value enters the
                // convective flux only on inflow. On outflow the upwind state
                // is the owner-cell value, exactly as in matrix assembly.
                // Diffusion still uses the prescribed boundary value in both
                // flow directions.
                const double convective_face_value = F < 0.0 ? bc.value : field(o);
                result(f) = F * convective_face_value
                    - gamma * area / d * (bc.value - field(o))
                    - nonorthogonal_diffusion_flux[f];
            } else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT) {
                psi_face=field(o)+bc.gradient*d;
                result(f)=F*psi_face - gamma*area*bc.gradient
                    - nonorthogonal_diffusion_flux[f];
            } else {
                // Zero-gradient: no diffusive boundary flux and psi_f=psi_P.
                psi_face=field(o);
                result(f)=F*psi_face;
            }
        }
    }
    return result;
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
        const double candidate_residual =
            std::abs(diagonal * candidate_value - equation.rhs(0));
        if (candidate_residual <= controls.tolerance) {
            // A scalar 1x1 system has an exact closed-form solution. Do not
            // under-relax an already converged linear solve: doing so makes
            // the returned state fail the requested residual even though the
            // reported linear residual is zero.
            solution(0) = candidate_value;
            return {
                cfdx::core::SolverStatus::CONVERGED,
                1, candidate_residual, candidate_residual
            };
        }
        solution(0) += controls.relaxation * (candidate_value - solution(0));
        const double residual =
            std::abs(diagonal * solution(0) - equation.rhs(0));
        return {
            cfdx::core::SolverStatus::MAX_ITER_REACHED,
            controls.max_iterations, residual, residual
        };
    }

    cfdx::core::Vector candidate = solution;

    const auto selected = cfdx::core::solve_linear_system(
        equation.matrix, equation.rhs, candidate, controls.problem,
        controls.linear_solver, controls.max_iterations, controls.tolerance);
    auto result = selected.result;

    const bool fallback_allowed = controls.linear_solver.allow_fallback ||
        (controls.linear_solver.krylov == cfdx::core::KrylovModel::Auto &&
         controls.linear_solver.preconditioner ==
             cfdx::core::PreconditionerModel::Auto);

    // Keep all retries anchored to the same nonlinear iterate; the accepted
    // predictor is updated only after a solver reports convergence.
    // Momentum matrices can move between nearly symmetric diffusion-dominated
    // states and mildly nonsymmetric convection-dominated states. BiCGStab may
    // either stagnate or break down on the former even though the linear
    // system is well posed. Retry from the original nonlinear iterate with
    // restarted GMRES, then with CG as a final robust path. Every retry starts
    // from the original iterate so no unconverged Krylov state is injected.
    if (fallback_allowed &&
        result.status != cfdx::core::SolverStatus::CONVERGED) {
        candidate = solution;
        result = cfdx::core::solve_gmres(
            equation.matrix, equation.rhs, candidate,
            64, controls.max_iterations, controls.tolerance);
    }

    if (fallback_allowed &&
        result.status != cfdx::core::SolverStatus::CONVERGED) {
        candidate = solution;
        result = cfdx::core::solve_gauss_seidel(
            equation.matrix, equation.rhs, candidate,
            controls.max_iterations, controls.tolerance);
    }

    if (fallback_allowed &&
        result.status != cfdx::core::SolverStatus::CONVERGED) {
        candidate = solution;
        result = cfdx::core::solve_cg(
            equation.matrix, equation.rhs, candidate,
            controls.max_iterations, controls.tolerance);
    }

    // Very small coupled momentum systems (notably the analytical
    // acceptance meshes) can legitimately defeat every Krylov/stationary
    // applicability criterion while remaining perfectly nonsingular. Use a
    // bounded dense LU fallback only for small systems; production-size
    // systems continue to use the sparse iterative path above.
    if (fallback_allowed &&
        result.status != cfdx::core::SolverStatus::CONVERGED &&
        solution.size() <= 256) {
        const std::size_t n = solution.size();
        std::vector<double> a(n * n, 0.0);
        std::vector<double> b(n, 0.0);
        const auto* row = equation.matrix.row_offsets_data();
        const auto* col = equation.matrix.columns_data();
        const auto* val = equation.matrix.values_data();
        for (std::size_t i = 0; i < n; ++i) {
            b[i] = equation.rhs(i);
            for (std::uint32_t k = row[i]; k < row[i + 1]; ++k)
                a[i * n + col[k]] += val[k];
        }

        bool singular = false;
        for (std::size_t k = 0; k < n && !singular; ++k) {
            std::size_t pivot = k;
            double pivot_abs = std::abs(a[k * n + k]);
            for (std::size_t i = k + 1; i < n; ++i) {
                const double candidate_abs = std::abs(a[i * n + k]);
                if (candidate_abs > pivot_abs) {
                    pivot = i;
                    pivot_abs = candidate_abs;
                }
            }
            const double scale = std::max(1.0, pivot_abs);
            if (!(pivot_abs > 100.0 * std::numeric_limits<double>::epsilon() * scale) ||
                !std::isfinite(pivot_abs)) {
                singular = true;
                break;
            }
            if (pivot != k) {
                for (std::size_t j = k; j < n; ++j)
                    std::swap(a[k * n + j], a[pivot * n + j]);
                std::swap(b[k], b[pivot]);
            }
            for (std::size_t i = k + 1; i < n; ++i) {
                const double factor = a[i * n + k] / a[k * n + k];
                if (!std::isfinite(factor)) {
                    singular = true;
                    break;
                }
                a[i * n + k] = 0.0;
                for (std::size_t j = k + 1; j < n; ++j)
                    a[i * n + j] -= factor * a[k * n + j];
                b[i] -= factor * b[k];
            }
        }

        if (!singular) {
            candidate = solution;
            for (std::size_t ii = n; ii-- > 0;) {
                double value = b[ii];
                for (std::size_t j = ii + 1; j < n; ++j)
                    value -= a[ii * n + j] * candidate(j);
                const double d = a[ii * n + ii];
                if (!(std::abs(d) > 0.0) || !std::isfinite(d)) {
                    singular = true;
                    break;
                }
                candidate(ii) = value / d;
                if (!std::isfinite(candidate(ii))) {
                    singular = true;
                    break;
                }
            }
        }

        if (!singular) {
            double residual = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                double ri = -equation.rhs(i);
                for (std::uint32_t k = row[i]; k < row[i + 1]; ++k)
                    ri += val[k] * candidate(col[k]);
                residual = std::max(residual, std::abs(ri));
            }
            if (std::isfinite(residual)) {
                double rhs_norm = 0.0;
                for (std::size_t i = 0; i < n; ++i)
                    rhs_norm = std::max(rhs_norm, std::abs(equation.rhs(i)));
                const double relative = residual / std::max(rhs_norm, 1e-300);
                result = {
                    relative <= controls.tolerance
                        ? cfdx::core::SolverStatus::CONVERGED
                        : cfdx::core::SolverStatus::MAX_ITER_REACHED,
                    0, residual, relative
                };
            }
        }
    }

    if (result.status == cfdx::core::SolverStatus::CONVERGED) {
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
    // A MAX_ITER_REACHED result is not an acceptable predictor solution.
    // Do not inject an unconverged Krylov iterate into the nonlinear solver:
    // doing so can create an apparently finite but numerically meaningless
    // state that subsequently corrupts the face mass flux.
    return result;
}

inline double scalar_equation_residual_inf(
    const ScalarEquation& equation,
    const cfdx::core::Vector& solution)
{
    if (solution.size() != equation.rhs.size())
        throw std::invalid_argument("scalar_equation_residual_inf: dimension mismatch");
    if (equation.matrix.row_offsets_data() == nullptr ||
        equation.matrix.columns_data() == nullptr ||
        equation.matrix.values_data() == nullptr) {
        return std::numeric_limits<double>::infinity();
    }
    double max_r = 0.0;
    for (std::size_t i = 0; i < solution.size(); ++i) {
        const double rhs_i = equation.rhs(i);
        if (!std::isfinite(rhs_i)) return std::numeric_limits<double>::infinity();
        double r = -rhs_i;
        if (!std::isfinite(r)) return std::numeric_limits<double>::infinity();
        const std::uint32_t row_begin = equation.matrix.row_offsets_data()[i];
        const std::uint32_t row_end = equation.matrix.row_offsets_data()[i + 1];
        for (std::uint32_t k = row_begin; k < row_end; ++k) {
            const std::uint32_t cidx = equation.matrix.columns_data()[k];
            if (cidx >= solution.size()) return std::numeric_limits<double>::infinity();
            const double mval = equation.matrix.values_data()[k];
            const double sval = solution(cidx);
            if (!std::isfinite(mval) || !std::isfinite(sval))
                return std::numeric_limits<double>::infinity();
            const double term = mval * sval;
            if (!std::isfinite(term)) return std::numeric_limits<double>::infinity();
            r += term;
            if (!std::isfinite(r)) return std::numeric_limits<double>::infinity();
        }
        const double ar = std::abs(r);
        if (!std::isfinite(ar)) return std::numeric_limits<double>::infinity();
        max_r = std::max(max_r, ar);
    }
    return max_r;
}

} // namespace cfdx::physics
