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
#include "cfdx/core/fvm/least_squares_gradient.h"
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
    QUICK_BOUNDED,
    // Append the new variant to preserve all existing enum numeric values.
    // Smooth van Leer limiter for A/B convergence studies; not auto-selected.
    TVD_VANLEER,
    // Smooth TVD limiter, explicit opt-in for stationary convergence studies.
    TVD_VANALBADA
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

// Production gradient contract for the deferred non-orthogonal diffusion flux.
// Unlike the historic internal-only LS substitution, this includes the same
// face-specific Dirichlet/Neumann constraints as the assembled equation.
// Only rank- and conditioning-certified SVD fits are used; an unsuccessful fit
// falls back to the boundary-aware conservative Green-Gauss gradient and is
// counted explicitly. A fallback is NOT evidence of spatial qualification.
struct BoundaryWlsDiagnostics {
    std::size_t robust_cells = 0;
    std::size_t fallback_cells = 0;
    std::size_t rank_deficient_cells = 0;
    std::size_t ill_conditioned_cells = 0;
    double maximum_condition = 0.0;
};

inline cfdx::core::Field<double, cfdx::core::Location::CELL>
compute_boundary_aware_wls_gradient(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& field,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarBoundaryConditions& boundary_conditions,
    const ScalarBoundaryFaceValues* face_values = nullptr,
    const ScalarBoundaryFaceConditions* face_conditions = nullptr,
    BoundaryWlsDiagnostics* diagnostics = nullptr,
    double condition_limit = 1.e8,
    bool require_certified = false)
{
    using namespace cfdx::core;
    if (!(condition_limit >= 1.0) || !std::isfinite(condition_limit))
        throw std::invalid_argument("boundary WLS: invalid condition limit");
    auto gradient = compute_scalar_gradient_with_boundary(
        field, mesh, geometry, boundary_conditions, face_values, face_conditions);
    BoundaryWlsDiagnostics result;
    // An EMPTY patch denotes a 2D extrusion with the collapsed direction
    // intentionally omitted from the scalar transport equations.
    int dimension = 3;
    for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p)
        if (mesh.boundary().patch(p).type == PatchType::EMPTY)
            dimension = 2;
    const auto& own = mesh.ownership();
    const auto* offsets = mesh.cells().offsets_data();
    const auto* faces = mesh.cells().faces_data();
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        std::vector<Vec3> sample_points;
        std::vector<double> sample_values;
        const Vec3 centre = geometry.cell_centres[cell];
        const double value = field(cell);
        for (Offset k = offsets[cell]; k < offsets[cell + 1]; ++k) {
            const std::size_t face = faces[k];
            const auto nr = own.neighbour(face);
            const bool is_owner = own.owner(face) == cell;
            if (nr >= 0) {
                const std::size_t other = is_owner
                    ? static_cast<std::size_t>(nr)
                    : static_cast<std::size_t>(own.owner(face));
                if (other >= mesh.n_cells())
                    throw std::runtime_error("boundary WLS: invalid neighbour");
                sample_points.push_back(geometry.cell_centres[other]);
                sample_values.push_back(field(other));
                continue;
            }
            const std::size_t patch = geometry.face_patch[face];
            if (patch < mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch).type == PatchType::EMPTY)
                continue;
            const auto bc = resolve_scalar_boundary_condition(
                mesh, geometry, face, boundary_conditions,
                face_values, face_conditions);
            const Vec3 sf = is_owner ? geometry.face_area_vectors[face]
                                     : geometry.face_area_vectors[face] * (-1.0);
            const double area = sf.mag();
            if (!(area > 0.0) || !std::isfinite(area))
                throw std::runtime_error("boundary WLS: degenerate face");
            const Vec3 normal = sf * (1.0 / area);
            const double normal_distance =
                (geometry.face_centres[face] - centre).dot(normal);
            if (!(normal_distance > 0.0) || !std::isfinite(normal_distance))
                throw std::runtime_error("boundary WLS: degenerate normal distance");
            if (bc.type == ScalarBoundaryType::FIXED_VALUE) {
                sample_points.push_back(geometry.face_centres[face]);
                sample_values.push_back(bc.value);
            } else {
                // A Neumann datum constrains the NORMAL derivative, not a
                // tangential extrapolation to an off-axis face centroid.
                sample_points.push_back(centre + normal * normal_distance);
                sample_values.push_back(value + normal_distance *
                    (bc.type == ScalarBoundaryType::FIXED_GRADIENT ? bc.gradient : 0.0));
            }
        }
        StencilQuality quality;
        Vec3 value_gradient;
        bool certified = false;
        try {
            value_gradient = weighted_least_squares_gradient_robust(
                centre, value, sample_points, sample_values,
                GradientWeighting::INVERSE_DISTANCE_SQUARED,
                dimension, &quality);
            certified = quality.full_rank && quality.well_conditioned(condition_limit);
        } catch (const std::runtime_error&) {
            // Rejected rank or non-finite input is never silently accepted.
        }
        if (std::isfinite(quality.condition_estimate))
            result.maximum_condition = std::max(
                result.maximum_condition, quality.condition_estimate);
        if (!quality.full_rank) ++result.rank_deficient_cells;
        else if (!certified) ++result.ill_conditioned_cells;
        if (!certified) {
            ++result.fallback_cells;
            if (require_certified)
                throw std::runtime_error(
                    "boundary WLS: rank/conditioning gate rejected cell " +
                    std::to_string(cell));
            continue;
        }
        ++result.robust_cells;
        gradient.component_data(0)[cell] = value_gradient.x;
        gradient.component_data(1)[cell] = value_gradient.y;
        gradient.component_data(2)[cell] = value_gradient.z;
    }
    if (diagnostics) *diagnostics = result;
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
    const ScalarBoundaryFaceConditions* face_conditions = nullptr,
    const std::vector<double>* boundary_face_diffusion = nullptr)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells(), nf = mesh.n_faces();
    if (field.dimension() != 1 || field.size() != nc ||
        geometry.cell_centres.size() != nc || geometry.cell_volumes.size() != nc ||
        geometry.face_centres.size() != nf || geometry.face_area_vectors.size() != nf ||
        geometry.face_patch.size() != nf ||
        (cell_diffusion && cell_diffusion->size() != nc) ||
        (boundary_face_diffusion && boundary_face_diffusion->size() != nf))
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

    const auto gradient = compute_boundary_aware_wls_gradient(
        field, mesh, geometry, boundary_conditions, face_values, face_conditions);
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

        double gamma_face = neighbour_raw < 0 && boundary_face_diffusion
            ? (*boundary_face_diffusion)[face]
            : (cell_diffusion ? (*cell_diffusion)[owner] : diffusion_coefficient);
        if (!(gamma_face >= 0.0) || !std::isfinite(gamma_face))
            throw std::invalid_argument(
                "compute_nonorthogonal_diffusion_face_fluxes: invalid face diffusion");
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
            const Vec3 centre_line =
                geometry.cell_centres[neighbour] - geometry.cell_centres[owner];
            const double d2_line = centre_line.mag2();
            if (!(d2_line > 0.0) || !std::isfinite(d2_line))
                throw std::runtime_error(
                    "compute_nonorthogonal_diffusion_face_fluxes: invalid gradient interpolation distance");
            // Project the true face centroid onto the centre line. Euclidean
            // face-to-cell distances bias strongly stretched skewed faces.
            const double neighbour_weight = std::clamp(
                (geometry.face_centres[face] - geometry.cell_centres[owner]).dot(centre_line)
                    / d2_line, 0.0, 1.0);
            const double w = 1.0 - neighbour_weight;
            const Vec3 grad_neighbour{
                gradient.component_data(0)[neighbour],
                gradient.component_data(1)[neighbour],
                gradient.component_data(2)[neighbour]};
            grad_face = grad_face * w + grad_neighbour * (1.0 - w);
        }

        // The WLS gradient above includes resolved Dirichlet/Neumann
        // boundary data and falls back to boundary-aware Green-Gauss if its
        // stencil is not certified. Do not replace it unconditionally by
        // Green-Gauss on boundary faces: that destroys affine exactness on
        // stretched, non-uniform and skewed wall grids.
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
    bool apply_nonorthogonal_correction = true,
    const std::vector<double>* boundary_face_diffusion = nullptr,
    double high_order_fraction = 1.0)
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
        (cell_diffusion && cell_diffusion->size() != nc) ||
        (boundary_face_diffusion && boundary_face_diffusion->size() != nf))
        throw std::invalid_argument("assemble_scalar_equation: field dimensions do not match mesh");
    if (boundary_face_diffusion)
        for (const double gamma : *boundary_face_diffusion)
            if (!(gamma >= 0.0) || !std::isfinite(gamma))
                throw std::invalid_argument("assemble_scalar_equation: invalid boundary face viscosity");
    if (!(diffusion_coefficient >= 0.0) || !std::isfinite(diffusion_coefficient))
        throw std::invalid_argument("assemble_scalar_equation: invalid diffusion coefficient");
    if (!std::isfinite(high_order_fraction) || high_order_fraction < 0.0 ||
        high_order_fraction > 1.0)
        throw std::invalid_argument("assemble_scalar_equation: high-order fraction must be in [0,1]");

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
        convection_scheme == ConvectionScheme::TVD ||
        convection_scheme == ConvectionScheme::TVD_VANLEER ||
                convection_scheme == ConvectionScheme::TVD_VANALBADA) {
        // Use the same boundary-aware rank-certified reconstruction as
        // deferred diffusion. Historical raw Gauss ignored physical BCs and
        // can lose affine accuracy on skewed production meshes.
        reconstructed_gradient = compute_boundary_aware_wls_gradient(
            *convected_field, mesh, geometry, boundary_conditions,
            face_values, face_conditions);
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
            boundary_conditions, cell_diffusion, face_values, face_conditions,
            boundary_face_diffusion);
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
                convection_scheme == ConvectionScheme::TVD_VANLEER ||
                convection_scheme == ConvectionScheme::TVD_VANALBADA ||
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
                } else if (convection_scheme == ConvectionScheme::TVD ||
                           convection_scheme == ConvectionScheme::TVD_VANLEER ||
                convection_scheme == ConvectionScheme::TVD_VANALBADA) {
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
                    const double limiter = std::max(0.0,
                        cfdx::core::limiter_psi(ratio,
                            convection_scheme == ConvectionScheme::TVD_VANALBADA
                            ? cfdx::core::LimiterType::VAN_ALBADA
                            : convection_scheme == ConvectionScheme::TVD_VANLEER
                                ? cfdx::core::LimiterType::VANLEER
                                : cfdx::core::LimiterType::MINMOD));
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

                // A conservative deferred-correction homotopy. The implicit
                // UPWIND matrix stays unchanged and only the signed face
                // correction is ramped; owner/neighbour RHS increments remain
                // exactly opposite at every fraction in [0,1].
                const double correction = high_order_fraction * F * (phi_high - phi_up);
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
                const double gamma_owner = boundary_face_diffusion
                    ? (*boundary_face_diffusion)[f]
                    : (cell_diffusion ? (*cell_diffusion)[o] : diffusion_coefficient);
                if (gamma_owner < 0.0)
                    throw std::invalid_argument("assemble_scalar_equation: negative boundary diffusion");
                const double D = gamma_owner * area / distance;
                const double a_owner = D + std::max(F, 0.0);
                const double a_boundary = D + std::max(-F, 0.0);
                diag[o] += a_owner;
                rhs[o] += a_boundary * bc.value;
                div_phi[o] += F;
            } else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT) {
                const double gamma_owner = boundary_face_diffusion
                    ? (*boundary_face_diffusion)[f]
                    : (cell_diffusion ? (*cell_diffusion)[o] : diffusion_coefficient);
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
    bool /*bounded_convection*/ = true,
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND,
    const ScalarBoundaryFaceConditions* face_conditions = nullptr,
    const std::vector<double>* cell_diffusion = nullptr,
    const ScalarBoundaryFaceValues* face_values = nullptr,
    const std::vector<double>* boundary_face_diffusion = nullptr,
    double high_order_fraction = 1.0)
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
    if (!std::isfinite(high_order_fraction) || high_order_fraction < 0.0 ||
        high_order_fraction > 1.0)
        throw std::invalid_argument("reconstruct_scalar_transport_flux: high-order fraction must be in [0,1]");
    if (cell_diffusion && cell_diffusion->size() != nc)
        throw std::invalid_argument("reconstruct_scalar_transport_flux: invalid cell diffusion");
    if (boundary_face_diffusion && boundary_face_diffusion->size() != nf)
        throw std::invalid_argument("reconstruct_scalar_transport_flux: invalid boundary face diffusion");

    Field<double, Location::FACE> result(nf, field.name()+"_flux", "", 1);
    const auto nonorthogonal_diffusion_flux = compute_nonorthogonal_diffusion_face_fluxes(
        mesh, geometry, field, diffusion_coefficient, boundary_conditions,
        cell_diffusion, face_values, face_conditions, boundary_face_diffusion);
    Field<double, Location::CELL> grad;
    Field<double, Location::FACE> reconstructed_face;
    // SOU and both TVD variants consume the cell gradient below.
    // CENTRAL/BLENDED/QUICK variants use reconstructed_face directly.
    if (convection_scheme == ConvectionScheme::SECOND_ORDER_UPWIND ||
        convection_scheme == ConvectionScheme::TVD ||
        convection_scheme == ConvectionScheme::TVD_VANLEER ||
                convection_scheme == ConvectionScheme::TVD_VANALBADA)
        // Match assemble_scalar_equation exactly: an independent flux audit
        // must use the SAME boundary-aware gradient and limiter inputs.
        grad = compute_boundary_aware_wls_gradient(
            field, mesh, geometry, boundary_conditions, face_values, face_conditions);
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
                convection_scheme == ConvectionScheme::TVD_VANLEER ||
                convection_scheme == ConvectionScheme::TVD_VANALBADA ||
                convection_scheme == ConvectionScheme::CENTRAL ||
                convection_scheme == ConvectionScheme::BLENDED ||
                convection_scheme == ConvectionScheme::QUICK ||
                convection_scheme == ConvectionScheme::QUICK_BOUNDED) {
                const std::size_t u = F >= 0.0 ? o : n;
                const std::size_t dcell = F >= 0.0 ? n : o;
                if (convection_scheme == ConvectionScheme::SECOND_ORDER_UPWIND ||
                    convection_scheme == ConvectionScheme::TVD ||
                    convection_scheme == ConvectionScheme::TVD_VANLEER ||
                convection_scheme == ConvectionScheme::TVD_VANALBADA) {
                    const double gx=grad.component_data(0)[u];
                    const double gy=grad.component_data(1)[u];
                    const double gz=grad.component_data(2)[u];
                    const auto& Cu=geometry.cell_centres[u];
                    const auto& Cf=geometry.face_centres[f];
                    const double increment=gx*(Cf.x-Cu.x)+gy*(Cf.y-Cu.y)+gz*(Cf.z-Cu.z);
                    double high=field(u)+increment;
                    if (convection_scheme == ConvectionScheme::TVD ||
                        convection_scheme == ConvectionScheme::TVD_VANLEER ||
                convection_scheme == ConvectionScheme::TVD_VANALBADA) {
                        const auto& Cd=geometry.cell_centres[dcell];
                        const double delta=field(dcell)-field(u);
                        const double full=gx*(Cd.x-Cu.x)+gy*(Cd.y-Cu.y)+gz*(Cd.z-Cu.z);
                        const double ratio=std::abs(delta)>1e-14 ? (2.0*full-delta)/delta : 0.0;
                        high=field(u)+std::max(0.0, limiter_psi(ratio,
                            convection_scheme == ConvectionScheme::TVD_VANALBADA
                            ? LimiterType::VAN_ALBADA
                            : convection_scheme == ConvectionScheme::TVD_VANLEER
                                ? LimiterType::VANLEER : LimiterType::MINMOD))*increment;
                    }
                    // SOU and both TVD face states are bounded by
                    // scheme contract, independently of bounded_convection.
                    // That option controls a separate bounded-divergence term.
                    psi_face = std::clamp(
                        high, std::min(field(u), field(dcell)),
                        std::max(field(u), field(dcell)));
                } else {
                    psi_face = reconstructed_face(f);
                }
            }
            // Reconstruct the SAME conservative fraction used in the matrix
            // assembly. The physical full-order audit requires fraction=1.
            psi_face = up + high_order_fraction * (psi_face - up);
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
            gamma = boundary_face_diffusion ? (*boundary_face_diffusion)[f]
                : (cell_diffusion ? (*cell_diffusion)[o] : diffusion_coefficient);
            if (!(gamma >= 0.0) || !std::isfinite(gamma))
