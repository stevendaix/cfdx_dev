#pragma once

// CFDX wall-distance service: explicit numerical provenance and fail-closed
// geometric qualification.  This layer does not alter any of the individual
// wall-distance discretizations in wall_distance.h.
#include "cfdx/physics/wall_distance.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

enum class WallDistanceCheckScope {
    ALL_FLUID_NODES,    // qualifies the returned discrete samples
    SAMPLED_MONITORING, // diagnostics only: never qualifies the entire field
    NONE                // no geometric accuracy claim
};

enum class WallDistanceFailureAction {
    REJECT,                  // return a result marked unusable
    EXPLICIT_FULL_BVH        // replace the WHOLE fluid field, never individual nodes
};

enum class WallDistanceDecision {
    EXACT_TRIANGULATED_SURFACE,
    VERIFIED_APPROXIMATION,
    MONITORING_ONLY,
    NOT_CHECKED,
    REJECTED_GEOMETRIC_ERROR,
    REJECTED_NUMERICAL_FAILURE,
    REJECTED_METRIC_MISMATCH,
    REJECTED_INVALID_FIELD,
    FULL_BVH_FALLBACK
};

struct WallDistanceMethodTraits {
    const char* family;
    bool direct_geometric;
    bool euclidean_distance_target;
    bool needs_active_pde_nodes;
};

inline WallDistanceMethodTraits wall_distance_method_traits(WallDistanceMethod method) {
    switch (method) {
    case WallDistanceMethod::EXACT_GEOMETRIC:
        return {"brute_force_point_triangle", true, true, false};
    case WallDistanceMethod::SEARCH_BASED:
        return {"bvh_point_triangle", true, true, false};
    case WallDistanceMethod::MESH_WAVE:
        return {"geometric_wall_point_propagation", false, true, false};
    case WallDistanceMethod::DIRECTIONAL_MESH_WAVE:
        // CFDX implements a penalized directional heuristic, not a standard
        // Euclidean distance and not OpenFOAM's specified-direction metric.
        return {"non_euclidean_directional_heuristic", false, false, false};
    case WallDistanceMethod::POISSON:
        return {"poisson_potential_reconstruction", false, true, true};
    case WallDistanceMethod::EIKONAL:
        return {"first_order_upwind_eikonal", false, true, true};
    case WallDistanceMethod::EIKONAL_SECOND_ORDER:
        return {"second_order_upwind_eikonal", false, true, true};
    case WallDistanceMethod::HAMILTON_JACOBI:
        return {"regularized_hamilton_jacobi", false, true, true};
    case WallDistanceMethod::ADVECTION_DIFFUSION:
        return {"regularized_advection_diffusion", false, true, true};
    case WallDistanceMethod::HYBRID_POISSON_EIKONAL:
        return {"poisson_hamilton_jacobi_hybrid", false, true, true};
    }
    throw std::invalid_argument("unknown wall-distance method");
}

struct WallDistanceServiceOptions {
    // Unlike the legacy generic compute_wall_distance overload, the service
    // defaults to the exact BVH geometry for physical wall-distance use.
    WallDistanceMethod method{WallDistanceMethod::SEARCH_BASED};
    std::size_t max_iterations{1000};
    WallDistanceCheckScope check_scope{WallDistanceCheckScope::ALL_FLUID_NODES};
    std::size_t sampling_stride{8};
    double absolute_tolerance{1e-12};
    double relative_tolerance{5e-3};
    WallDistanceFailureAction on_failure{WallDistanceFailureAction::REJECT};
};

struct WallDistanceQuality {
    WallDistanceDecision decision{WallDistanceDecision::NOT_CHECKED};
    std::string source_method;
    std::string selected_method;
    std::string method_family;
    std::size_t fluid_nodes{0};
    std::size_t geometry_checked_nodes{0};
    std::size_t geometry_tolerance_violations{0};
    std::size_t invalid_fluid_nodes{0};
    std::size_t active_pde_nodes{0};
    double relative_l2_error_checked{std::numeric_limits<double>::quiet_NaN()};
    double absolute_linf_error_checked{std::numeric_limits<double>::quiet_NaN()};
    bool solver_converged{false};
    bool solver_actually_exercised{false};
    bool full_geometry_check{false};
    bool geometry_within_tolerance{false};
    bool algorithm_qualified{false};
    bool approved_for_euclidean_wall_distance{false};
    bool entirely_replaced_by_bvh{false};
};

struct WallDistanceServiceResult {
    // Raw is always the actual requested algorithm, even after a fallback.
    WallDistanceResult raw;
    WallDistanceResult selected;
    WallDistanceQuality quality;
};

inline void validate_wall_distance_service_options(const WallDistanceServiceOptions& o) {
    (void)wall_distance_method_traits(o.method);
    if (o.max_iterations == 0 || o.sampling_stride == 0 ||
        !std::isfinite(o.absolute_tolerance) || o.absolute_tolerance < 0.0 ||
        !std::isfinite(o.relative_tolerance) || o.relative_tolerance < 0.0)
        throw std::invalid_argument("invalid wall-distance service controls");
    switch (o.check_scope) {
    case WallDistanceCheckScope::ALL_FLUID_NODES:
    case WallDistanceCheckScope::SAMPLED_MONITORING:
    case WallDistanceCheckScope::NONE: break;
    default: throw std::invalid_argument("unknown wall-distance check scope");
    }
    switch (o.on_failure) {
    case WallDistanceFailureAction::REJECT:
    case WallDistanceFailureAction::EXPLICIT_FULL_BVH: break;
    default: throw std::invalid_argument("unknown wall-distance fallback mode");
    }
}

inline void validate_wall_distance_service_grid(const WallDistanceGrid& grid) {
    const auto max=std::numeric_limits<std::size_t>::max();
    if (grid.nx<3 || grid.ny<3 || grid.nz<3 ||
        grid.ny>max/grid.nx || grid.nz>max/(grid.nx*grid.ny) ||
        !std::isfinite(grid.spacing.x) || !std::isfinite(grid.spacing.y) ||
        !std::isfinite(grid.spacing.z) || grid.spacing.x<=0.0 ||
        grid.spacing.y<=0.0 || grid.spacing.z<=0.0 ||
        !std::isfinite(grid.origin.x) || !std::isfinite(grid.origin.y) ||
        !std::isfinite(grid.origin.z) ||
        grid.points.size()!=grid.nx*grid.ny*grid.nz ||
        grid.solid.size()!=grid.points.size() ||
        grid.neighbours.size()!=grid.points.size())
        throw std::invalid_argument("invalid wall-distance service grid layout");
    for (std::size_t i=0; i<grid.points.size(); ++i) {
        const auto& p=grid.points[i];
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
            throw std::invalid_argument("wall-distance service grid has non-finite point");
        for (const auto next : grid.neighbours[i])
            if (next>=grid.points.size())
                throw std::invalid_argument("wall-distance service grid has invalid neighbour");
    }
}

// Also exported for regression tests, offline qualification and postprocessing.
// No part of this function silently modifies the requested algorithm's field.
inline WallDistanceServiceResult assess_wall_distance_field(
    WallDistanceResult raw, const WallDistanceMethod method,
    const WallSurface& surface, const WallDistanceGrid& grid,
    const WallDistanceServiceOptions& controls)
{
    validate_wall_distance_service_options(controls);
    if (method != controls.method)
        throw std::invalid_argument("wall-distance service method mismatch");
    validate_wall_distance_surface(surface);
    validate_wall_distance_service_grid(grid);
    const std::size_t n=grid.points.size();
    if (n==0 || grid.solid.size()!=n || raw.distance.size()!=n || raw.valid.size()!=n)
        throw std::invalid_argument("wall-distance service field/grid sizes mismatch");

    const auto traits=wall_distance_method_traits(method);
    WallDistanceServiceResult out;
    out.raw=std::move(raw);
    out.selected=out.raw;
    auto& q=out.quality;
    q.source_method=wall_distance_method_name(method);
    q.selected_method=out.selected.method;
    q.method_family=traits.family;
    q.solver_converged=out.raw.converged;
    q.active_pde_nodes=out.raw.active_pde_samples;

    for (std::size_t i=0;i<n;++i) {
        if (grid.solid[i]) continue;
        ++q.fluid_nodes;
        if (!out.raw.valid[i] || !std::isfinite(out.raw.distance[i]) ||
            out.raw.distance[i]<0.0) ++q.invalid_fluid_nodes;
    }
    if (q.fluid_nodes==0)
        throw std::invalid_argument("wall-distance service requires fluid nodes");

    q.solver_actually_exercised=!traits.needs_active_pde_nodes || q.active_pde_nodes!=0;
    q.algorithm_qualified=q.solver_converged && q.solver_actually_exercised &&
        q.invalid_fluid_nodes==0 && traits.euclidean_distance_target;

    const bool direct=traits.direct_geometric;
    q.full_geometry_check=direct || controls.check_scope==WallDistanceCheckScope::ALL_FLUID_NODES;
    q.geometry_within_tolerance=direct && q.invalid_fluid_nodes==0;

    if (!direct && controls.check_scope!=WallDistanceCheckScope::NONE) {
        const WallDistanceBvh oracle(surface);
        long double err2=0.0L, ref2=0.0L;
        long double max_error=0.0L;
        std::size_t ordinal=0;
        for (std::size_t i=0;i<n;++i) {
            if (grid.solid[i]) continue;
            const bool check=controls.check_scope==WallDistanceCheckScope::ALL_FLUID_NODES ||
                (ordinal % controls.sampling_stride)==0 || ordinal+1==q.fluid_nodes;
            ++ordinal;
            if (!check) continue;
            ++q.geometry_checked_nodes;
            const double exact=oracle.nearest_distance(grid.points[i]);
            if (!std::isfinite(exact) || exact<0.0)
                throw std::runtime_error("invalid wall-distance BVH oracle");
            const bool valid=out.raw.valid[i] && std::isfinite(out.raw.distance[i]) &&
                out.raw.distance[i]>=0.0;
            const long double diff=valid
                ? std::abs(static_cast<long double>(out.raw.distance[i])-exact)
                : std::numeric_limits<long double>::infinity();
            const long double tol=static_cast<long double>(controls.absolute_tolerance) +
                static_cast<long double>(controls.relative_tolerance)*exact;
            if (!(diff<=tol)) ++q.geometry_tolerance_violations;
            err2+=diff*diff;
            ref2+=static_cast<long double>(exact)*exact;
            max_error=std::max(max_error,diff);
        }
        q.relative_l2_error_checked=ref2>0.0L
            ? static_cast<double>(std::sqrt(err2/ref2))
            : (err2==0.0L ? 0.0 : std::numeric_limits<double>::infinity());
        q.absolute_linf_error_checked=static_cast<double>(max_error);
        q.geometry_within_tolerance=q.invalid_fluid_nodes==0 &&
            q.geometry_tolerance_violations==0 && q.geometry_checked_nodes!=0;
    } else if (direct) {
        q.geometry_checked_nodes=q.fluid_nodes;
        q.relative_l2_error_checked=0.0;
        q.absolute_linf_error_checked=0.0;
    }

    if (direct && q.geometry_within_tolerance && out.raw.converged) {
        q.decision=WallDistanceDecision::EXACT_TRIANGULATED_SURFACE;
        q.approved_for_euclidean_wall_distance=true;
    } else if (q.invalid_fluid_nodes!=0) {
        q.decision=WallDistanceDecision::REJECTED_INVALID_FIELD;
    } else if (!traits.euclidean_distance_target) {
        q.decision=WallDistanceDecision::REJECTED_METRIC_MISMATCH;
    } else if (!q.solver_converged || !q.solver_actually_exercised) {
        q.decision=WallDistanceDecision::REJECTED_NUMERICAL_FAILURE;
    } else if (controls.check_scope==WallDistanceCheckScope::NONE) {
        q.decision=WallDistanceDecision::NOT_CHECKED;
    } else if (!q.full_geometry_check) {
        q.decision=WallDistanceDecision::MONITORING_ONLY;
    } else if (!q.geometry_within_tolerance) {
        q.decision=WallDistanceDecision::REJECTED_GEOMETRIC_ERROR;
    } else {
        q.decision=WallDistanceDecision::VERIFIED_APPROXIMATION;
        q.approved_for_euclidean_wall_distance=true;
    }

    if (!q.approved_for_euclidean_wall_distance &&
        controls.on_failure==WallDistanceFailureAction::EXPLICIT_FULL_BVH) {
        // Explicit, WHOLE-field fallback. Raw remains untouched and retains
        // its original residual, flags, and method. Not a PDE validation.
        out.selected=compute_wall_distance(WallDistanceMethod::SEARCH_BASED,
                                           surface,grid,controls.max_iterations);
        if (!out.selected.converged || !out.selected.geometrically_exact_on_input_mesh)
            throw std::runtime_error("exact BVH fallback did not produce a valid field");
        for (std::size_t i=0;i<n;++i) if (!grid.solid[i] &&
            (!out.selected.valid[i] || !std::isfinite(out.selected.distance[i]) ||
             out.selected.distance[i]<0.0))
            throw std::runtime_error("exact BVH fallback returned invalid fluid distance");
        q.decision=WallDistanceDecision::FULL_BVH_FALLBACK;
        q.entirely_replaced_by_bvh=true;
        q.selected_method="search_based_bvh_full_fallback";
        out.selected.method=q.selected_method;
        q.approved_for_euclidean_wall_distance=true;
    }
    return out;
}

inline WallDistanceServiceResult audit_wall_distance_service(
    const WallSurface& surface,const WallDistanceGrid& grid,
    const WallDistanceServiceOptions& controls={})
{
    validate_wall_distance_service_options(controls);
    validate_wall_distance_surface(surface);
    validate_wall_distance_service_grid(grid);
    return assess_wall_distance_field(
        compute_wall_distance(controls.method,surface,grid,controls.max_iterations),
        controls.method,surface,grid,controls);
}

// Production-facing fail-closed API: an unqualified approximate field cannot
// leak to SST or y+ consumers even if its internal linear solver converged.
inline WallDistanceServiceResult require_qualified_wall_distance(
    const WallSurface& surface,const WallDistanceGrid& grid,
    const WallDistanceServiceOptions& controls={})
{
    auto r=audit_wall_distance_service(surface,grid,controls);
    if (!r.quality.approved_for_euclidean_wall_distance)
        throw std::domain_error("wall-distance field is not geometrically qualified: "+
            r.quality.source_method);
    return r;
}

} // namespace cfdx::physics
