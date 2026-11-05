#include "cfdx/physics/wall_distance_service.h"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace cfdx::physics;

static void require(bool ok,const char* message) {
    if (!ok) throw std::runtime_error(message);
}

static WallSurface cube() {
    WallSurface s;
    s.points={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
              {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    constexpr std::array<std::array<std::size_t,4>,6> faces={{{{0,1,2,3}},
        {{4,5,6,7}},{{0,1,5,4}},{{1,2,6,5}},{{2,3,7,6}},{{3,0,4,7}}}};
    for (const auto& f:faces) {
        s.triangles.push_back({{f[0],f[1],f[2]}});
        s.triangles.push_back({{f[0],f[2],f[3]}});
    }
    return s;
}

int main() {
    const auto surface=cube();
    const auto grid=make_wall_distance_grid(9,9,9,{-0.5,-0.5,-0.5},
        {0.25,0.25,0.25},[](WallDistanceVec3 p) {
            return p.x>0&&p.x<1&&p.y>0&&p.y<1&&p.z>0&&p.z<1;
        });

    // Production default: exact BVH, no ambiguity with legacy Eikonal default.
    const auto legacy_default=compute_wall_distance(surface,grid,20);
    require(legacy_default.method=="search_based" &&
            legacy_default.geometrically_exact_on_input_mesh,
            "generic dispatcher did not switch to the exact default");
    const auto safe=require_qualified_wall_distance(surface,grid);
    require(safe.quality.approved_for_euclidean_wall_distance,"default not approved");
    require(safe.quality.decision==WallDistanceDecision::EXACT_TRIANGULATED_SURFACE,
            "default not exact");
    require(safe.quality.method_family=="bvh_point_triangle","wrong default route");
    require(safe.selected.geometrically_exact_on_input_mesh,"missing exact provenance");
    require(!safe.quality.entirely_replaced_by_bvh,"direct BVH must not be called fallback");

    // A raw Eikonal field and a BVH field are DIFFERENT algorithms. Neither
    // a tiny PDE residual nor the geometric oracle may rewrite raw outputs.
    WallDistanceServiceOptions eik;
    eik.method=WallDistanceMethod::EIKONAL;
    eik.max_iterations=500;
    eik.absolute_tolerance=0.0;
    eik.relative_tolerance=1e-12;
    const auto unqualified=audit_wall_distance_service(surface,grid,eik);
    require(unqualified.quality.fluid_nodes>0,"empty fluid test");
    require(unqualified.quality.geometry_checked_nodes==unqualified.quality.fluid_nodes,
            "not a full-field check");
    require(unqualified.quality.geometry_tolerance_violations>0,"no Eikonal accuracy test");
    require(!unqualified.quality.approved_for_euclidean_wall_distance,"inaccurate PDE accepted");
    require(!unqualified.quality.entirely_replaced_by_bvh,"unexpected BVH substitution");
    require(unqualified.raw.distance==unqualified.selected.distance,"raw distances were rewritten");
    require(unqualified.quality.active_pde_nodes>0,"Eikonal active nodes not reported");
    bool threw=false;
    try { (void)require_qualified_wall_distance(surface,grid,eik); }
    catch (const std::domain_error&) { threw=true; }
    require(threw,"production route accepted inaccurate Eikonal");

    eik.on_failure=WallDistanceFailureAction::EXPLICIT_FULL_BVH;
    const auto fallback=require_qualified_wall_distance(surface,grid,eik);
    require(fallback.quality.entirely_replaced_by_bvh,"missing explicit fallback");
    require(fallback.quality.decision==WallDistanceDecision::FULL_BVH_FALLBACK,
            "incorrect fallback provenance");
    require(!fallback.raw.geometrically_exact_on_input_mesh,
            "fallback falsely certified original method");
    require(fallback.selected.geometrically_exact_on_input_mesh,
            "full fallback was not geometrically exact");
    require(fallback.selected.method=="search_based_bvh_full_fallback",
            "fallback is not labelled");
    require(fallback.quality.geometry_tolerance_violations>0,
            "original errors erased from diagnostics");
    for (std::size_t i=0;i<grid.points.size();++i) if (!grid.solid[i])
        require(std::abs(fallback.selected.distance[i]-safe.selected.distance[i])<1e-12,
                "fallback did not replace WHOLE field");

    // Sampled monitoring never qualifies a whole field, even if the sampled
    // points happen to be exact and an interior unsampled node is wrong.
    WallDistanceResult fake=safe.raw;
    fake.method="eikonal_synthetic_field";
    fake.active_pde_samples=1;
    fake.converged=true;
    std::size_t fluid_ordinal=0,corrupt=grid.points.size();
    for (std::size_t i=0;i<grid.points.size();++i) if (!grid.solid[i]) {
        if (fluid_ordinal==1) { corrupt=i; break; }
        ++fluid_ordinal;
    }
    require(corrupt<grid.points.size(),"no unsampled node");
    fake.distance[corrupt]+=0.5;
    WallDistanceServiceOptions monitor=eik;
    monitor.on_failure=WallDistanceFailureAction::REJECT;
    monitor.check_scope=WallDistanceCheckScope::SAMPLED_MONITORING;
    monitor.sampling_stride=3;
    const auto sampled=assess_wall_distance_field(fake,WallDistanceMethod::EIKONAL,
                                                   surface,grid,monitor);
    require(sampled.quality.geometry_checked_nodes<sampled.quality.fluid_nodes,
            "monitor was not sampled");
    require(sampled.quality.geometry_tolerance_violations==0,
            "corrupt node unexpectedly sampled");
    require(sampled.quality.decision==WallDistanceDecision::MONITORING_ONLY,
            "monitoring falsely certified output");
    require(!sampled.quality.approved_for_euclidean_wall_distance,
            "sampled monitoring approved SST field");
    monitor.check_scope=WallDistanceCheckScope::ALL_FLUID_NODES;
    const auto exhaustive=assess_wall_distance_field(fake,WallDistanceMethod::EIKONAL,
                                                      surface,grid,monitor);
    require(exhaustive.quality.geometry_tolerance_violations==1,
            "full check missed injected geometry error");
    require(exhaustive.quality.decision==WallDistanceDecision::REJECTED_GEOMETRIC_ERROR,
            "wrong full-check verdict");

    // Non-Euclidean directional heuristics must not be accidentally promoted
    // to y+ / SST even if they match BVH at every tested point by coincidence.
    WallDistanceServiceOptions directional;
    directional.method=WallDistanceMethod::DIRECTIONAL_MESH_WAVE;
    auto directional_fake=safe.raw;
    directional_fake.method="directional_synthetic";
    directional_fake.converged=true;
    directional_fake.geometrically_exact_on_input_mesh=false;
    const auto mismatch=assess_wall_distance_field(
        directional_fake,directional.method,surface,grid,directional);
    require(mismatch.quality.geometry_tolerance_violations==0,"bad synthetic reference");
    require(mismatch.quality.decision==WallDistanceDecision::REJECTED_METRIC_MISMATCH,
            "directional metric was promoted to Euclidean physics");

    // Full geometry agreement cannot assert PDE qualification if no PDE node
    // actually participated or the raw iteration did not converge.
    eik.on_failure=WallDistanceFailureAction::REJECT;
    auto only_seeds=safe.raw;
    only_seeds.geometrically_exact_on_input_mesh=false;
    only_seeds.active_pde_samples=0;
    only_seeds.converged=true;
    const auto unsolved=assess_wall_distance_field(only_seeds,eik.method,surface,grid,eik);
    require(unsolved.quality.geometry_within_tolerance,"synthetic exact seed wrong");
    require(!unsolved.quality.solver_actually_exercised,"empty PDE mislabeled as solved");
    require(unsolved.quality.decision==WallDistanceDecision::REJECTED_NUMERICAL_FAILURE,
            "all-BVH seeds certified an Eikonal solve");

    only_seeds.active_pde_samples=1;
    only_seeds.converged=false;
    const auto nonconverged=assess_wall_distance_field(only_seeds,eik.method,surface,grid,eik);
    require(nonconverged.quality.decision==WallDistanceDecision::REJECTED_NUMERICAL_FAILURE,
            "low geometry error hid numerical failure");

    // A corrupted UNSAMPLED value must still fail field-integrity validation.
    fake.distance[corrupt]=std::numeric_limits<double>::quiet_NaN();
    monitor.check_scope=WallDistanceCheckScope::SAMPLED_MONITORING;
    const auto nonfinite=assess_wall_distance_field(fake,eik.method,surface,grid,monitor);
    require(nonfinite.quality.invalid_fluid_nodes==1,"non-finite node was overlooked");
    require(nonfinite.quality.decision==WallDistanceDecision::REJECTED_INVALID_FIELD,
            "non-finite unsampled node was accepted");

    // Explicitly unsupported inputs must fail closed, not turn into an oracle.
    auto invalid=eik;
    invalid.relative_tolerance=std::numeric_limits<double>::infinity();
    threw=false;
    try { (void)audit_wall_distance_service(surface,grid,invalid); }
    catch (const std::invalid_argument&) { threw=true; }
    require(threw,"infinite tolerance accepted");
    invalid=eik;
    invalid.check_scope=WallDistanceCheckScope::SAMPLED_MONITORING;
    invalid.sampling_stride=0;
    threw=false;
    try { (void)audit_wall_distance_service(surface,grid,invalid); }
    catch (const std::invalid_argument&) { threw=true; }
    require(threw,"zero sampling stride accepted");
    invalid=eik;
    invalid.check_scope=WallDistanceCheckScope::NONE;
    const auto unchecked=assess_wall_distance_field(only_seeds,eik.method,surface,grid,invalid);
    require(!unchecked.quality.approved_for_euclidean_wall_distance,
            "unchecked field certified");

    // The dispatcher must reject a damaged grid BEFORE a graph/PDE method
    // can dereference a neighbour or index outside its arrays.
    auto damaged=grid;
    damaged.neighbours[0].push_back(grid.points.size()+1);
    threw=false;
    try { (void)audit_wall_distance_service(surface,damaged,eik); }
    catch (const std::invalid_argument&) { threw=true; }
    require(threw,"invalid grid neighbour reached the PDE dispatcher");
    damaged=grid;
    damaged.nx+=1;
    threw=false;
    try { (void)audit_wall_distance_service(surface,damaged,eik); }
    catch (const std::invalid_argument&) { threw=true; }
    require(threw,"wrong grid dimensions reached the PDE dispatcher");

    std::cout<<"PASS: exact default, raw Eikonal preserved, full BVH fallback, "
                "sampled diagnostics not certified, numerical & metric gates\n";
    return 0;
}
