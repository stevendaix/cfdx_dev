#include "cfdx/physics/wall_distance_precision.h"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cfdx::physics;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

static WallSurface make_cube() {
    WallSurface s;
    s.points={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
              {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    const std::array<std::array<std::size_t,4>,6> faces = {{
        {{0,3,2,1}},{{4,5,6,7}},{{0,1,5,4}},
        {{1,2,6,5}},{{2,3,7,6}},{{3,0,4,7}}}};
    for (const auto& f : faces) {
        s.triangles.push_back({{f[0],f[1],f[2]}});
        s.triangles.push_back({{f[0],f[2],f[3]}});
    }
    return s;
}

int main() {
    const auto surface=make_cube();
    const auto g=make_wall_distance_grid(
        17,17,17,{-0.5,-0.5,-0.5},{0.125,0.125,0.125},
        [](const WallDistanceVec3& p){
            return p.x>0 && p.x<1 && p.y>0 && p.y<1 && p.z>0 && p.z<1;
        });
    const auto reference=exact_reference(surface,g);
    // Independent analytical cube oracle: outside distance to the clamped
    // point on [0,1]^3; inside distance to the nearest of its six planes.
    // This verifies BVH+triangle distance against an implementation that
    // does not call point_triangle_distance2() or exact_reference().
    const WallDistanceBvh independent_bvh(surface);
    for (std::size_t id=0; id<g.points.size(); ++id) {
        const auto& p=g.points[id];
        const auto outside=[](double x) { return x<0.0 ? -x : (x>1.0 ? x-1.0 : 0.0); };
        const double dx=outside(p.x), dy=outside(p.y), dz=outside(p.z);
        const double outer=std::sqrt(dx*dx+dy*dy+dz*dz);
        const double analytic=outer>0.0 ? outer : std::min({
            p.x,1.0-p.x,p.y,1.0-p.y,p.z,1.0-p.z});
        require(std::abs(independent_bvh.nearest_distance(p)-analytic)<1e-12,
                "BVH does not match independent analytical cube distance");
    }
    const auto ordinary=compute_wall_distance(WallDistanceMethod::EIKONAL,surface,g,100);
    std::size_t sweeps=0;
    const auto improved=second_order_eikonal_fast_sweep(surface,g,100,&sweeps);
    const auto e1=compare_wall_distance(g,reference,ordinary.distance,0.5);
    const auto e2=compare_wall_distance(g,reference,improved,0.5);
    require(sweeps>0 && sweeps<100,"second-order Eikonal did not finish its sweeps");
    require(e2.invalid_samples==0,"second-order Eikonal produced invalid samples");
    require(e1.l2_relative>0.07 && e1.l2_relative<0.08,"baseline Eikonal reference has changed");
    require(e2.l2_relative<0.03,"second-order Eikonal cube accuracy regression");
    require(e2.l2_relative < 0.5*e1.l2_relative,"second order lost its measured advantage");

    // The accuracy gate must test EVERY fluid point against the true nearest
    // triangle, not against the scalar PDE residual or a sampled subset.
    const std::array<WallDistanceMethod,9> methods = {{
        WallDistanceMethod::EXACT_GEOMETRIC,
        WallDistanceMethod::SEARCH_BASED,
        WallDistanceMethod::MESH_WAVE,
        WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
        WallDistanceMethod::POISSON,
        WallDistanceMethod::EIKONAL,
        WallDistanceMethod::HAMILTON_JACOBI,
        WallDistanceMethod::ADVECTION_DIFFUSION,
        WallDistanceMethod::HYBRID_POISSON_EIKONAL}};
    for (const auto method : methods) {
        const auto r=compute_wall_distance_accuracy_checked(method,surface,g,1000,{0,1e-10,true});
        require(r.accuracy_qualified,"precision guard did not qualify output");
        require(r.fluid_samples>4000,"precision guard did not examine all fluid samples");
        require(r.corrected_samples==0 || r.corrected_samples==r.fluid_samples,
                "BVH fallback must replace the entire field, not selected nodes");
        require(r.accepted_linf_absolute <= 1e-10,"precision guard returned inaccurate output");
        require(r.accepted.converged,"corrected geometric output must be accepted");
        if (method == WallDistanceMethod::EXACT_GEOMETRIC ||
            method == WallDistanceMethod::SEARCH_BASED) {
            require(r.violations==0,"exact distance generated unexpected corrections");
        } else if (method == WallDistanceMethod::MESH_WAVE) {
            // A wall-point wave can be exact on a cube aligned with its grid:
            // do not REQUIRE an error merely to prove an algorithm is distinct.
            require(r.violations==0,"aligned cube wall-point propagation regressed");
        } else {
            require(r.violations>0,"approximate method was not independently exercised");
            require(r.entire_field_replaced_by_bvh && r.corrected_samples==r.fluid_samples &&
                    r.accepted.method=="search_based_bvh_full_fallback" &&
                    !r.original_algorithm_qualified,
                    "full BVH fallback or provenance was misreported");
        }
    }
    const auto fail=compute_wall_distance_accuracy_checked(
        WallDistanceMethod::POISSON,surface,g,1000,{0,1e-10,false});
    require(!fail.accuracy_qualified && !fail.accepted.converged && fail.violations>0,
            "bad Poisson field must fail closed without exact correction");
    require(fail.accepted.stopping_reason=="geometry_accuracy_unqualified",
            "failed precision gate must be reported separately");

    for (double tolerance : {
        -1.0, std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::quiet_NaN()}) {
        bool rejected=false;
        try { (void)compute_wall_distance_accuracy_checked(
                    WallDistanceMethod::EIKONAL,surface,g,20,{tolerance,1e-8,true}); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"invalid absolute tolerance accepted");
        rejected=false;
        try { (void)compute_wall_distance_accuracy_checked(
                    WallDistanceMethod::EIKONAL,surface,g,20,{0.0,tolerance,true}); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"invalid relative tolerance accepted");
    }
    // A wall absent from the grid must not be silently labelled as an
    // accurate second-order solution.
    const auto far=make_wall_distance_grid(5,5,5,{100,100,100},{1,1,1},
                                         [](const WallDistanceVec3&){return false;});
    bool rejected=false;
    try { (void)second_order_eikonal_fast_sweep(surface,far,20); }
    catch (const std::invalid_argument&) { rejected=true; }
    require(rejected,"second-order Eikonal accepted missing wall seeds");

    std::cout<<"PASS: second-order Eikonal, 9 exhaustive BVH precision gates, fail-closed and invalid parameters\n";
}
