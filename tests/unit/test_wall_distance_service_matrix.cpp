#include "cfdx/physics/wall_distance_service.h"
#include <array>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>

using namespace cfdx::physics;

static void check(bool c,const char* message) {
    if (!c) throw std::runtime_error(message);
}

int main() {
    WallSurface s;
    s.points={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
              {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    constexpr std::array<std::array<std::size_t,4>,6> faces={{{{0,1,2,3}},
        {{4,5,6,7}},{{0,1,5,4}},{{1,2,6,5}},{{2,3,7,6}},{{3,0,4,7}}}};
    for (const auto& f:faces) {
        s.triangles.push_back({{f[0],f[1],f[2]}});
        s.triangles.push_back({{f[0],f[2],f[3]}});
    }
    const auto g=make_wall_distance_grid(13,13,13,{-0.5,-0.5,-0.5},
        {1.0/6,1.0/6,1.0/6},[](WallDistanceVec3 p) {
            return p.x>0 && p.x<1 && p.y>0 && p.y<1 && p.z>0 && p.z<1;
        });
    constexpr std::array<WallDistanceMethod,10> methods={{
        WallDistanceMethod::EXACT_GEOMETRIC,WallDistanceMethod::SEARCH_BASED,
        WallDistanceMethod::MESH_WAVE,WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
        WallDistanceMethod::POISSON,WallDistanceMethod::EIKONAL,
        WallDistanceMethod::EIKONAL_SECOND_ORDER,WallDistanceMethod::HAMILTON_JACOBI,
        WallDistanceMethod::ADVECTION_DIFFUSION,
        WallDistanceMethod::HYBRID_POISSON_EIKONAL
    }};
    std::cout<<"method,family,pde_active,solver_converged,checked_nodes,"
                "violations,raw_l2_relative,decision,approved,full_bvh_fallback\n";
    for (const auto method:methods) {
        WallDistanceServiceOptions options;
        options.method=method;
        options.max_iterations=500;
        options.absolute_tolerance=1e-12;
        options.relative_tolerance=1e-4;
        const auto r=audit_wall_distance_service(s,g,options);
        check(r.quality.source_method==wall_distance_method_name(method),
              "method selector did not preserve requested method");
        check(r.raw.distance==r.selected.distance,"REJECT mode silently replaced raw field");
        check(r.quality.geometry_checked_nodes==r.quality.fluid_nodes,
              "full control skipped fluid nodes");
        check(!r.quality.entirely_replaced_by_bvh,"unexpected exact fallback");
        if (method==WallDistanceMethod::EXACT_GEOMETRIC ||
            method==WallDistanceMethod::SEARCH_BASED)
            check(r.quality.approved_for_euclidean_wall_distance,
                  "geometric reference rejected");
        if (method==WallDistanceMethod::DIRECTIONAL_MESH_WAVE)
            check(!r.quality.approved_for_euclidean_wall_distance,
                  "directional metric incorrectly qualified");
        if (r.quality.approved_for_euclidean_wall_distance)
            check(r.quality.full_geometry_check && r.quality.geometry_within_tolerance,
                  "unverified method approved");
        std::cout<<r.quality.source_method<<','<<r.quality.method_family<<','
                 <<r.quality.active_pde_nodes<<','<<r.quality.solver_converged<<','
                 <<r.quality.geometry_checked_nodes<<','
                 <<r.quality.geometry_tolerance_violations<<','
                 <<std::setprecision(10)<<r.quality.relative_l2_error_checked<<','
                 <<static_cast<int>(r.quality.decision)<<','
                 <<r.quality.approved_for_euclidean_wall_distance<<','
                 <<r.quality.entirely_replaced_by_bvh<<'\n';
    }
}
