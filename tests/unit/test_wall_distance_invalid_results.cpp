#include "cfdx/physics/wall_distance.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using namespace cfdx::physics;

int main() {
    int failures=0, checks=0;
    const auto expect=[&](bool ok,const char* label) {
        ++checks;
        if(!ok) { ++failures; std::cerr << "FAIL: " << label << '\n'; }
    };
    const auto grid=make_wall_distance_grid(4,4,4,{10,0,0},{1,1,1},
        [](const WallDistanceVec3&){ return false; });
    std::vector<double> reference(grid.points.size(),1.0), distance=reference;
    auto metrics=compare_wall_distance(grid,reference,distance,2.0);
    expect(metrics.samples==grid.points.size() && metrics.l2_relative==0.0,
           "valid identical fields retain zero error");
    distance[0]=std::numeric_limits<double>::quiet_NaN();
    metrics=compare_wall_distance(grid,reference,distance,2.0);
    expect(std::isinf(metrics.l2_relative) && std::isinf(metrics.linf_relative),
           "one NaN cannot be excluded from the accuracy verdict");
    distance.assign(grid.points.size(),std::numeric_limits<double>::infinity());
    metrics=compare_wall_distance(grid,reference,distance,2.0);
    expect(metrics.samples==0 && std::isinf(metrics.l2_relative) &&
           std::isinf(metrics.near_wall_l2_relative),
           "a wholly invalid field cannot have zero global or near-wall error");
    distance=reference;
    reference[1]=std::numeric_limits<double>::quiet_NaN();
    metrics=compare_wall_distance(grid,reference,distance,2.0);
    expect(std::isinf(metrics.l2_relative),"a corrupt oracle cannot certify accuracy");
    reference=distance;
    distance[2]=-1.0;
    metrics=compare_wall_distance(grid,reference,distance,2.0);
    expect(std::isinf(metrics.l2_relative),"negative distance is an invalid field");
    distance=reference;
    metrics=compare_wall_distance(grid,reference,distance,0.5);
    expect(std::isinf(metrics.near_wall_l2_relative),
           "an empty near-wall sample band cannot report zero error");

    WallSurface wall;
    wall.points={{0,-4,-4},{0,4,-4},{0,4,4},{0,-4,4}};
    wall.triangles={{{0,1,2}},{{0,2,3}}};
    // The wall is outside the grid and its seed band. Geometric distance is
    // well defined, but propagation methods cannot reach any of these nodes.
    for(const auto method:{WallDistanceMethod::MESH_WAVE,
                           WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
                           WallDistanceMethod::EIKONAL,
                           WallDistanceMethod::HAMILTON_JACOBI,
                           WallDistanceMethod::ADVECTION_DIFFUSION}) {
        const auto result=compute_wall_distance(method,wall,grid,20);
        expect(!result.converged && std::isinf(result.residual_inf),
               "a propagation result with no reachable sample cannot converge");
    }
    auto solid=grid;
    solid.solid.assign(grid.points.size(),1);
    const auto empty_result=compute_wall_distance(WallDistanceMethod::EXACT_GEOMETRIC,wall,solid,20);
    expect(!empty_result.converged,"an all-solid domain is not an accuracy certificate");
    metrics=compare_wall_distance(solid,reference,distance,2.0);
    expect(metrics.samples==0 && std::isinf(metrics.l2_relative),
           "an all-solid comparison has no accuracy evidence");

    distance.assign(grid.points.size(),1.0);
    distance[0]=std::numeric_limits<double>::quiet_NaN();
    std::vector<unsigned char> fixed(grid.points.size(),1);
    expect(std::isinf(wall_distance_pde_residual_inf(WallDistanceMethod::EIKONAL,distance,grid,&fixed)),
           "invalid fixed seeds invalidate the PDE residual");
    distance.assign(grid.points.size(),0.0);
    expect(wall_distance_pde_residual_inf(WallDistanceMethod::ADVECTION_DIFFUSION,distance,grid)>0.99,
           "zero gradient does not erase the unit source residual");
    std::cout << "checks=" << checks << " failures=" << failures << '\n';
    return failures==0 ? 0 : 1;
}
