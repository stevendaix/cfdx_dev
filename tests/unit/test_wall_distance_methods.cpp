#include "cfdx/physics/wall_distance.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace cfdx::physics;

static void require(bool ok,const char* msg) {
    if(!ok) throw std::runtime_error(msg);
}

int main() {
    WallSurface s;
    s.points={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
              {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    const std::array<std::array<std::size_t,4>,6> q={{
        {{0,3,2,1}},{{4,5,6,7}},{{0,1,5,4}},{{1,2,6,5}},{{2,3,7,6}},{{3,0,4,7}}}};
    for(const auto& f:q) {
        s.triangles.push_back({{f[0],f[1],f[2]}});
        s.triangles.push_back({{f[0],f[2],f[3]}});
    }

    require(std::abs(exact_point_distance(s,{0.5,0.5,0.5})-0.5)<1e-12,
            "exact point-to-triangle distance is wrong");

    const auto g=make_wall_distance_grid(
        9,9,9,{-1,-1,-1},{0.25,0.25,0.25},
        [](const WallDistanceVec3& p){ return p.x>0.0 && p.x<1.0 && p.y>0.0 && p.y<1.0 && p.z>0.0 && p.z<1.0; });
    const auto ref=exact_reference(s,g);
    const std::array<WallDistanceMethod,9> methods={{
        WallDistanceMethod::EXACT_GEOMETRIC,
        WallDistanceMethod::SEARCH_BASED,
        WallDistanceMethod::MESH_WAVE,
        WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
        WallDistanceMethod::POISSON,
        WallDistanceMethod::EIKONAL,
        WallDistanceMethod::HAMILTON_JACOBI,
        WallDistanceMethod::ADVECTION_DIFFUSION,
        WallDistanceMethod::HYBRID_POISSON_EIKONAL}};
    for(const auto method:methods) {
        const auto r=compute_wall_distance(method,s,g,20);
        require(r.distance.size()==g.points.size(),"distance size mismatch");
        for(std::size_t i=0;i<r.distance.size();++i)
            if(!g.solid[i]) require(r.valid[i] && std::isfinite(r.distance[i]) && r.distance[i]>=0.0,
                                    "invalid wall distance");
    }
    const auto exact=compute_wall_distance(WallDistanceMethod::EXACT_GEOMETRIC,s,g,20);
    const auto m=compare_wall_distance(g,ref,exact.distance,0.5);
    require(m.l2_relative<1e-14,"exact benchmark regression failed");
    std::cout<<"wall-distance method contract: PASS\n";
    return 0;
}
