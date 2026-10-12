#include "cfdx/physics/wall_distance.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace cfdx::physics;

static void require(bool ok,const char* msg) {
    if(!ok) throw std::runtime_error(msg);
}

// N11 backlog F-W1 (#480): the nearest_normal() {1,0,0} fallback for a
// degenerate (zero-area) nearest triangle must stay behaviourally unchanged
// and become attributable through the degenerate_normal_fallbacks() counter
// and the WallDistanceResult::nearest_normal_fallbacks field.
int main() {
    // Healthy cube surface (same patch as test_wall_distance_methods.cpp).
    WallSurface s;
    s.points={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
              {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    const std::array<std::array<std::size_t,4>,6> q={{
        {{0,3,2,1}},{{4,5,6,7}},{{0,1,5,4}},{{1,2,6,5}},{{2,3,7,6}},{{3,0,4,7}}}};
    for(const auto& f:q) {
        s.triangles.push_back({{f[0],f[1],f[2]}});
        s.triangles.push_back({{f[0],f[2],f[3]}});
    }

    // Contract 1: a healthy surface never increments the counter and
    // nearest_normal returns a unit geometric normal.
    {
        const WallDistanceBvh bvh(s);
        const std::array<WallDistanceVec3,4> probes={{{-0.2,0.4,0.4},{0.5,0.5,-0.3},
            {1.2,0.2,0.7},{0.12,0.34,0.91}}};
        for(const auto& p:probes) {
            const auto n=bvh.nearest_normal(p);
            require(std::abs(wd_norm(n)-1.0)<1e-13,
                    "healthy nearest_normal is not unit norm");
        }
        require(bvh.degenerate_normal_fallbacks()==0,
                "healthy surface must not report degenerate-normal fallbacks");
    }

    // Degenerate surface: a single collinear (zero-area) triangle. The cross
    // product is the zero vector, so nearest_normal must keep the historical
    // {1,0,0} fallback -- and must now count it (F-W1).
    WallSurface deg;
    deg.points={{0,0,0},{1,0,0},{2,0,0}};
    deg.triangles={{{0,1,2}}};
    {
        const WallDistanceBvh bvh(deg);
        const auto n=bvh.nearest_normal({0.5,0.25,0.0});
        require(n.x==1.0 && n.y==0.0 && n.z==0.0,
                "degenerate nearest_normal changed the historical fallback value");
        require(bvh.degenerate_normal_fallbacks()==1,
                "degenerate nearest_normal did not increment the fallback counter");
        (void)bvh.nearest_normal({0.5,0.25,0.0});
        require(bvh.degenerate_normal_fallbacks()==2,
                "fallback counter is not monotonic across calls");
        require(std::isfinite(bvh.nearest_distance({0.5,0.25,0.0})),
                "degenerate surface broke nearest_distance");
    }

    // Contract 2: WallDistanceResult carries the counter for the methods
    // that query the wall normal.
    const auto g=make_wall_distance_grid(
        9,9,9,{-1,-1,-1},{0.25,0.25,0.25},
        [](const WallDistanceVec3& p){ return p.x>0.0 && p.x<1.0 && p.y>0.0 && p.y<1.0 && p.z>0.0 && p.z<1.0; });
    const auto dir_healthy=compute_wall_distance(WallDistanceMethod::DIRECTIONAL_MESH_WAVE,s,g,20);
    require(dir_healthy.nearest_normal_fallbacks==0,
            "healthy directional mesh wave reported degenerate-normal fallbacks");
    const auto pois_healthy=compute_wall_distance(WallDistanceMethod::POISSON,s,g,20);
    require(pois_healthy.nearest_normal_fallbacks==0,
            "healthy Poisson solve reported degenerate-normal fallbacks");

    // Degenerate line-segment "wall" (x in [0,2] on the y=z=0 axis): every
    // nearest_normal() call of the directional metric falls back, so the
    // result counter must be strictly positive while the distance field
    // stays finite and valid (the fallback is counted, not fatal).
    const auto dg=make_wall_distance_grid(
        9,9,9,{-1,-1,-1},{0.5,0.5,0.5},
        [](const WallDistanceVec3& p){ return p.x<0.0; });
    const auto dir_deg=compute_wall_distance(WallDistanceMethod::DIRECTIONAL_MESH_WAVE,deg,dg,20);
    require(dir_deg.nearest_normal_fallbacks>0,
            "degenerate surface did not surface fallbacks in WallDistanceResult");
    for(std::size_t i=0;i<dir_deg.distance.size();++i) if(!dg.solid[i])
        require(dir_deg.valid[i] && std::isfinite(dir_deg.distance[i]) && dir_deg.distance[i]>=0.0,
                "directional mesh wave produced an invalid distance on the degenerate surface");

    std::cout<<"N11 F-W1 nearest-normal fallback reporting: PASS
";
    return 0;
}
