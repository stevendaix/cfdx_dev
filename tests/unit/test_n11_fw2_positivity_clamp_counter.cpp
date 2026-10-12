#include "cfdx/physics/wall_distance.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace cfdx::physics;

static void require(bool ok,const char* msg) {
    if(!ok) throw std::runtime_error(msg);
}

// N11 backlog F-W2 (#480): the pseudo-time positivity clamp
// nd = max(0, candidate) of the explicit transport steppers
// (Hamilton-Jacobi and hybrid Poisson+HJ) must stay behaviourally unchanged
// and become attributable through the positivity_clamp_activations counter
// (stepper out-parameter + WallDistanceResult field).
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
    const auto g=make_wall_distance_grid(
        9,9,9,{-1,-1,-1},{0.25,0.25,0.25},
        [](const WallDistanceVec3& p){ return p.x>0.0 && p.x<1.0 && p.y>0.0 && p.y<1.0 && p.z>0.0 && p.z<1.0; });

    // Contract 1: whatever the natural clamp activity is, the produced
    // distance field must stay finite and non-negative (the clamp does its
    // job) and the reported counter must be deterministic across runs.
    std::size_t hj_clamps=0;
    {
        const auto r1=compute_wall_distance(WallDistanceMethod::HAMILTON_JACOBI,s,g,20);
        const auto r2=compute_wall_distance(WallDistanceMethod::HAMILTON_JACOBI,s,g,20);
        require(r1.positivity_clamp_activations==r2.positivity_clamp_activations,
                "HJ positivity clamp counter is not deterministic");
        for(std::size_t i=0;i<r1.distance.size();++i) if(!g.solid[i])
            require(r1.valid[i] && std::isfinite(r1.distance[i]) && r1.distance[i]>=0.0,
                    "HJ produced a negative or invalid distance despite the positivity clamp");
        hj_clamps=r1.positivity_clamp_activations;
    }

    // Contract 2: the stepper out-parameter and the WallDistanceResult field
    // must agree (the dispatcher wires the counter, it does not invent it).
    {
        std::size_t used=0,clamps=0; double resid=0.0;
        const auto d=hamilton_jacobi_distance(s,g,20,0.25,0.7,&used,&resid,&clamps);
        require(clamps==hj_clamps,
                "stepper clamp counter disagrees with the WallDistanceResult wiring");
        require(d.size()==g.points.size(),
                "HJ stepper returned a wrong field size");
        const auto rh=compute_wall_distance(WallDistanceMethod::HYBRID_POISSON_EIKONAL,s,g,20);
        std::size_t hused=0,hclamps=0,hiter=0; double hresid=0.0,hpois=0.0;
        (void)hybrid_poisson_hamilton_jacobi_distance(s,g,20,0.25,0.9,&hused,&hresid,
            1.0,&hiter,&hpois,&hclamps);
        require(hclamps==rh.positivity_clamp_activations,
                "hybrid stepper clamp counter disagrees with the result wiring");
        for(std::size_t i=0;i<rh.distance.size();++i) if(!g.solid[i])
            require(rh.valid[i] && std::isfinite(rh.distance[i]) && rh.distance[i]>=0.0,
                    "hybrid Poisson+HJ produced a negative or invalid distance");
    }

    // Contract 3: monotonic budget -- a larger iteration budget executes a
    // superset of the deterministic pseudo-time steps, so the clamp counter
    // cannot decrease when the budget grows.
    {
        const auto small=compute_wall_distance(WallDistanceMethod::HAMILTON_JACOBI,s,g,20);
        const auto large=compute_wall_distance(WallDistanceMethod::HAMILTON_JACOBI,s,g,80);
        require(large.positivity_clamp_activations>=small.positivity_clamp_activations,
                "positivity clamp counter decreased with a larger iteration budget");
    }

    std::cout<<"N11 F-W2 positivity clamp reporting: PASS\n";
    return 0;
}
