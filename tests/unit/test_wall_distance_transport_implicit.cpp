// Regression for transport solvers (raw PDE values, never BVH replacement).
#include "cfdx/physics/wall_distance.h"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
static WallSurface cube() {
    WallSurface s;
    s.points={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
              {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    const std::array<std::array<std::size_t,4>,6> faces={{{{0,3,2,1}},
       {{4,5,6,7}},{{0,1,5,4}},{{1,2,6,5}},{{2,3,7,6}},{{3,0,4,7}}}};
    for (auto f:faces) {
        s.triangles.push_back({{f[0],f[1],f[2]}});
        s.triangles.push_back({{f[0],f[2],f[3]}});
    }
    return s;
}
int main() {
    const auto s=cube();
    constexpr std::size_t n=33;
    const double h=2.0/static_cast<double>(n-1);
    const auto g=make_wall_distance_grid(n,n,n,{-0.5,-0.5,-0.5},{h,h,h},
        [](auto p){return p.x>0 && p.x<1 && p.y>0 && p.y<1 && p.z>0 && p.z<1;});
    const WallDistanceBvh bvh(s);
    for (const auto method:{WallDistanceMethod::ADVECTION_DIFFUSION,
                            WallDistanceMethod::HYBRID_POISSON_EIKONAL}) {
        const std::size_t budget=method==WallDistanceMethod::ADVECTION_DIFFUSION ? 42 : 120;
        auto r=compute_wall_distance(method,s,g,budget);
        if (r.geometrically_exact_on_input_mesh || r.active_pde_samples<100 ||
            r.iterations>budget || !r.converged || !(r.residual_inf<1e-8))
            throw std::runtime_error("transport PDE failed strict raw residual gate");
        double error2=0.0, reference2=0.0;
        for (std::size_t id=0;id<g.points.size();++id) {
            if(g.solid[id]) continue;
            const double ref=bvh.nearest_distance(g.points[id]);
            if (!(std::isfinite(r.distance[id]) && r.distance[id]>=0.0))
                throw std::runtime_error("nonfinite or negative transported distance");
            error2+=(r.distance[id]-ref)*(r.distance[id]-ref);
            reference2+=ref*ref;
        }
        const double l2=std::sqrt(error2/reference2);
        if (!(l2<0.08)) throw std::runtime_error("transport L2 geometry regression");
        std::cout<<r.method<<" iterations="<<r.iterations<<" residual="
                 <<r.residual_inf<<" l2="<<l2<<" active="<<r.active_pde_samples<<'\n';
    }
    // Genuine 3-D high-resolution stress: the old fixed Poisson-direction
    // blending stagnates near outer corners here, with O(1) PDE residual.
    // Retain strict true-operator residual; never approve by geometry L2 alone.
    {
        constexpr std::size_t fine_n=65, fine_budget=200;
        const double hh=2.0/static_cast<double>(fine_n-1);
        const auto fine=make_wall_distance_grid(fine_n,fine_n,fine_n,
            {-0.5,-0.5,-0.5},{hh,hh,hh},
            [](auto p){return p.x>0 && p.x<1 && p.y>0 && p.y<1 && p.z>0 && p.z<1;});
        const auto r=compute_wall_distance(WallDistanceMethod::HYBRID_POISSON_EIKONAL,
                                           s,fine,fine_budget);
        if(!r.converged || r.iterations>=fine_budget || r.residual_inf>=1e-8 ||
           r.active_pde_samples<1000 || r.geometrically_exact_on_input_mesh)
            throw std::runtime_error("65^3 corner hybrid nonconvergence regression");
        std::cout<<"hybrid_cube_65 iter="<<r.iterations<<" residual="
                 <<r.residual_inf<<" active="<<r.active_pde_samples<<'\n';
    }
}