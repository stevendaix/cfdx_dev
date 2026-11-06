// Regression audit on hard wall geometry. No accuracy replacement by BVH.
#include "cfdx/physics/wall_distance.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
static void check(bool ok,const char* msg){if(!ok)throw std::runtime_error(msg);}
static void quad(WallSurface& s,WallDistanceVec3 a,WallDistanceVec3 b,
                 WallDistanceVec3 c,WallDistanceVec3 d){
    auto i=s.points.size();s.points.insert(s.points.end(),{a,b,c,d});
    s.triangles.push_back({{{i,i+1,i+2}}});s.triangles.push_back({{{i,i+2,i+3}}});
}
int main(){
    // A thin, rotated slab: the reference is analytic distance to planes.
    const auto nrm=WallDistanceVec3{std::cos(.63),std::sin(.63),0};
    const auto t=WallDistanceVec3{-nrm.y,nrm.x,0};
    const double half_width=.135;
    WallSurface slab;
    for(double sign:{-1.,1.}) {
        const auto p=nrm*(sign*half_width);
        quad(slab,p+t*(-3.)+WallDistanceVec3{0,0,-3},
                  p+t*3.+WallDistanceVec3{0,0,-3},
                  p+t*3.+WallDistanceVec3{0,0,3},
                  p+t*(-3.)+WallDistanceVec3{0,0,3});
    }
    const std::size_t n=17;
    const double h=1.8/double(n-1);
    const auto grid=make_wall_distance_grid(n,n,n,{-.9,-.9,-.9},{h,h,h},
       [&](auto p){return std::abs(wd_dot(p,nrm))>half_width;});
    const WallDistanceBvh bvh(slab);
    std::size_t fluid=0;
    for(std::size_t id=0;id<grid.points.size();++id)if(!grid.solid[id])++fluid;
    const auto seeds=wall_seed_nodes(bvh,grid,1.6*h);
    check(fluid>0 && seeds.size()==fluid,"test setup must have all points BVH seeded");
    for(auto m:{WallDistanceMethod::EIKONAL,WallDistanceMethod::EIKONAL_SECOND_ORDER,
                WallDistanceMethod::HAMILTON_JACOBI,WallDistanceMethod::ADVECTION_DIFFUSION,
                WallDistanceMethod::HYBRID_POISSON_EIKONAL}) {
        const auto r=compute_wall_distance(m,slab,grid,100);
        check(r.converged && r.geometrically_exact_on_input_mesh &&
              r.convergence_status==WallDistanceConvergenceStatus::DIRECT &&
              r.stopping_reason=="all_fluid_nodes_exact_bvh_seeds_or_outer_no_pde_solve" &&
              r.residual_inf==0.0 && r.active_pde_samples==0 && r.exact_wall_seed_samples==fluid,
              "all-BVH-seeded grid must be diagnosed as direct geometric data, not PDE convergence");
        for(std::size_t id=0;id<grid.points.size();++id)if(!grid.solid[id])
           check(std::abs(r.distance[id] - (half_width-std::abs(wd_dot(grid.points[id],nrm))))<1e-12,
                 "thin gap analytic distance mismatch");
    }
    // Nontrivial concave geometry: the higher-order solver can stagnate at
    // a corner. In that case it MUST NOT claim convergence to the PDE.
    WallSurface corner;
    quad(corner,{0,-3,-3},{0,0,-3},{0,0,3},{0,-3,3});
    quad(corner,{-3,0,-3},{0,0,-3},{0,0,3},{-3,0,3});
    const std::size_t nc=17;const double hc=1.8/double(nc-1);
    const auto cg=make_wall_distance_grid(nc,nc,nc,{-.9,-.9,-.9},{hc,hc,hc},
             [](auto p){return p.x < -1e-12 && p.y < -1e-12;});
    const auto rr=compute_wall_distance(WallDistanceMethod::EIKONAL_SECOND_ORDER,corner,cg,100);
    check(!rr.geometrically_exact_on_input_mesh && rr.exact_wall_seed_samples>0,
          "partly propagated corner field was incorrectly certified as exact BVH");
    check(rr.converged==(std::isfinite(rr.residual_inf)&&rr.residual_inf<1e-6),
          "corner residual and convergence flag disagree");
    std::cout<<"Difficult-wall tests passed: no PDE credit for BVH-only seed band; corner fail-closed\n";
}
