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
        std::cerr << "wall-distance method: " << wall_distance_method_name(method) << "\\n";
        require(r.distance.size()==g.points.size(),"distance size mismatch");
        for(std::size_t i=0;i<r.distance.size();++i)
            if(!g.solid[i]) require(r.valid[i] && std::isfinite(r.distance[i]) && r.distance[i]>=0.0,
                                    "invalid wall distance");
    }

    const WallDistanceBvh bvh(s);
    require(bvh.triangle_count()==s.triangles.size(),"BVH triangle count mismatch");
    require(bvh.node_count()>0,"BVH is empty");
    const std::array<WallDistanceVec3,8> probes={{{-0.2,0.4,0.4},{0.5,0.5,-0.3},
        {1.2,0.2,0.7},{0.25,0.8,1.4},{-0.4,-0.2,0.3},{0.7,1.3,0.1},
        {1.4,1.4,1.4},{0.12,0.34,0.91}}};
    for(const auto& p:probes) {
        const double brute=exact_point_distance(s,p);
        const double accelerated=bvh.nearest_distance(p);
        require(std::abs(brute-accelerated)<1e-12,"BVH nearest distance mismatch");
    }
    const auto search=compute_wall_distance(WallDistanceMethod::SEARCH_BASED,s,g,20);
    for(std::size_t i=0;i<g.points.size();++i) if(!g.solid[i]) {
        const double brute=ref[i];
        if(brute<=1.0) require(std::abs(search.distance[i]-brute)<1e-12,
                               "search-based near-wall distance mismatch");
    }

    const auto mesh=compute_wall_distance(WallDistanceMethod::MESH_WAVE,s,g,40);
    const auto directional=compute_wall_distance(WallDistanceMethod::DIRECTIONAL_MESH_WAVE,s,g,40);
    const auto eik=compute_wall_distance(WallDistanceMethod::EIKONAL,s,g,40);
    const auto hj=compute_wall_distance(WallDistanceMethod::HAMILTON_JACOBI,s,g,40);
    const auto adv=compute_wall_distance(WallDistanceMethod::ADVECTION_DIFFUSION,s,g,40);
    bool directional_differs=false, hj_differs=false, adv_differs=false;
    for(std::size_t i=0;i<g.points.size();++i) if(!g.solid[i]) {
        directional_differs |= std::abs(mesh.distance[i]-directional.distance[i])>1e-13;
        hj_differs |= std::abs(eik.distance[i]-hj.distance[i])>1e-13;
        adv_differs |= std::abs(eik.distance[i]-adv.distance[i])>1e-13;
    }
    require(directional_differs,"directional mesh wave collapsed to isotropic mesh wave");
    require(hj_differs,"Hamilton-Jacobi collapsed to Eikonal");
    require(adv_differs,"advection-diffusion collapsed to Eikonal");

    // Canonical analytical test: an infinite planar wall at x=0 has
    // d(x,y,z)=x for x>=0, |grad d|=1 and laplacian(d)=0. This isolates
    // the PDE operators from geometry/corner errors.
    WallSurface plane;
    plane.points={{0,-1,-1},{0,1,-1},{0,1,1},{0,-1,1}};
    plane.triangles={{{0,1,2}},{{0,2,3}}};
    const auto pg=make_wall_distance_grid(
        9,5,5,{0,-1,-1},{0.25,0.5,0.5},
        [](const WallDistanceVec3& p){ return p.x<0.0; });
    const std::size_t pc=pg.index(4,2,2);
    std::vector<double> linear(pg.points.size(),0.0);
    for(std::size_t id=0;id<linear.size();++id) linear[id]=pg.points[id].x;
    require(std::abs(godunov_gradient_at(linear,pg,pc)-1.0)<1e-13,
            "Godunov gradient is not exact for d=x");
    const std::array<WallDistanceMethod,4> plane_methods={{
        WallDistanceMethod::EIKONAL,
        WallDistanceMethod::HAMILTON_JACOBI,
        WallDistanceMethod::ADVECTION_DIFFUSION,
        WallDistanceMethod::HYBRID_POISSON_EIKONAL}};
    for(const auto method:plane_methods) {
        const auto r=compute_wall_distance(method,plane,pg,120);
        for(std::size_t id=0;id<r.distance.size();++id) {
            const std::size_t k=id/(pg.nx*pg.ny), rem=id%(pg.nx*pg.ny), j=rem/pg.nx, i=rem%pg.nx;
            // The analytical d=x solution is used only where the Cartesian
            // stencil is complete. The outer boundary has a separate
            // numerical boundary condition and is not part of this operator
            // consistency test.
            if(pg.solid[id] || i==0 || i+1==pg.nx || j==0 || j+1==pg.ny || k==0 || k+1==pg.nz) continue;
            const double expected=pg.points[id].x;
            require(r.valid[id] && std::isfinite(r.distance[id]),
                    "analytical plane produced an invalid distance");
            require(std::abs(r.distance[id]-expected)<2e-10,
                    "wall-distance PDE failed the analytical planar-wall solution");
        }
    }

    // Poisson operator check independent of the seed-band initialization:
    // phi(x)=x*L-x^2/2 satisfies phi''=-1 for 0<x<L.
    const double L=2.0;
    std::vector<double> phi(pg.points.size(),0.0);
    std::vector<unsigned char> fixed(pg.points.size(),0);
    for(std::size_t id=0;id<phi.size();++id) {
        const double x=pg.points[id].x;
        phi[id]=x*L-0.5*x*x;
    }
    const double poisson_res=poisson_residual_inf(phi,pg,fixed);
    require(poisson_res<1e-13,
            "discrete Poisson operator is not exact for the analytical planar solution");

    // The hybrid method must solve its transport equation; it must not
    // reduce to an algebraic blend of two completed distance fields.
    const auto hybrid=compute_wall_distance(WallDistanceMethod::HYBRID_POISSON_EIKONAL,s,g,40);
    const auto poisson=compute_wall_distance(WallDistanceMethod::POISSON,s,g,40);
    const auto eikonal=compute_wall_distance(WallDistanceMethod::EIKONAL,s,g,40);
    bool hybrid_is_distinct=false;
    for(std::size_t i=0;i<g.points.size();++i) if(!g.solid[i]) {
        const double blended=0.35*poisson.distance[i]+0.65*eikonal.distance[i];
        if(std::abs(hybrid.distance[i]-blended)>1e-8) { hybrid_is_distinct=true; break; }
    }
    require(hybrid_is_distinct,"hybrid method is still an algebraic Poisson/Eikonal blend");

    const auto exact=compute_wall_distance(WallDistanceMethod::EXACT_GEOMETRIC,s,g,20);
    const auto m=compare_wall_distance(g,ref,exact.distance,0.5);
    require(m.l2_relative<1e-14,"exact benchmark regression failed");
    std::cout<<"wall-distance method contract: PASS\n";
    return 0;
}
