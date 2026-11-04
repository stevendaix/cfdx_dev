// Independent smooth-sphere analytical oracle AND exact triangulated-surface
// reference. Raw CFDX methods only: never calls the precision guard.
// Build: c++ -std=c++17 -O2 -Isrc tests/benchmarks/wall_distance_curved_raw.cpp -o wall_distance_curved_raw
#include "cfdx/physics/wall_distance.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;

static constexpr double radius=0.42;
static const WallDistanceVec3 center{0.5,0.5,0.5};
static WallDistanceVec3 unit(WallDistanceVec3 q) {
    const auto norm=wd_norm(q);
    return q*(1.0/norm);
}
static WallSurface sphere(int levels) {
    WallSurface s;
    const std::array<WallDistanceVec3,6> v={{{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}}};
    const std::array<std::array<std::size_t,3>,8> faces={{{{4,0,2}},{{4,2,1}},{{4,1,3}},{{4,3,0}},{{5,2,0}},{{5,1,2}},{{5,3,1}},{{5,0,3}}}};
    std::vector<std::array<WallDistanceVec3,3>> triangles;
    for(auto f:faces) triangles.push_back({v[f[0]],v[f[1]],v[f[2]]});
    for(int lev=0;lev<levels;++lev) {
        std::vector<std::array<WallDistanceVec3,3>> next;
        next.reserve(triangles.size()*4);
        for(auto t:triangles) {
            const auto ab=unit((t[0]+t[1])*0.5);
            const auto bc=unit((t[1]+t[2])*0.5);
            const auto ca=unit((t[2]+t[0])*0.5);
            next.push_back({t[0],ab,ca});
            next.push_back({ab,t[1],bc});
            next.push_back({ca,bc,t[2]});
            next.push_back({ab,bc,ca});
        }
        triangles.swap(next);
    }
    for(auto t:triangles) {
        const auto first=s.points.size();
        for(auto q:t) s.points.push_back(center+q*radius);
        s.triangles.push_back({{{first,first+1,first+2}}});
    }
    return s;
}
int main() {
    std::cout<<"shape,n,method,l2_triangle_percent,l2_smooth_percent,active_l2_triangle_percent,triangle_linf,fluid,active,anchored,converged,residual\n";
    const auto s=sphere(3); // 512 planar facets; exact CAD sphere is intentionally DIFFERENT.
    const WallDistanceBvh bvh(s);
    for(const std::size_t n:{std::size_t{9},std::size_t{17},std::size_t{25}}) {
        const double h=2.0/static_cast<double>(n-1);
        const auto g=make_wall_distance_grid(n,n,n,{-0.5,-0.5,-0.5},{h,h,h},[](auto q){return wd_norm(q-center)<radius;});
        std::vector<unsigned char> anchored(g.points.size(),0);
        for(auto id:wall_seed_nodes(bvh,g,1.6*h)) anchored[id]=1;
        for(std::size_t k=0;k<g.nz;++k) for(std::size_t j=0;j<g.ny;++j) for(std::size_t i=0;i<g.nx;++i) {
            const auto id=g.index(i,j,k);
            if(i==0||j==0||k==0||i+1==g.nx||j+1==g.ny||k+1==g.nz) anchored[id]=1;
        }
        std::vector<double> reference(g.points.size(),0.0),smooth(g.points.size(),0.0);
        for(std::size_t i=0;i<g.points.size();++i) if(!g.solid[i]) {
            reference[i]=bvh.nearest_distance(g.points[i]);
            smooth[i]=std::max(0.0,wd_norm(g.points[i]-center)-radius);
        }
        for(auto m:{WallDistanceMethod::EXACT_GEOMETRIC,WallDistanceMethod::SEARCH_BASED,
                    WallDistanceMethod::MESH_WAVE,WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
                    WallDistanceMethod::POISSON,WallDistanceMethod::EIKONAL,
                    WallDistanceMethod::EIKONAL_SECOND_ORDER,WallDistanceMethod::HAMILTON_JACOBI,
                    WallDistanceMethod::ADVECTION_DIFFUSION,WallDistanceMethod::HYBRID_POISSON_EIKONAL}) {
            const auto r=compute_wall_distance(m,s,g,1000);
            double et2=0,es2=0,ea2=0,rt2=0,rs2=0,ra2=0,maxerr=0;
            std::size_t fluid=0,active=0;
            for(std::size_t i=0;i<g.points.size();++i) if(!g.solid[i]) {
                const double d=r.distance[i];
                if(!std::isfinite(d)) {std::cerr<<"NONFINITE: "<<wall_distance_method_name(m)<<" n="<<n<<"\n";return 1;}
                ++fluid;
                const double e=d-reference[i];
                et2+=e*e;rt2+=reference[i]*reference[i];
                const double se=d-smooth[i];es2+=se*se;rs2+=smooth[i]*smooth[i];
                maxerr=std::max(maxerr,std::abs(e));
                if(!anchored[i]) {++active;ea2+=e*e;ra2+=reference[i]*reference[i];}
            }
            const double l2_triangle=std::sqrt(et2/rt2);
            if(fluid==0 || active==0 || active>=fluid)
                throw std::runtime_error("invalid benchmark active/anchored partition");
            if(m==WallDistanceMethod::SEARCH_BASED || m==WallDistanceMethod::EXACT_GEOMETRIC) {
                if(l2_triangle>1e-12 || !r.geometrically_exact_on_input_mesh)
                    throw std::runtime_error("direct geometric distance differs from independent triangle BVH");
                if(std::sqrt(es2/rs2)>0.01)
                    throw std::runtime_error("unexpected triangulated-vs-smooth sphere geometry error");
            }
            if(m==WallDistanceMethod::EIKONAL_SECOND_ORDER) {
                if(!r.converged || r.residual_inf>1e-6 || l2_triangle>0.015)
                    throw std::runtime_error("second-order eikonal accuracy/residual regression on sphere");
            }
            if(m==WallDistanceMethod::MESH_WAVE && l2_triangle>0.03)
                throw std::runtime_error("wall-point propagation lost its curved-wall accuracy");
            std::cout<<"sphere_triangulated,"<<n<<","<<wall_distance_method_name(m)<<","<<std::setprecision(11)
                <<100*std::sqrt(et2/rt2)<<","<<100*std::sqrt(es2/rs2)<<","<<100*std::sqrt(ea2/ra2)
                <<","<<maxerr<<","<<fluid<<","<<active<<","<<fluid-active<<","<<r.converged<<","<<r.residual_inf<<"\n";
        }
    }
}
