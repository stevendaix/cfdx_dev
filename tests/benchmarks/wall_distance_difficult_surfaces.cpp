// Difficult wall geometries benchmark: no BVH corrections to method outputs.
// Independent analytic distances for a concave 270 degree fluid corner and tilted slits.
#include "cfdx/physics/wall_distance.h"
#include <array>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
using namespace cfdx::physics;
static void add_quad(WallSurface& s, WallDistanceVec3 a, WallDistanceVec3 b,
                     WallDistanceVec3 c, WallDistanceVec3 d) {
    const auto j=s.points.size();
    s.points.insert(s.points.end(),{a,b,c,d});
    s.triangles.push_back({{{j,j+1,j+2}}});
    s.triangles.push_back({{{j,j+2,j+3}}});
}
static double dot(WallDistanceVec3 a,WallDistanceVec3 b){return wd_dot(a,b);}
int main(int argc,char** argv) {
    std::cout<<"case,n,method,l2_percent,active_l2_percent,linf,fluid,active,invalid,converged,residual,seed_count,outer_count,active_pde,stopping_reason\n";
    for (const std::string name:{"concave_corner","tilted_slab","micro_gap"}) {
        WallSurface s;
        WallDistanceVec3 normal{std::cos(.63),std::sin(.63),0};
        WallDistanceVec3 tang{-normal.y,normal.x,0};
        const double w=name=="micro_gap"?.135:.4;
        if(name=="concave_corner") {
            // Solid occupies x<0,y<0. The fluid has a 270 deg reentrant corner.
            add_quad(s,{0,-3,-3},{0,0,-3},{0,0,3},{0,-3,3});
            add_quad(s,{-3,0,-3},{0,0,-3},{0,0,3},{-3,0,3});
        } else {
            for (double sign:{-1.,1.}) {
                const auto p=normal*(sign*w);
                add_quad(s,p+tang*(-3.)+WallDistanceVec3{0,0,-3},
                           p+tang*(3.)+WallDistanceVec3{0,0,-3},
                           p+tang*(3.)+WallDistanceVec3{0,0,3},
                           p+tang*(-3.)+WallDistanceVec3{0,0,3});
            }
        }
        const WallDistanceBvh bvh(s);
        for(std::size_t n:{std::size_t{9},std::size_t{17},std::size_t{25}}) {
            if(argc>1 && n!=static_cast<std::size_t>(std::atoi(argv[1]))) continue;
            const double h=1.8/double(n-1);
            auto solid=[&](WallDistanceVec3 p) {
                if(name=="concave_corner") return p.x < -1e-12 && p.y < -1e-12;
                return std::abs(dot(p,normal))>w;
            };
            auto exact=[&](WallDistanceVec3 p) {
                if(name!="concave_corner") return w-std::abs(dot(p,normal));
                if(p.x>=0 && p.y>=0) return std::hypot(p.x,p.y);
                if(p.x>=0) return p.x;
                if(p.y>=0) return p.y;
                return 0.;
            };
            const auto g=make_wall_distance_grid(n,n,n,{-0.9,-0.9,-0.9},{h,h,h},solid);
            std::vector<unsigned char> anchored(g.points.size(),0);
            for(const auto id:wall_seed_nodes(bvh,g,1.6*h)) anchored[id]=1;
            for(std::size_t k=0;k<n;++k)for(std::size_t j=0;j<n;++j)for(std::size_t i=0;i<n;++i)
                if(i==0 || j==0 || k==0 || i+1==n || j+1==n || k+1==n)
                    anchored[g.index(i,j,k)]=1;
            for(auto method:{WallDistanceMethod::EXACT_GEOMETRIC,WallDistanceMethod::SEARCH_BASED,
                  WallDistanceMethod::MESH_WAVE,WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
                  WallDistanceMethod::POISSON,WallDistanceMethod::EIKONAL,
                  WallDistanceMethod::EIKONAL_SECOND_ORDER,WallDistanceMethod::HAMILTON_JACOBI,
                  WallDistanceMethod::ADVECTION_DIFFUSION,WallDistanceMethod::HYBRID_POISSON_EIKONAL}) {
                auto r=compute_wall_distance(method,s,g,1000);
                double sum=0,den=0,asum=0,aden=0,linf=0;
                std::size_t fluid=0,active=0,invalid=0;
                for(std::size_t id=0;id<g.points.size();++id) {
                    if(g.solid[id]) continue;
                    ++fluid;
                    const auto ref=exact(g.points[id]);
                    if(!std::isfinite(r.distance[id])) {++invalid;continue;}
                    auto err=r.distance[id]-ref;
                    sum+=err*err;den+=ref*ref;
                    linf=std::max(linf,std::abs(err));
                    if(!anchored[id]) {++active;asum+=err*err;aden+=ref*ref;}
                }
                std::cout<<name<<','<<n<<','<<wall_distance_method_name(method)<<','
                  <<std::setprecision(10)<<(invalid?std::numeric_limits<double>::infinity():100*std::sqrt(sum/std::max(den,1e-30)))<<','
                  <<(active>0?100*std::sqrt(asum/std::max(aden,1e-30)):std::numeric_limits<double>::quiet_NaN())<<','
                  <<linf<<','<<fluid<<','<<active<<','<<invalid<<','<<r.converged<<','<<r.residual_inf<<','
                  <<r.exact_wall_seed_samples<<','<<r.exact_outer_boundary_samples<<','
                  <<r.active_pde_samples<<','<<r.stopping_reason<<'\n';
            }
        }
    }
}
