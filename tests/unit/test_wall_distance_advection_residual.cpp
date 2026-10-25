#include "cfdx/physics/wall_distance.h"
#include <cmath>
#include <iostream>
using namespace cfdx::physics;
int main() {
    int failed=0;
    const auto expect=[&](bool ok,const char* message){if(!ok){++failed;std::cerr<<message<<'\n';}};
    // Manufactured stationary ridge for the actual discrete transport rule.
    // Opposing derivatives cancel in U. With h=1/4, gamma=1/20 and
    // d_E=d_W=1/4, gamma*(2*d_P-d_E-d_W)/h^2=1 gives d_P=7/8.
    auto g=make_wall_distance_grid(5,3,3,{0,-1,-1},{.25,1,1},
        [](const WallDistanceVec3& p){return p.x<=0 || p.x>=1;});
    std::vector<double> d(g.points.size(),0.0);
    std::vector<unsigned char> fixed(g.points.size(),1);
    for(std::size_t id=0;id<d.size();++id)
        if(!g.solid[id]) d[id]=g.points[id].x==.5 ? .875 : .25;
    const auto centre=g.index(2,1,1);fixed[centre]=0;
    const double stationary=wall_distance_pde_residual_inf(
        WallDistanceMethod::ADVECTION_DIFFUSION,d,g,&fixed,.05);
    expect(stationary<1e-13,"stationary transport row was measured with another operator");
    // Zero field must retain the unit source even when U is undefined/zero.
    std::fill(d.begin(),d.end(),0.0);
    expect(std::abs(wall_distance_pde_residual_inf(
        WallDistanceMethod::ADVECTION_DIFFUSION,d,g,&fixed,.05)-1.0)<1e-13,
        "zero gradient erased the unit source");

    // Fractional wall at x=.1: at P=.25, d=(x-.1)^2 gives d_P=.0225,
    // d_E=.16, delta_W=.15. U=1, transport=.15 and the assembled diffusion
    // Laplacian is -.0225/(.25*.15)+(.16-.0225)/.25^2=1.6.
    // Hence |.15-1-.05*1.6|=.93, versus .95 with Shortley-Weller.
    WallSurface wall;wall.points={{.1,-4,-4},{.1,4,-4},{.1,4,4},{.1,-4,4}};
    wall.triangles={{{0,1,2}},{{0,2,3}}};WallDistanceBvh bvh(wall);
    auto cut=make_wall_distance_grid(5,3,3,{0,-1,-1},{.25,1,1},
        [](const WallDistanceVec3& p){return p.x<.1;});
    std::vector<double> quadratic(cut.points.size(),0.0);
    std::vector<unsigned char> mask(cut.points.size(),1);
    for(std::size_t id=0;id<quadratic.size();++id)
        if(!cut.solid[id])quadratic[id]=std::pow(cut.points[id].x-.1,2);
    mask[cut.index(1,1,1)]=0;
    const double fractional=wall_distance_pde_residual_inf(
        WallDistanceMethod::ADVECTION_DIFFUSION,quadratic,cut,&mask,.05,&bvh);
    expect(std::abs(fractional-.93)<1e-13,"fractional-wall residual changed the assembled diffusion stencil");
    // A stationary fractional-cut row with nonzero U also distinguishes the
    // operators without relying on the old zero-gradient skip: d_P=.2,
    // d_E=.95 gives transport=4/3 and the assembled Laplacian=20/3.
    // Thus transport-gamma*lap=1, while the S-W residual is 1/12.
    for(std::size_t id=0;id<quadratic.size();++id)
        if(!cut.solid[id])quadratic[id]=.2+3.0*(cut.points[id].x-.25);
    const double stationary_cut=wall_distance_pde_residual_inf(
        WallDistanceMethod::ADVECTION_DIFFUSION,quadratic,cut,&mask,.05,&bvh);
    expect(stationary_cut<1e-13,"stationary cut row was certified with a different diffusion operator");
    std::cout<<"stationary_residual="<<stationary<<" fractional_residual="<<fractional
             <<" stationary_cut_residual="<<stationary_cut<<" failed="<<failed<<'\n';
    return failed==0?0:1;
}
