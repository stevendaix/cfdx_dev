// Strict numerical tests on the two historically stagnating V52 geometries.
// The same WLS operator, coefficients, AMG and tolerance are used throughout.
#include "cfdx/physics/dom_holo_amg_global_v53.h"
#include "thermal_3d_mesh_v36.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
static void run(double depth,bool recover){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(8,8,10,.018,0.,depth);
    for(std::size_t i=0;i<mesh.points().size();++i){
        const double x=mesh.points().x_data()[i],y=mesh.points().y_data()[i],z=mesh.points().z_data()[i];
        const double xr=std::cos(.41)*x-std::sin(.41)*y;
        const double yr=std::sin(.41)*x+std::cos(.41)*y;
        mesh.points().set(i,std::cos(.33)*xr+std::sin(.33)*z,yr,
                          -std::sin(.33)*xr+std::cos(.33)*z);
    }
    const auto geo=build_fv_geometry(mesh);
    const auto dirs=make_dom_tensor_gauss_legendre_v39(4,16);
    DomWallBoundaryConditions walls;
    for(std::size_t i=0;i<mesh.boundary().n_patches();++i){
        const auto name=mesh.boundary().patch(i).name;
        walls[name]={name=="bottom"?0.:name=="west"?.45:.85,
                     name=="east"?515.:name=="west"?385.:420.};
    }
    const auto n=mesh.n_cells();
    std::vector<DomSmmSecondMomentV48::Tensor> E(n);
    std::vector<double> ka(n),kt(n),ib(n);
    for(std::size_t i=0;i<n;++i){
        const double x=geo.cell_centres[i].x;
        E[i]={1./3.+.009*x,1./3.-.009*x,1./3.,.014,.006,.007};
        ka[i]=.3+.02*std::abs(x);
        kt[i]=3.+.1*std::abs(x);
        ib[i]=blackbody_intensity(420.+10.*x);
    }
    DomSmmGrayWallV49 lo(mesh,geo,dirs,walls,E,kt,ka);
    DomSmmAnisotropicV52Options historical;
    const auto initial=solve_smm_anisotropic_v52(mesh,geo,lo,E,kt,ka,ib,historical);
    if(initial.linear_strict_pass || !initial.accepted_J.empty())
        throw std::runtime_error("V53: historical strict negative baseline changed");
    DomV53Options options;
    options.max_krylov=300;
    options.wide_budget=450;
    options.wide_restart=120;
    const auto improved=solve_smm_amg_corrected_v53(lo,lo.rhs(ib),ib,options);
    if(improved.strict!=recover || improved.J.empty()==recover)
        throw std::runtime_error("V53: strict original-WLS gate or fail-closed contract broken");
    std::cout<<"V53_EXTREME depth="<<depth
             <<" historical_L2="<<initial.true_linear_l2
             <<" improved_L2="<<improved.true_l2
             <<" strict_pass="<<improved.strict
             <<" wide_used="<<improved.wide_used
             <<" iterations="<<improved.iterations<<'\n';
}
int main(){
    try{run(.006,true);run(.003,false);std::cout<<"V53_EXTREME_PASS\n";}
    catch(const std::exception& e){std::cerr<<"V53_EXTREME_FAIL "<<e.what()<<'\n';return 1;}
}
