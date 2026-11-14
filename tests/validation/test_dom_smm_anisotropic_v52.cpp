#include "cfdx/physics/dom_smm_anisotropic_gate_v52.h"
#include "thermal_3d_mesh_v36.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
void require(bool good,const char* msg){if(!good)throw std::runtime_error(msg);}
Mesh grid(double depth){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(8,8,10,.018,0.,depth);
    for(std::size_t i=0;i<mesh.points().size();++i){
        const double x=mesh.points().x_data()[i],y=mesh.points().y_data()[i];
        const double z=mesh.points().z_data()[i];
        const double xr=std::cos(.41)*x-std::sin(.41)*y;
        const double yr=std::sin(.41)*x+std::cos(.41)*y;
        mesh.points().set(i,std::cos(.33)*xr+std::sin(.33)*z,
                            yr,-std::sin(.33)*xr+std::cos(.33)*z);
    }
    return mesh;
}
void check_case(double depth,bool strictly_qualified){
    auto mesh=grid(depth);const auto geo=build_fv_geometry(mesh);
    const auto nc=mesh.n_cells();require(nc==640,"V52: stress grid cells");
    const auto directions=make_dom_tensor_gauss_legendre_v39(4,16);
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
        const auto name=mesh.boundary().patch(p).name;
        walls[name]={name=="bottom"?0.:name=="west"?.45:.85,
            name=="east"?515.:name=="west"?385.:420.};
    }
    using Tensor=DomSmmSecondMomentV48::Tensor;
    std::vector<Tensor> E(nc);
    std::vector<double> ka(nc),kt(nc),gas_ib(nc);
    for(std::size_t c=0;c<nc;++c){
        const double x=geo.cell_centres[c].x;
        E[c]={1./3.+.009*x,1./3.-.009*x,1./3.,.014,.006,.007};
        ka[c]=.3+.02*std::abs(x);
        kt[c]=3.+.1*std::abs(x);
        gas_ib[c]=blackbody_intensity(420.+10.*x);
    }
    DomSmmGrayWallV49 lo(mesh,geo,directions,walls,E,kt,ka);
    // Deliberately mismatched physical audit coefficients must fail closed;
    // energy conservation for a DIFFERENT operator cannot certify this solve.
    if(strictly_qualified) {
        bool refused=false;
        auto wrong_ka=ka;
        wrong_ka[0]+=1.e-3;
        try {
            (void)solve_smm_anisotropic_v52(mesh,geo,lo,E,kt,wrong_ka,gas_ib);
        } catch(const std::invalid_argument&) { refused=true; }
        require(refused,"V52: mismatched audit absorption was accepted");
        refused=false;
        auto wrong_E=E;
        wrong_E[0][3]+=1.e-4;
        try {
            (void)solve_smm_anisotropic_v52(mesh,geo,lo,wrong_E,kt,ka,gas_ib);
        } catch(const std::invalid_argument&) { refused=true; }
        require(refused,"V52: mismatched Eddington tensor was accepted");
    }
    const auto result=solve_smm_anisotropic_v52(mesh,geo,lo,E,kt,ka,gas_ib);
    std::cout<<"V52_AMG cells="<<nc<<" depth="<<depth
        <<" linear_gate="<<result.linear_strict_pass
        <<" physical_gate="<<result.physical_pass
        <<" qualified="<<result.qualified
        <<" iters="<<result.iterations
        <<" levels="<<result.amg_levels
        <<" lin_l2="<<result.true_linear_l2
        <<" lin_max="<<result.true_linear_max
        <<" cell="<<result.max_cell_balance
        <<" global="<<result.global_energy_balance
        <<" seconds="<<result.total_seconds
        <<" diagnostic="<<result.diagnostic<<'\n';
    require(result.amg_setup_ok,"V52: AMG setup failed");
    require(result.physical_pass,"V52: physical gate failed at 640 cells");
    if(strictly_qualified){
        require(result.qualified,"V52: strict AMG gate failed");
        require(result.accepted_J.size()==nc,"V52: qualified state not published");
    }else{
        require(!result.linear_strict_pass,"V52: extreme-grid negative regression disappeared");
        require(result.accepted_J.empty(),"V52: unqualified state published");
    }
    // Negative tolerance test, no changes to the historical default gate.
    if(strictly_qualified){
        DomSmmAnisotropicV52Options strict;
        strict.physical_tolerance=1e-30;
        const auto denied=solve_smm_anisotropic_v52(mesh,geo,lo,E,kt,ka,gas_ib,strict);
        require(!denied.qualified && denied.accepted_J.empty(),
                "V52: impossible physical threshold published an iterate");
    }
}
}
int main(int argc,char** argv){
    try{
        const std::string mode=argc==1?"--quick":argc==2?argv[1]:"invalid";
        if(mode=="--quick") check_case(.06,true);
        else if(mode=="--stress") {check_case(.06,true);check_case(.03,true);check_case(.012,true);check_case(.006,false);check_case(.003,false);}
        else throw std::invalid_argument("V52: --quick or --stress");
        std::cout<<"V52_AMG_GATE_PASS\n";
    }catch(const std::exception& e){std::cerr<<"V52_AMG_GATE_FAIL: "<<e.what()<<'\n';return 1;}
}
