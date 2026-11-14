// V52: verify fresh HO/LO face corrections at 640 cells on the ACTUAL angular
// Schur transport, including reflecting gray walls and variable optics.
// Deliberately does not claim a coupled SMM/AMG acceleration at 640 cells.
#include "cfdx/physics/dom_holo_global_v51.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
void require(bool pass,const char* message){if(!pass)throw std::runtime_error(message);}
}
int main(){try{
    auto mesh=cfdx::testing::thermal_v36::make_mesh(8,8,10,.018,0.,.06);
    for(std::size_t i=0;i<mesh.points().size();++i){
        const double x=mesh.points().x_data()[i],y=mesh.points().y_data()[i];
        const double z=mesh.points().z_data()[i];
        const double xr=std::cos(.41)*x-std::sin(.41)*y;
        const double yr=std::sin(.41)*x+std::cos(.41)*y;
        mesh.points().set(i,std::cos(.33)*xr+std::sin(.33)*z,
                            yr,-std::sin(.33)*xr+std::cos(.33)*z);
    }
    const auto geo=build_fv_geometry(mesh);
    const auto directions=make_dom_tensor_gauss_legendre_v39(4,16);
    const auto nc=mesh.n_cells();
    DomMomentSchurSystemV44::CellField T(nc,"T","K",1);
    std::vector<double> ka(nc),ks(nc);
    for(std::size_t c=0;c<nc;++c){
        T(c)=420.+5.*geo.cell_centres[c].x;
        ka[c]=1.2+.2*std::abs(geo.cell_centres[c].x);
        ks[c]=.35+.1*std::abs(geo.cell_centres[c].z);
    }
    RadiationTransportControls control;
    control.absorption=1.2;control.scattering=.35;
    control.tolerance=1e-9;control.linear_tolerance=1e-14;
    control.linear_max_iterations=2000;
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
        const auto name=mesh.boundary().patch(p).name;
        walls[name]={name=="bottom"?0.:name=="west"?.45:.85,
                     420.+(name=="east"?20.:name=="west"?-12.:0.)};
    }
    DomMomentSchurSystemV44 ho(mesh,geo,T,directions,walls,control,ka,ks);
    std::vector<double> x(ho.dimension(),blackbody_intensity(420.));
    for(const auto& w:ho.wall_faces())
        if(w.reflect_index<ho.n_reflecting_faces())
            x[nc+w.reflect_index]=w.ib;
    const auto snapshot_and_closure=[&](const std::vector<double>& trial){
        auto snap=evaluate_dom_variable_v50(ho,trial);
        std::vector<double> kt(nc);
        for(std::size_t c=0;c<nc;++c)kt[c]=ka[c]+ks[c];
        DomSmmGrayWallV49 lo(mesh,geo,directions,walls,snap.eddington,kt,ka);
        std::vector<double> projected(snap.projected.begin(),snap.projected.begin()+static_cast<std::ptrdiff_t>(nc));
        return std::make_pair(snap.true_relative,
                   lo.consistent_source(projected,snap.gas_ib,snap.face_current));
    };
    const auto before=snapshot_and_closure(x);
    DomV51Options options;
    options.base.max_outer=70;
    options.base.physical_tolerance=1e-9;
    options.base.linear_tolerance=1e-12;
    options.use_low_order=false; // SMM/AMG at 640 still independent gate.
    options.anderson_depth=6;
    auto result=solve_dom_holo_global_v51(ho,x,options);
    require(result.converged,"V52 HO: 640-cell full angular transport not converged");
    require(result.final_true<=1e-9,"V52 HO: physical residual rejected");
    const auto after=snapshot_and_closure(result.final_moments);
    require(before.second.max_discrete_identity_relative<1e-10,
        "V52 HO: initial face identity not conservative");
    require(after.second.max_discrete_identity_relative<1e-10,
        "V52 HO: accepted face identity not conservative");
    double delta=0.;
    const auto& a=before.second.face_ho_minus_lo_watts;
    const auto& b=after.second.face_ho_minus_lo_watts;
    require(a.size()==b.size(),"V52 HO: face set changed");
    for(std::size_t f=0;f<a.size();++f)
        delta=std::max(delta,std::abs(a[f]-b[f]));
    require(delta>1e-9,"V52 HO: face correction never refreshed");
    require(after.second.max_ho_energy_residual_relative<=1e-9,
        "V52 HO: accepted cell energy defect not closed");
    std::cout<<"V52_HO_REFRESH cells="<<nc
             <<" states="<<result.accepted_steps
             <<" final="<<result.final_true
             <<" face_change="<<delta
             <<" closure_before="<<before.second.max_discrete_identity_relative
             <<" closure_after="<<after.second.max_discrete_identity_relative
             <<" ho_energy_after="<<after.second.max_ho_energy_residual_relative
             <<" smm_amg_coupled=NO\n";
    std::cout<<"V52_HO_REFRESH_PASS\n";
}catch(const std::exception& e){
    std::cerr<<"V52_HO_REFRESH_FAIL: "<<e.what()<<'\n';return 1;
}}
