#include "cfdx/physics/dom_independent_angular_audit_v54.h"
#include "cfdx/physics/dom_holo_certificate_v54.h"
#include "thermal_3d_mesh_v36.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace cfdx::core;
using namespace cfdx::physics;
static void check(bool b,const char* s){if(!b)throw std::runtime_error(s);}
static void check_audit(const DomV54AngularAudit& a,double tol){
    std::cout<<"V54_ANGULAR ordinate="<<a.max_ordinate_equation
        <<" wall_R="<<a.max_wall_reflection<<" cell="<<a.max_cell_energy
        <<" global="<<a.global_energy<<" fixed="<<a.max_fixed_moment
        <<" partition="<<a.face_partition_error<<" accepted="<<a.accepted<<'\n';
    check(a.accepted && a.max_ordinate_equation<=tol &&
          a.max_wall_reflection<=tol && a.face_partition_error<=tol,
          "V54 independent angular equations not accepted");
}
static void run(bool stress){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(stress?8:3,stress?8:3,
                                                    stress?10:3,.018,0.,stress?.06:.03);
    for(std::size_t i=0;i<mesh.points().size();++i){
        const double x=mesh.points().x_data()[i],y=mesh.points().y_data()[i];
        const double z=mesh.points().z_data()[i];
        const double xr=std::cos(.41)*x-std::sin(.41)*y;
        const double yr=std::sin(.41)*x+std::cos(.41)*y;
        mesh.points().set(i,std::cos(.33)*xr+std::sin(.33)*z,yr,
                            -std::sin(.33)*xr+std::cos(.33)*z);
    }
    const auto geo=build_fv_geometry(mesh);
    const auto dirs=make_dom_tensor_gauss_legendre_v39(4,16);
    const std::size_t nc=mesh.n_cells();
    DomMomentSchurSystemV44::CellField T(nc,"temperature","K",1);
    std::vector<double> ka(nc),ks(nc);
    for(std::size_t c=0;c<nc;++c){
        T(c)=420.+5.*geo.cell_centres[c].x;
        ka[c]=1.2+.2*std::abs(geo.cell_centres[c].x);
        ks[c]=.35+.1*std::abs(geo.cell_centres[c].z);
    }
    DomWallBoundaryConditions bc;
    for(std::size_t i=0;i<mesh.boundary().n_patches();++i){
        const auto name=mesh.boundary().patch(i).name;
        bc[name]={name=="bottom"?0.:name=="west"?.45:.85,
                  420.+(name=="east"?20.:name=="west"?-12.:0.)};
    }
    RadiationTransportControls transport;
    transport.absorption=1.2;transport.scattering=.35;
    transport.tolerance=1e-9;transport.linear_tolerance=1e-14;
    transport.linear_max_iterations=2000;
    DomMomentSchurSystemV44 ho(mesh,geo,T,dirs,bc,transport,ka,ks);
    std::vector<double> initial(ho.dimension(),blackbody_intensity(420.));
    for(const auto& w:ho.wall_faces())if(w.reflect_index<ho.n_reflecting_faces())
        initial[nc+w.reflect_index]=w.ib;
    DomV53Options opts;
    opts.base.base.max_outer=70;opts.base.base.max_backtracks=9;
    opts.base.base.physical_tolerance=1e-9;
    const auto converged=solve_dom_holo_amg_v53(ho,initial,opts);
    check(converged.converged && converged.coupled_amg,"V54 V53 base not converged");
    const auto certified=certify_dom_holo_v54(ho,converged,1e-9,7);
    check(certified.published && certified.accepted_state==converged.state,
          "V54 independent certificate did not publish verified state");
    for(const auto p:{1u,2u,7u,32u}){
        const auto a=audit_dom_angular_v54(ho,converged.state,1e-9,p);
        check_audit(a,1e-9);
    }
    const auto snapshot=evaluate_dom_variable_v50(ho,converged.state);
    std::vector<double> kt(nc),J(nc);
    for(std::size_t c=0;c<nc;++c){kt[c]=ka[c]+ks[c];J[c]=snapshot.projected[c];}
    DomSmmGrayWallV49 lo(mesh,geo,dirs,bc,snapshot.eddington,kt,ka);
    const auto identity=lo.consistent_source(J,snapshot.gas_ib,snapshot.face_current);
    check(identity.max_discrete_identity_relative<1e-10,"V54 HO-LO identity failed");
    auto corrupt=converged.state;corrupt[0]*=1.02;
    auto corrupted_converged=converged;
    corrupted_converged.state=corrupt;
    const auto denied=certify_dom_holo_v54(ho,std::move(corrupted_converged),1e-9,3);
    check(!denied.published && !denied.pilot.converged &&
          denied.pilot.state.empty() && denied.accepted_state.empty(),
          "V54 invalid state leaked through final independent certificate");
    const auto negative=audit_dom_angular_v54(ho,corrupt,1e-9,3);
    check(!negative.accepted && negative.max_fixed_moment>1e-9,
          "V54 false physical acceptance after moment corruption");
    const auto perturbed=evaluate_dom_variable_v50(ho,corrupt);
    for(std::size_t c=0;c<nc;++c)J[c]=perturbed.projected[c];
    DomSmmGrayWallV49 disturbed_lo(mesh,geo,dirs,bc,perturbed.eddington,kt,ka);
    const auto disturbed_identity=disturbed_lo.consistent_source(
        J,perturbed.gas_ib,perturbed.face_current);
    check(disturbed_identity.max_discrete_identity_relative<1e-10,
          "V54 correction identity lost on HO-perturbed state");
    std::cout<<"V54_IDENTITY_VS_ORACLE cells="<<nc
        <<" identity_bad_state="<<disturbed_identity.max_discrete_identity_relative
        <<" independent_fixed_bad_state="<<negative.max_fixed_moment
        <<" independent_rejected="<<!negative.accepted<<"\n";
}
static void isothermal_analytic(bool stress){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(stress?8:3,stress?8:3,
                                                    stress?10:3,.018,0.,stress?.06:.03);
    for(std::size_t i=0;i<mesh.points().size();++i){
        const double x=mesh.points().x_data()[i],y=mesh.points().y_data()[i];
        const double z=mesh.points().z_data()[i];
        const double xr=std::cos(.41)*x-std::sin(.41)*y;
        const double yr=std::sin(.41)*x+std::cos(.41)*y;
        mesh.points().set(i,std::cos(.33)*xr+std::sin(.33)*z,yr,
                            -std::sin(.33)*xr+std::cos(.33)*z);
    }
    const auto geo=build_fv_geometry(mesh);
    const auto dirs=make_dom_tensor_gauss_legendre_v39(4,16);
    const std::size_t nc=mesh.n_cells();
    DomMomentSchurSystemV44::CellField T(nc,"T","K",1);T.fill(420.);
    std::vector<double> ka(nc),ks(nc);
    for(std::size_t c=0;c<nc;++c){
        ka[c]=.15+2.*std::abs(geo.cell_centres[c].x);
        ks[c]=2.+3.*std::abs(geo.cell_centres[c].z);
    }
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
        const auto name=mesh.boundary().patch(p).name;
        walls[name]={name=="bottom"?0.:name=="west"?.45:.85,420.};
    }
    RadiationTransportControls ctrl;
    ctrl.absorption=.15;ctrl.scattering=2.;
    ctrl.tolerance=1e-9;ctrl.linear_tolerance=1e-14;
    ctrl.linear_max_iterations=2000;
    DomMomentSchurSystemV44 ho(mesh,geo,T,dirs,walls,ctrl,ka,ks);
    const double black=blackbody_intensity(420.);
    std::vector<double> known(ho.dimension(),black);
    const auto audit=audit_dom_angular_v54(ho,known,1e-9,7);
    check_audit(audit,1e-9);
    check(audit.maximum_wall_abs_watts<1e-6,
          "V54 isothermal exact reference leaked gray-wall energy");
    std::cout<<"V54_ANALYTIC_ISOTHERMAL cells="<<nc
             <<" blackbody="<<black
             <<" max_wall_watts="<<audit.maximum_wall_abs_watts
             <<" accepted="<<audit.accepted<<'\n';
}
int main(int argc,char** argv){
    try{const bool stress=argc==2 && std::string(argv[1])=="--stress-640";
        const bool analytic=argc==2 && std::string(argv[1])=="--analytic";
        if(argc>2 || (argc==2 && !stress && !analytic))
            throw std::invalid_argument("bad args");
        if(analytic){isothermal_analytic(false);isothermal_analytic(true);}
        else run(stress);
        std::cout<<"V54_INDEPENDENT_PASS\n";return 0;
    }catch(const std::exception& e){std::cerr<<"V54_INDEPENDENT_FAIL "<<e.what()<<'\n';return 1;}
}
