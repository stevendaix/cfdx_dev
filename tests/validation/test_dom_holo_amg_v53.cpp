#include "cfdx/physics/dom_holo_amg_global_v53.h"
#include "thermal_3d_mesh_v36.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
static void check(bool x,const char* msg){if(!x)throw std::runtime_error(msg);}
static void run(std::size_t nx,std::size_t ny,std::size_t nz,double depth){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(nx,ny,nz,.018,0.,depth);
    for(std::size_t i=0;i<mesh.points().size();++i){
        const double x=mesh.points().x_data()[i],y=mesh.points().y_data()[i],z=mesh.points().z_data()[i];
        const double xr=std::cos(.41)*x-std::sin(.41)*y;
        const double yr=std::sin(.41)*x+std::cos(.41)*y;
        mesh.points().set(i,std::cos(.33)*xr+std::sin(.33)*z,yr,
                          -std::sin(.33)*xr+std::cos(.33)*z);
    }
    const auto geo=build_fv_geometry(mesh);
    const auto dirs=make_dom_tensor_gauss_legendre_v39(4,16);
    const auto nc=mesh.n_cells();
    DomMomentSchurSystemV44::CellField T(nc,"T","K",1);
    std::vector<double> ka(nc),ks(nc);
    for(std::size_t i=0;i<nc;++i){
        T(i)=420.+5.*geo.cell_centres[i].x;
        ka[i]=1.2+.2*std::abs(geo.cell_centres[i].x);
        ks[i]=.35+.1*std::abs(geo.cell_centres[i].z);
    }
    DomWallBoundaryConditions walls;
    for(std::size_t i=0;i<mesh.boundary().n_patches();++i){
        const auto name=mesh.boundary().patch(i).name;
        walls[name]={name=="bottom"?0.:name=="west"?.45:.85,
                     420.+(name=="east"?20.:name=="west"?-12.:0.)};
    }
    RadiationTransportControls ctrl;
    ctrl.absorption=1.2;ctrl.scattering=.35;
    ctrl.tolerance=1e-9;ctrl.linear_tolerance=1e-14;
    ctrl.linear_max_iterations=2000;
    DomMomentSchurSystemV44 ho(mesh,geo,T,dirs,walls,ctrl,ka,ks);
    std::vector<double> start(ho.dimension(),blackbody_intensity(420.));
    for(const auto& w:ho.wall_faces())
        if(w.reflect_index<ho.n_reflecting_faces())start[nc+w.reflect_index]=w.ib;
    DomV53Options options;
    options.base.base.max_outer=70;
    options.base.base.max_backtracks=9;
    options.base.base.physical_tolerance=1e-9;
    options.base.anderson_depth=6;
    const auto solved=solve_dom_holo_amg_v53(ho,start,options);
    std::cout<<"V53 cells="<<nc<<" convergence="<<solved.converged
             <<" AMG="<<solved.amg_strict<<" LO_steps="<<solved.lo_accepted
             <<" Anderson="<<solved.anderson_accepted<<" Picard="<<solved.picard_accepted
             <<" HO="<<solved.final_ho<<" LO="<<solved.final_lo
             <<" global="<<solved.final_global<<" last_linear="<<solved.last_linear
             <<" refreshes="<<solved.refreshes<<std::endl;
    check(solved.converged && solved.coupled_amg,"coupled HO+AMG did not converge");
    check(solved.lo_accepted>0 && solved.lo_accepted==solved.accepted,
          "not all steps used SMM/AMG");
    check(solved.final_ho<=1e-9 && solved.final_lo<=1e-9 &&
          solved.final_global<=1e-9,"not both physical gates closed");
    check(solved.last_linear<=1e-13,"strict AMG linear residual not closed");
    check(solved.state.size()==ho.dimension(),"unqualified state published");
    check(solved.refreshes>=solved.accepted+1,"correction not rebuilt");
    const auto audit=evaluate_dom_variable_v50(ho,solved.state);
    check(audit.true_relative<=1e-9,"independent HO check failed");
    // Verify the face correction is actually recomputed with the new HO E.
    const auto initial_snap=evaluate_dom_variable_v50(ho,start);
    const auto closure=[&](const DomV50Snapshot& s){
        std::vector<double> kt(nc),J(nc);
        for(std::size_t i=0;i<nc;++i){
            kt[i]=ka[i]+ks[i];J[i]=s.projected[i];
        }
        DomSmmGrayWallV49 lo(mesh,geo,dirs,walls,s.eddington,kt,ka);
        return lo.consistent_source(J,s.gas_ib,s.face_current);
    };
    const auto before=closure(initial_snap);
    const auto after=closure(audit);
    double change=0.;
    for(std::size_t i=0;i<mesh.n_faces();++i)
        change=std::max(change,std::abs(before.face_ho_minus_lo_watts[i]-
                                        after.face_ho_minus_lo_watts[i]));
    check(change>1e-9 && before.max_discrete_identity_relative<=1e-10 &&
          after.max_discrete_identity_relative<=1e-10,
          "face closure not refreshed between accepted states");
    std::cout<<"V53_FACE_REFRESH change_W="<<change
             <<" final_identity="<<after.max_discrete_identity_relative<<'\n';
    if(nc==27){
        auto impossible=options;
        impossible.base.base.physical_tolerance=1e-30;
        impossible.base.base.max_outer=2;
        const auto bad=solve_dom_holo_amg_v53(ho,start,impossible);
        check(!bad.converged && bad.state.empty(),"impossible physical gate published");
        auto impossible_linear=options;
        impossible_linear.strict_linear=1e-30;
        impossible_linear.max_krylov=2;
        impossible_linear.wide_budget=2;
        impossible_linear.iterative_refinements=0;
        impossible_linear.base.use_anderson=false;
        impossible_linear.base.base.allow_picard_fallback=false;
        impossible_linear.base.base.max_outer=2;
        const auto no_amg=solve_dom_holo_amg_v53(ho,start,impossible_linear);
        check(!no_amg.converged && !no_amg.amg_strict && no_amg.state.empty(),
              "nonconverged AMG allowed final publish");
    }
}
int main(int argc,char** argv){
    try{
        const std::string mode=argc==1?"--quick":argc==2?argv[1]:"invalid";
        if(mode=="--quick")run(3,3,3,.03);
        else if(mode=="--stress-640")run(8,8,10,.06);
        else throw std::invalid_argument("V53: invalid test mode");
        std::cout<<"V53_PASS\n";
    }catch(const std::exception& e){
        std::cerr<<"V53_FAIL "<<e.what()<<'\n';return 1;
    }
}
