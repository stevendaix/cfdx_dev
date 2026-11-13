// V56 research HO-SMM/AMG with persistent split state and mandatory independent
// angular auditor. Negative tests fail-closed on a low-component perturbation.
#include "cfdx/physics/dom_hilo_holo_v56.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
static void check(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
static void run(std::size_t nx,std::size_t ny,std::size_t nz,double depth,
                unsigned outer,bool require_certified) {
    auto mesh=cfdx::testing::thermal_v36::make_mesh(nx,ny,nz,.018,0.,depth);
    for(std::size_t i=0;i<mesh.points().size();++i){
        const double x=mesh.points().x_data()[i];
        const double y=mesh.points().y_data()[i];
        const double z=mesh.points().z_data()[i];
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
    std::vector<DomHiLoV55> start(ho.dimension(),DomHiLoV55(blackbody_intensity(420.)));
    for(const auto& w:ho.wall_faces())
        if(w.reflect_index<ho.n_reflecting_faces())start[nc+w.reflect_index]=DomHiLoV55(w.ib);
    auto high_start=std::vector<double>(start.size());
    for(std::size_t i=0;i<start.size();++i)high_start[i]=start[i].hi;
    const auto snapshot=evaluate_dom_variable_v50(ho,high_start);
    std::vector<double> kt(nc);
    for(std::size_t i=0;i<nc;++i)kt[i]=ka[i]+ks[i];
    DomSmmGrayWallV49 lo(mesh,geo,dirs,walls,snapshot.eddington,kt,ka);
    DomStableWlsV54 stable(mesh,geo,lo,snapshot.eddington,kt,ka);
    std::vector<DomHiLoV55> projected(nc);
    for(std::size_t i=0;i<nc;++i)projected[i]=DomHiLoV55(snapshot.projected[i]);
    const auto source=build_dom_hilo_face_closure_v56(mesh,geo,lo,stable,
         projected,snapshot.gas_ib,ka,snapshot.face_current);
    check(source.identity<=1e-13,"V56: exact HO/LO ledger identity failure");
    check(source.corrected_rhs.size()==nc,"V56: split RHS missing");
    const auto reference=audit_dom_angular_v54(ho,high_start,1e-9,2);
    check(std::isfinite(reference.max_ordinate_equation),"V56: angular auditor unavailable");
    auto altered_faces=snapshot.face_current;
    altered_faces[0]+=.125;
    const auto perturbed_source=build_dom_hilo_face_closure_v56(mesh,geo,lo,stable,
         projected,snapshot.gas_ib,ka,altered_faces);
    const auto owner=mesh.ownership().owner(0);
    check(std::abs(dom_to_long_double_v55(perturbed_source.corrected_rhs[owner]-
         source.corrected_rhs[owner]))>.12L,
         "V56: HO face refresh ignored by split source");
    check(perturbed_source.physical_ho>source.physical_ho,
          "V56: altered HO current escaped physical ledger");
    auto poisoned=projected;
    poisoned[0].lo+=.02*poisoned[0].hi;
    const auto low_dirty=stable.apply_hilo_v55(poisoned);
    const auto low_clean=stable.apply_hilo_v55(projected);
    long double difference=0.;
    for(std::size_t c=0;c<nc;++c)
        difference=std::max(difference,std::abs(dom_to_long_double_v55(
             low_dirty[c]-low_clean[c])));
    check(difference>1e-6L,"V56: low-component poisoning invisible to LO");
    if(outer==0){
        std::cout<<"V56_EXTREME_PREFLIGHT cells="<<nc<<" depth="<<depth
                 <<" identity="<<source.identity
                 <<" ho_cell="<<source.physical_ho
                 <<" angular="<<reference.max_ordinate_equation
                 <<" deliberately_no_coupled_solve=1\n";
        return;
    }
    DomV56Options options;
    options.v53.base.base.max_outer=outer;
    options.v53.base.base.max_backtracks=9;
    options.v53.base.base.physical_tolerance=1e-9;
    options.v53.base.base.required_contraction=.999;
    options.v53.base.base.allow_picard_fallback=true;
    options.linear.strict=1e-13;
    options.simulated_partitions=2;
    const auto result=solve_dom_holo_hilo_v56(ho,start,options);
    std::cout<<"V56_HO_SPLIT cells="<<nc<<" depth="<<depth
             <<" certified="<<result.experimental_certified
             <<" production="<<result.production_published
             <<" accepted="<<result.accepted
             <<" lo_accepted="<<result.lo_accepted
             <<" lo_strict="<<result.lo_strict
             <<" tensor_refresh="<<result.max_tensor_refresh
             <<" face_refresh="<<result.max_face_correction_refresh
             <<" refreshed="<<result.refreshed
             <<" krylov="<<result.krylov
             <<" ho="<<result.final_ho
             <<" lo="<<result.final_lo
             <<" original_double_corrected="<<result.final_original_double_corrected
             <<" angular="<<result.final_angular
             <<" identity="<<result.last_identity
             <<" last_linear="<<result.last_linear
             <<" diagnostic="<<result.diagnostic<<'\n';
    check(!result.production_published && result.published_double_state.empty(),
          "V56: split state incorrectly published into production double field");
    check(result.refreshed>=result.accepted+1,
          "V56: coefficients not refreshed on accepted states");
    check(result.last_identity<=1e-9,
          "V56: final HO/LO face identity was not reassembled");
    check(std::isfinite(result.final_original_double_corrected),
          "V56: historical V48 double operator was not reassembled");
    check(result.ho_history.size()==result.lo_history.size(),
          "V56: inconsistent physical audit history");
    if(require_certified){
        check(result.experimental_certified,"V56: independently certified HOLO not converged");
        check(result.lo_accepted>0 && result.lo_strict>0,
              "V56: certified solution skipped compensated AMG");
        check(result.max_tensor_refresh>1e-12 &&
              result.max_face_correction_refresh>1e-12,
              "V56: accepted states did not refresh E or HO/LO faces");
        check(result.final_ho<=1e-9 && result.final_lo<=1e-9 &&
              result.final_angular<=1e-9 && result.last_linear<=1e-13,
              "V56: claimed convergence without all strict gates");
        check(result.experimental_state.size()==ho.dimension(),
              "V56: certified split state missing");
        auto rounded=std::vector<double>(result.experimental_state.size());
        for(std::size_t i=0;i<rounded.size();++i)rounded[i]=result.experimental_state[i].hi;
        const auto independent=audit_dom_angular_v54(ho,rounded,1e-9,7);
        check(independent.accepted,"V56: independent final angular rerun rejected");
    }else check(!result.experimental_certified ||
                (result.final_angular<=1e-9 && result.final_lo<=1e-9),
               "V56: false research certificate");
    // Refusal of altered controls must be independent of the outer budget.
    auto invalid=options;
    invalid.linear.strict=1.1e-13;
    bool refused=false;
    try{(void)solve_dom_holo_hilo_v56(ho,start,invalid);}
    catch(const std::invalid_argument&){refused=true;}
    check(refused,"V56: relaxed 1e-13 precision gate was accepted");
    if(nc==27){
        auto poison_initial=start;
        poison_initial[0].lo+=.02*poison_initial[0].hi;
        bool hidden_rejected=false;
        try{(void)solve_dom_holo_hilo_v56(ho,poison_initial,options);}
        catch(const std::invalid_argument&){hidden_rejected=true;}
        check(hidden_rejected,"V56: HO-invisible 2% low component was accepted");
    }
}
int main(int argc,char** argv){
    try{
        const std::string mode=argc==1?"--quick":argc==2?argv[1]:"invalid";
        if(mode=="--quick")run(3,3,3,.03,70,true);
        else if(mode=="--extreme")run(8,8,10,.003,6,false);
        else if(mode=="--extreme-one")run(8,8,10,.003,1,false);
        else if(mode=="--extreme-preflight")run(8,8,10,.003,0,false);
        else if(mode=="--closure")run(3,3,3,.006,2,false);
        else throw std::invalid_argument("V56: unknown mode");
        std::cout<<"V56_HO_SPLIT_PASS "<<mode<<'\n';return 0;
    }catch(const std::exception& e){
        std::cerr<<"V56_HO_SPLIT_FAIL "<<e.what()<<'\n';return 1;
    }
}
