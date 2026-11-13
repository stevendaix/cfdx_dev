// V57 HO angular, face current, E tensor and SMM corrector in two-double
// arithmetic; selected acceptance always requires independent angular audit.
#include "cfdx/physics/dom_holo_hilo_v57.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
static void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
static void run(std::size_t nx,std::size_t ny,std::size_t nz,
                double depth,std::size_t outer,bool should_converge){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(nx,ny,nz,.018,0.,depth);
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
    const auto nc=mesh.n_cells();
    DomMomentSchurSystemV44::CellField T(nc,"T","K",1);
    std::vector<double> ka(nc),ks(nc);
    for(std::size_t c=0;c<nc;++c){
        T(c)=420.+5.*geo.cell_centres[c].x;
        ka[c]=1.2+.2*std::abs(geo.cell_centres[c].x);
        ks[c]=.35+.1*std::abs(geo.cell_centres[c].z);
    }
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
        const auto name=mesh.boundary().patch(p).name;
        walls[name]={name=="bottom"?0.:name=="west"?.45:.85,
                     420.+(name=="east"?20.:name=="west"?-12.:0.)};
    }
    RadiationTransportControls ctrl;
    ctrl.absorption=1.2;ctrl.scattering=.35;
    ctrl.tolerance=1e-9;ctrl.linear_tolerance=1e-14;
    ctrl.linear_max_iterations=2000;
    DomMomentSchurSystemV44 ho(mesh,geo,T,dirs,walls,ctrl,ka,ks);
    std::vector<DomHiLoV55> start(ho.dimension(),
        DomHiLoV55(blackbody_intensity(420.)));
    for(const auto& w:ho.wall_faces())
        if(w.reflect_index<ho.n_reflecting_faces())
            start[nc+w.reflect_index]=DomHiLoV55(w.ib);
    DomV57Options opt;
    opt.v56.v53.base.base.max_outer=outer;
    opt.v56.v53.base.base.max_backtracks=9;
    opt.v56.v53.base.base.physical_tolerance=1e-9;
    opt.v56.v53.base.base.required_contraction=.999;
    opt.v56.v53.base.base.allow_picard_fallback=true;
    opt.v56.linear.strict=1e-13;
    opt.v56.simulated_partitions=2;
    if(nc==640)opt.progress=[depth](std::size_t i,double ho,double lo) {
        std::cout<<"V57_PROGRESS depth="<<depth<<" accepted="<<i
                 <<" HO="<<ho<<" LO="<<lo<<std::endl;
    };
    const auto result=solve_dom_holo_hilo_v57(ho,start,opt);
    std::cout<<"V57_HO_SPLIT cells="<<nc<<" depth="<<depth
        <<" certified="<<result.experimental_certified
        <<" published="<<result.production_published
        <<" accepted="<<result.accepted
        <<" rejected="<<result.rejected
        <<" lo_accepted="<<result.lo_accepted
        <<" lo_strict="<<result.lo_strict
        <<" angular_sweeps="<<result.angular_sweeps
        <<" E_refresh="<<result.max_tensor_refresh
        <<" face_refresh="<<result.max_face_refresh
        <<" krylov="<<result.krylov
        <<" ho="<<result.final_ho
        <<" lo="<<result.final_lo
        <<" angular="<<result.final_angular
        <<" energy="<<result.final_energy
        <<" identity="<<result.final_identity
        <<" last_linear="<<result.last_linear
        <<" original_double="<<result.original_double_corrected
        <<" angular_seconds="<<result.angular_seconds
        <<" refresh_seconds="<<result.refresh_seconds
        <<" amg_seconds="<<result.amg_seconds
        <<" reason="<<result.diagnostic<<std::endl;
    require(!result.production_published && result.production_double_state.empty(),
            "V57: double production state was published");
    require(result.refreshed>=1 && result.angular_sweeps>=result.refreshed,
            "V57: angular intensities were not refreshed");
    require(result.final_identity<=1e-13,
            "V57: HO/LO shared-face ledger inconsistent");
    require(result.max_angular_csr<=1e-13,
            "V57: split angular CSR residual exceeded strict 1e-13");
    require(std::isfinite(result.original_double_corrected),
            "V57: original V48 double operator not evaluated");
    require(result.ho_history.size()==result.lo_history.size() &&
            result.ho_history.size()==result.energy_history.size(),
            "V57: inconsistent non-linear history");
    if(should_converge){
        require(result.experimental_certified,"V57: global HO+LO not converged");
        require(result.lo_accepted>0 && result.lo_strict>0,
                "V57: certification bypassed AMG");
        require(result.final_ho<=1e-9 && result.final_lo<=1e-9 &&
                result.final_angular<=1e-9 && result.final_energy<=1e-9 &&
                result.last_linear<=1e-13,
                "V57: false coupled convergence certificate");
        require(result.certified_split_state.size()==ho.dimension(),
                "V57: split state not retained");
        require(result.max_tensor_refresh>1e-12 && result.max_face_refresh>1e-12,
                "V57: E / face correction not refreshed");
        DomAngularHiloTransportV57 independent_transport(ho);
        const auto reread=independent_transport.solve(result.certified_split_state);
        const auto audit=audit_dom_angular_hilo_v57(ho,
            result.certified_split_state,reread,1e-9);
        require(audit.accepted,"V57: final standalone angular audit rejected");
    }else{
        require(!result.experimental_certified ||
                (result.final_ho<=1e-9 && result.final_lo<=1e-9),
                "V57: false publication on nonconverged pilot");
    }
    if(nc==27){
        auto invalid=opt;invalid.v56.linear.strict=1.1e-13;
        bool blocked=false;
        try{(void)solve_dom_holo_hilo_v57(ho,start,invalid);}
        catch(const std::invalid_argument&){blocked=true;}
        require(blocked,"V57: loosened linear threshold accepted");
        invalid=opt;
        invalid.v56.v53.base.base.physical_tolerance=1.1e-9;
        blocked=false;
        try{(void)solve_dom_holo_hilo_v57(ho,start,invalid);}
        catch(const std::invalid_argument&){blocked=true;}
        require(blocked,"V57: loosened physical threshold accepted");
        invalid=opt;invalid.v56.v53.base.base.max_outer=1;
        const auto truncated=solve_dom_holo_hilo_v57(ho,start,invalid);
        require(!truncated.experimental_certified &&
                truncated.certified_split_state.empty(),
                "V57: premature global certificate");
    }
}
int main(int argc,char** argv){
    try{
        const std::string mode=argc==2?argv[1]:"--quick";
        if(mode=="--quick")run(3,3,3,.03,70,true);
        else if(mode=="--extreme-one")run(8,8,10,.003,1,false);
        else if(mode=="--extreme-two")run(8,8,10,.003,2,false);
        else if(mode=="--extreme-three")run(8,8,10,.003,3,false);
        else if(mode=="--extreme")run(8,8,10,.003,70,true);
        else if(mode=="--family-one")
            for(const double d:{.06,.03,.012,.006,.003})run(8,8,10,d,1,false);
        else if(mode=="--family")
            for(const double d:{.06,.03,.012,.006,.003})run(8,8,10,d,70,true);
        else throw std::invalid_argument("V57: unknown option");
        std::cout<<"V57_HO_SPLIT_PASS "<<mode<<std::endl;
        return 0;
    }catch(const std::exception& e){
        std::cerr<<"V57_HO_SPLIT_FAIL "<<e.what()<<std::endl;
        return 1;
    }
}
