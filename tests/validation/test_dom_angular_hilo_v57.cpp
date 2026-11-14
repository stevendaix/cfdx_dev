#include "cfdx/physics/dom_angular_hilo_v57.h"
#include "thermal_3d_mesh_v36.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
static void require(bool condition,const char* description){
    if(!condition)throw std::runtime_error(description);
}
static void one(std::size_t nx,std::size_t ny,std::size_t nz,double depth,bool analytic=false) {
    auto mesh=cfdx::testing::thermal_v36::make_mesh(nx,ny,nz,.018,0.,depth);
    for(std::size_t i=0;i<mesh.points().size();++i) {
        const double x=mesh.points().x_data()[i];
        const double y=mesh.points().y_data()[i];
        const double z=mesh.points().z_data()[i];
        const double xr=std::cos(.41)*x-std::sin(.41)*y;
        const double yr=std::sin(.41)*x+std::cos(.41)*y;
        mesh.points().set(i,std::cos(.33)*xr+std::sin(.33)*z,
                          yr,-std::sin(.33)*xr+std::cos(.33)*z);
    }
    const auto geo=build_fv_geometry(mesh);
    const auto dirs=make_dom_tensor_gauss_legendre_v39(4,16);
    const auto nc=mesh.n_cells();
    DomMomentSchurSystemV44::CellField T(nc,"T","K",1);
    std::vector<double> ka(nc),ks(nc);
    for(std::size_t c=0;c<nc;++c){
        T(c)=analytic?420.:420.+5.*geo.cell_centres[c].x;
        ka[c]=1.2+.2*std::abs(geo.cell_centres[c].x);
        ks[c]=.35+.1*std::abs(geo.cell_centres[c].z);
    }
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
        const auto name=mesh.boundary().patch(p).name;
        walls[name]={name=="bottom"?0.:name=="west"?.45:.85,
                     analytic?420.:420.+(name=="east"?20.:name=="west"?-12.:0.)};
    }
    RadiationTransportControls ctrl;
    ctrl.absorption=1.2;ctrl.scattering=.35;
    ctrl.tolerance=1e-9;ctrl.linear_tolerance=1e-14;
    ctrl.linear_max_iterations=2000;
    DomMomentSchurSystemV44 ho(mesh,geo,T,dirs,walls,ctrl,ka,ks);
    DomAngularHiloTransportV57 hilo(ho);
    std::vector<DomHiLoV55> x(ho.dimension(),DomHiLoV55(blackbody_intensity(420.)));
    for(const auto& w:ho.wall_faces())if(w.reflect_index<ho.n_reflecting_faces())
        x[nc+w.reflect_index]=DomHiLoV55(w.ib);
    const auto solution=hilo.solve(x);
    const auto audit=audit_dom_angular_hilo_v57(ho,x,solution,1e-9);
    require(solution.max_true_csr_relative<=1e-13,"V57 angular CSR true residual");
    require(audit.max_ordinate_equation<=1e-9,"V57 angular independent face residual");
    if(analytic){
        require(audit.accepted,"V57 isothermal exact equilibrium failed");
        require(audit.max_fixed_moment<=1e-11 && audit.global_energy<=1e-11,
                "V57 isothermal analytic oracle inaccurate");
    }
    const auto high_only=ho.project([&]{
        std::vector<double> hi(x.size());
        for(std::size_t i=0;i<hi.size();++i)hi[i]=x[i].hi;
        return hi;
    }(),true);
    double same=0.;
    for(std::size_t i=0;i<high_only.size();++i)
        same=std::max(same,std::abs(high_only[i]-solution.moment[i].hi));
    require(same<1e-8,"V57 split angular departs from V44 on simple input");
    auto low_poison=x;
    low_poison[0].lo+=.02*x[0].hi;
    const auto perturbed=hilo.solve(low_poison);
    const long double delta=std::abs(dom_to_long_double_v55(
        perturbed.moment[0]-solution.moment[0]));
    require(delta>1e-8L,"V57 angular transport ignores low component");
    const auto bad=audit_dom_angular_hilo_v57(ho,low_poison,perturbed,1e-9);
    require(!bad.accepted,"V57 bad HO fixed point was accepted");
    auto faulty=solution;
    faulty.intensities[0][0].lo+=.02*faulty.intensities[0][0].hi;
    const auto wrong_intensity=audit_dom_angular_hilo_v57(ho,x,faulty,1e-9);
    require(!wrong_intensity.accepted &&
            wrong_intensity.max_ordinate_equation>1e-9,
            "V57 independent FV audit failed to detect faulty intensity");
    faulty=solution;
    faulty.face_current[0].lo+=.125;
    const auto wrong_face=audit_dom_angular_hilo_v57(ho,x,faulty,1e-9);
    require(!wrong_face.accepted && wrong_face.max_face_recompute>1e-9,
            "V57 independent face audit failed to detect stale current");
    std::cout<<"V57_ANGULAR cells="<<nc<<" depth="<<depth
             <<" csr="<<solution.max_true_csr_relative
             <<" face="<<audit.max_ordinate_equation
             <<" fixed="<<audit.max_fixed_moment
             <<" low_delta="<<static_cast<double>(delta)
             <<" baseline="<<same<<" analytic="<<analytic<<" passed=1\n";
}
int main(int argc,char** argv){
    try{
        const std::string mode=argc==2?argv[1]:"--quick";
        if(mode=="--quick")one(3,3,3,.03);
        else if(mode=="--extreme")one(8,8,10,.003);
        else if(mode=="--analytic") {
            one(3,3,3,.03,true);one(8,8,10,.003,true);
        }
        else if(mode=="--family")
            for(const double depth:{.06,.03,.012,.006,.003})one(8,8,10,depth);
        else throw std::invalid_argument("unknown V57 mode");
        return 0;
    }catch(const std::exception& e){
        std::cerr<<"V57_ANGULAR_FAIL "<<e.what()<<'\n';return 1;
    }
}
