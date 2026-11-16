// Non-isothermal *walls* and isothermal emitting/absorbing participating medium.
// Independent 1D pure-absorption RTE exponential solution is the DOM oracle.
// Lateral mesh faces are EMPTY only for DOM, modelling a homogeneous 1D slab;
// the same 3D cell geometry is used for the P1/Marshak approximation.
#include "cfdx/physics/p1_marshak_v39.h"
#include "cfdx/physics/dom_quadrature_v39.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
void req(bool b,const char* m){if(!b)throw std::runtime_error(m);}
constexpr double PI=3.14159265358979323846264338327950288;
const double a=.45,bulk=500.,leftT=850.,rightT=300.;
inline double B(double t){return blackbody_intensity(t);}
inline double exponential_G(double x){
    // Simpson angular integral, 32768 midpoint widths: angular oracle only,
    // independent of the CFDX finite-volume DOM/Marshak implementations.
    const int N=32768;
    long double accum=0.;
    for(int k=0;k<N;++k){
        const double mu=(k+.5)/N;
        const double il=B(bulk)+(B(leftT)-B(bulk))*std::exp(-a*x/mu);
        const double ir=B(bulk)+(B(rightT)-B(bulk))*std::exp(-a*(1.-x)/mu);
        accum+=il+ir;
    }
    return static_cast<double>(2.*PI*accum/N);
}
struct Results{double dom=0,p1=0,wallmoment=0;};
Results run(std::size_t nx,std::size_t polar){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(nx,3,3,0.,0.,.2);
    auto geometry=build_fv_geometry(mesh);
    Field<double,Location::CELL> T(mesh.n_cells(),"T","K",1);
    Field<double,Location::CELL> Gp(mesh.n_cells(),"Gp","W/m2",1);
    Field<double,Location::CELL> Qp(mesh.n_cells(),"Qp","W/m3",1);
    Field<double,Location::CELL> Gd(mesh.n_cells(),"Gd","W/m2",1);
    Field<double,Location::CELL> Qd(mesh.n_cells(),"Qd","W/m3",1);
    T.fill(bulk);Gp.fill(4*blackbody_emissive_power(bulk));
    Gd=Gp;Qp.fill(0);Qd.fill(0);
    P1MarshakWalls pm;
    DomWallBoundaryConditions dw;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
        const auto name=mesh.boundary().patch(p).name;
        if(name=="west") {pm[name]={1.,leftT};dw[name]={1.,leftT};}
        else if(name=="east") {pm[name]={1.,rightT};dw[name]={1.,rightT};}
        else {pm[name]={0.,bulk};dw[name]={1.,bulk};}
    }
    std::vector<double>abs(mesh.n_cells(),a),scat(mesh.n_cells(),0.);
    const auto pr=solve_p1_radiation_marshak(mesh,geometry,T,Gp,Qp,abs,scat,pm);
    req(pr.converged,"P1 1D black wall slab not converged");
    auto dirs=make_dom_tensor_gauss_legendre_v39(polar,8);
    // Tensor GL has polar z. Rotate axes to integrate the wall-normal mu=x
    // as a Gauss-Legendre ordinate for a clean one-dimensional comparison.
    for(auto& d:dirs){std::swap(d.dx,d.dz);}
    const auto moment=dom_hemisphere_moments(dirs,{1.,0.,0.});
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
        const auto n=mesh.boundary().patch(p).name;
        if(n!="west" && n!="east") mesh.boundary().patch(p).type=PatchType::EMPTY;
    }
    RadiationTransportControls c;
    c.absorption=a;c.scattering=0.;
    c.max_iterations=100;c.tolerance=1.e-9;c.linear_tolerance=1.e-12;
    c.intensity_relaxation=1.;
    const auto dr=solve_participating_radiation_diffuse_gray_walls(
        mesh,geometry,T,Gd,Qd,dirs,c,dw);
    req(dr.converged,"DOM 1D black wall slab not converged");
    long double ed=0.,ep=0.,w=0.;
    for(std::size_t i=0;i<mesh.n_cells();++i){
        const double ref=exponential_G(geometry.cell_centres[i].x);
        ed+=std::pow((Gd(i)-ref)/ref,2)*geometry.cell_volumes[i];
        ep+=std::pow((Gp(i)-ref)/ref,2)*geometry.cell_volumes[i];
        w+=geometry.cell_volumes[i];
        req(std::abs(Qd(i)-a*(4*blackbody_emissive_power(bulk)-Gd(i)))<1e-9,
            "DOM source not physically consistent");
    }
    const double errd=std::sqrt(static_cast<double>(ed/w));
    const double errp=std::sqrt(static_cast<double>(ep/w));
    std::cout<<std::scientific<<std::setprecision(11)
      <<"RTE_SLAB nx="<<nx<<" polar="<<polar<<" angles="<<dirs.size()
      <<" dom_L2_relative="<<errd<<" p1_L2_relative="<<errp
      <<" hemi_pi_error="<<moment.relative_pi_error
      <<" dom_iterations="<<dr.iterations<<"\n";
    return {errd,errp,moment.relative_pi_error};
}
}
int main(){
    try{
        const auto a=run(12,8),b=run(24,12),c=run(48,16);
        req(c.dom<b.dom && b.dom<a.dom,"DOM RTE does not contract against exponential oracle");
        req(c.dom<.03,"DOM nonisothermal RTE oracle not accurate enough");
        // P1 is a diffusion approximation, not the transport-equation oracle.
        req(c.p1>c.dom,"P1 unexpectedly beats high-angular DOM in thin/intermediate slab");
        std::cout<<"RTE_SLAB_V40_CERTIFICATE PASS\n";return 0;
    }catch(const std::exception& e){
        std::cerr<<"RTE_SLAB_V40_CERTIFICATE FAIL "<<e.what()<<"\n";return 1;
    }
}
