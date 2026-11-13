// V52: independent two-layer absorbing-slab characteristic oracle against
// ACTUAL CFDX V50 angular FV matrices, with per-cell discontinuous absorption.
// The 3D lateral dimensions are enlarged and central x-normal faces sampled;
// transverse black walls remain finite and their influence is audited, not
// silently equated to an exact infinite 1D slab.
#include "cfdx/physics/dom_holo_variable_v50.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
constexpr double pi=3.14159265358979323846;
void req(bool ok,const char* what){if(!ok)throw std::runtime_error(what);}
double exact_discrete_flux(const std::vector<DiscreteDirection>& dirs,
                           double k1,double k2,double B1,double B2,
                           double left,double right){
    long double flux=0.;
    for(const auto& d:dirs) {
        const double mu=d.dx;
        if(mu>0.) {
            const double I=B1+(left-B1)*std::exp(-k1*.45/mu);
            flux+=static_cast<long double>(d.weight)*mu*I;
        } else if(mu<0.) {
            const double m=-mu;
            const double I=B2+(right-B2)*std::exp(-k2*.55/m);
            flux-=static_cast<long double>(d.weight)*m*I;
        }
    }
    return static_cast<double>(flux);
}
double exact_continuous_flux(double k1,double k2,double B1,double B2,
                             double left,double right){
    constexpr unsigned n=32768;
    long double total=0.;
    for(unsigned a=0;a<n;++a) {
        const double mu=(a+.5)/double(n);
        const double il=B1+(left-B1)*std::exp(-k1*.45/mu);
        const double ir=B2+(right-B2)*std::exp(-k2*.55/mu);
        total+=mu*(il-ir);
    }
    return 2.*pi*static_cast<double>(total)/n;
}
struct Observed {double spatial=0.,total=0.,angle=0.,ho=0.;};
Observed run(unsigned nx,unsigned polar,double width,double k1,double k2){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(nx,3,3,0.,0.,1.);
    for(std::size_t p=0;p<mesh.points().size();++p)
        mesh.points().set(p,mesh.points().x_data()[p],
                           mesh.points().y_data()[p]*width,
                           mesh.points().z_data()[p]*width);
    const auto geo=build_fv_geometry(mesh);
    auto dirs=make_dom_tensor_gauss_legendre_v39(polar,8);
    // Polar GL abscissas must coincide with the slab's x-normal cosine.
    for(auto& d:dirs)std::swap(d.dx,d.dz);
    const std::size_t nc=mesh.n_cells();
    DomMomentSchurSystemV44::CellField T(nc,"T","K",1);
    std::vector<double> ka(nc),ks(nc,0.);
    constexpr double T1=460.,T2=600.,Tleft=860.,Tright=290.;
    const double B1=blackbody_intensity(T1),B2=blackbody_intensity(T2);
    const double L=blackbody_intensity(Tleft),R=blackbody_intensity(Tright);
    for(std::size_t c=0;c<nc;++c) {
        const bool first=geo.cell_centres[c].x<.45;
        T(c)=first?T1:T2; ka[c]=first?k1:k2;
    }
    DomWallBoundaryConditions bcs;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& name=mesh.boundary().patch(p).name;
        bcs[name]={1.,name=="west"?Tleft:name=="east"?Tright:T1};
    }
    RadiationTransportControls ctrl;
    ctrl.absorption=k1;ctrl.scattering=0.;ctrl.tolerance=1.e-10;
    ctrl.linear_tolerance=1.e-12;ctrl.linear_max_iterations=5000;
    DomMomentSchurSystemV44 ho(mesh,geo,T,dirs,bcs,ctrl,ka,ks);
    std::vector<DomMomentSchurSystemV44::CellField> I;
    auto state=ho.project(std::vector<double>(ho.dimension(),0.),true,&I);
    const auto audit=evaluate_dom_variable_v50(ho,state);
    req(I.size()==dirs.size(),"V52 layered: missing angular intensities");
    req(audit.true_relative<1.e-8,"V52 layered: real 3D FV transport not closed");
    bool found=false;
    double numerical=0.;
    for(std::size_t f=0;f<mesh.n_faces();++f) {
        const auto neighbour=mesh.ownership().neighbour(f);
        if(neighbour<0)continue;
        const auto fc=geo.face_centres[f];
        if(std::abs(fc.x-.45)>1.e-9 || std::abs(fc.y-width/2.)>1.e-9 ||
           std::abs(fc.z-width/2.)>1.e-9)continue;
        const auto sf=geo.face_area_vectors[f];
        req(std::abs(sf.x)>0.,"V52 layered: target face not x-normal");
        const auto owner=mesh.ownership().owner(f);
        long double watts=0.;
        for(std::size_t m=0;m<dirs.size();++m) {
            const double dot=dirs[m].dx*sf.x+dirs[m].dy*sf.y+dirs[m].dz*sf.z;
            watts+=dirs[m].weight*dot*(dot>=0.?I[m](owner):I[m](static_cast<std::size_t>(neighbour)));
        }
        numerical=static_cast<double>(watts)/sf.x; // +x W/m2
        found=true;break;
    }
    req(found,"V52 layered: central material interface not found");
    const auto discrete=exact_discrete_flux(dirs,k1,k2,B1,B2,L,R);
    const auto exact=exact_continuous_flux(k1,k2,B1,B2,L,R);
    const double scale=std::max(1.,std::abs(exact));
    const Observed errors{std::abs(numerical-discrete)/scale,
                          std::abs(numerical-exact)/scale,
                          std::abs(discrete-exact)/scale,
                          audit.true_relative};
    std::cout<<"V52_LAYER_REAL nx="<<nx<<" polar="<<polar<<" width="<<width
             <<" k1="<<k1<<" k2="<<k2<<" numerical="<<numerical
             <<" oracle_discrete="<<discrete<<" oracle_exact="<<exact
             <<" spatial="<<errors.spatial<<" angular="<<errors.angle
             <<" total="<<errors.total<<" ho="<<errors.ho<<'\n';
    return errors;
}
}
int main(){try{
    const auto a=run(20,8,200.,.25,3.0);
    const auto b=run(40,8,200.,.25,3.0);
    const auto c=run(80,8,200.,.25,3.0);
    req(b.spatial<a.spatial && c.spatial<b.spatial,
        "V52 layered: real DOM spatial error did not contract");
    req(a.spatial<.02 && c.spatial<.005,
        "V52 layered: real DOM spatial accuracy outside qualification gate");
    req(c.ho<1.e-8,"V52 layered: production-FV transport not converged");
    const auto finer_angles=run(40,16,200.,.25,3.0);
    req(finer_angles.angle<b.angle && finer_angles.angle<.002,
        "V52 layered: angular quadrature not improved/accurate");
    const auto strong_a=run(20,8,200.,.03,14.);
    const auto strong_b=run(40,8,200.,.03,14.);
    req(strong_b.spatial<strong_a.spatial && strong_b.spatial<.001,
        "V52 layered: high-contrast spatial error did not contract/qualify");
    const auto wide=run(40,8,400.,.25,3.0);
    req(std::abs(wide.spatial-b.spatial)<.002,
        "V52 layered: unexpected lateral-wall contamination");
    std::cout<<"V52_LAYERED_REAL_TRANSPORT_PASS\n";
}catch(const std::exception& e){std::cerr<<"V52_LAYERED_REAL_TRANSPORT_FAIL: "<<e.what()<<'\n';return 1;}}
