// V48: spatially variable E, heterogeneous LO extinction and exactly
// conservative HO/LO currents on rotated, anisotropic 3D finite volumes.
// The production HO DOM still uses constant optical properties; the variable
// coefficients tested here belong to the opt-in LO study only.
#include "cfdx/physics/dom_smm_second_moment_v48.h"
#include "cfdx/physics/dom_quadrature_v39.h"
#include "thermal_3d_mesh_v36.h"
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
void require(bool result,const char* why) {if(!result)throw std::runtime_error(why);}
Mesh make_grid(std::size_t nx,std::size_t ny,std::size_t nz,double depth) {
    auto mesh=cfdx::testing::thermal_v36::make_mesh(nx,ny,nz,.018,0.,depth);
    constexpr double ax=.41,az=.33;
    for(std::size_t p=0;p<mesh.points().size();++p) {
        const double x=mesh.points().x_data()[p],y=mesh.points().y_data()[p];
        const double z=mesh.points().z_data()[p];
        const double xr=std::cos(ax)*x-std::sin(ax)*y;
        const double yr=std::sin(ax)*x+std::cos(ax)*y;
        mesh.points().set(p,std::cos(az)*xr+std::sin(az)*z,
                            yr,-std::sin(az)*xr+std::cos(az)*z);
    }
    return mesh;
}
using Tensor=DomSmmSecondMomentV48::Tensor;
void manufactured(std::size_t nx,std::size_t ny,std::size_t nz,double depth) {
    const auto mesh=make_grid(nx,ny,nz,depth);
    const auto g=build_fv_geometry(mesh);
    DomSmmSecondMomentV48 smm(mesh,g);
    const std::size_t nc=mesh.n_cells(),nf=mesh.n_faces();
    constexpr double PI4=12.56637061435917295385;
    std::vector<double> J(nc,1.7),sigma(nc);
    std::vector<Tensor> E(nc);
    for(std::size_t c=0;c<nc;++c){
        const double x=g.cell_centres[c].x;
        E[c]={1./3.+.02*x,1./3.-.02*x,1./3.,.04,.02,.01};
        sigma[c]=2.+.2*x;
    }
    // J constant, E variable: div(E J) must contain J*grad(E).
    // An implementation using only E*grad(J) would incorrectly give ZERO.
    std::vector<double> boundary(nf,0.);
    for(std::size_t f=0;f<nf;++f) {
        if(mesh.ownership().neighbour(f)<0) {
            const auto o=mesh.ownership().owner(f);
            boundary[f]=-PI4/sigma[o] * 1.7*.02 * g.face_area_vectors[f].x;
        }
    }
    const auto flux=smm.evaluate(J,E,sigma,boundary);
    double max_error=0.,max_scale=1.,sum=0.,exterior=0.;
    for(std::size_t f=0;f<nf;++f) {
        const auto o=mesh.ownership().owner(f);
        const auto n=mesh.ownership().neighbour(f);
        if(n>=0) {
            const auto j=static_cast<std::size_t>(n);
            const double sh=2./(1./sigma[o]+1./sigma[j]);
            const double exact=-PI4/sh * 1.7*.02*g.face_area_vectors[f].x;
            max_error=std::max(max_error,std::abs(flux.face_watts[f]-exact));
            max_scale=std::max(max_scale,std::abs(exact));
        } else exterior+=flux.face_watts[f];
    }
    for(double v:flux.divergence_watts)sum+=v;
    require(max_error/max_scale<2e-9,"V48: missing variable-E term in div(EJ)");
    require(std::abs(sum-exterior)<1e-9,"V48: nonconservative face current");
    require(std::abs(flux.internal_cancellation_watts)<1e-12,
            "V48: owner-neighbour current failed independent long-double sum");
    // Constant E and affine J -> exact tensor diffusion, including cross terms.
    const Tensor tensor={.45,.35,.20,.065,.025,.014};
    const Vec3 grad={.7,-.4,.2};
    for(std::size_t c=0;c<nc;++c) {
        const auto p=g.cell_centres[c];
        J[c]=2.8+grad.x*p.x+grad.y*p.y+grad.z*p.z;
        E[c]=tensor;
    }
    const Vec3 eg={tensor[0]*grad.x+tensor[3]*grad.y+tensor[4]*grad.z,
                   tensor[3]*grad.x+tensor[1]*grad.y+tensor[5]*grad.z,
                   tensor[4]*grad.x+tensor[5]*grad.y+tensor[2]*grad.z};
    for(std::size_t f=0;f<nf;++f) {
        if(mesh.ownership().neighbour(f)<0) {
            const auto o=mesh.ownership().owner(f);
            boundary[f]=-PI4/sigma[o]*g.face_area_vectors[f].dot(eg);
        }
    }
    const auto affine=smm.evaluate(J,E,sigma,boundary);
    double affine_error=0.,affine_scale=1.;
    for(std::size_t f=0;f<nf;++f) {
        if(mesh.ownership().neighbour(f)<0)continue;
        const auto o=mesh.ownership().owner(f);
        const auto j=static_cast<std::size_t>(mesh.ownership().neighbour(f));
        const double sh=2./(1./sigma[o]+1./sigma[j]);
        const double exact=-PI4/sh*g.face_area_vectors[f].dot(eg);
        affine_error=std::max(affine_error,std::abs(affine.face_watts[f]-exact));
        affine_scale=std::max(affine_scale,std::abs(exact));
    }
    require(affine_error/affine_scale<2e-9,"V48: rotated affine tensor oracle failed");
    // The HO correction MUST be refreshed with J; zeroing it is not allowed.
    std::vector<double> reference=affine.face_watts;
    for(std::size_t f=0;f<nf;++f)reference[f]+=0.05*std::sin(double(f)*.17);
    const auto reconcile=smm.reconcile(affine,reference);
    require(reconcile.max_face_identity_relative<1e-12 &&
            reconcile.max_cell_identity_relative<1e-12,
            "V48: HO/LO face-correction identity failed");
    for(double& j:J)j*=1.02;
    const auto changed=smm.evaluate(J,E,sigma,boundary);
    double delta=0.;
    for(std::size_t f=0;f<nf;++f)
        delta=std::max(delta,std::abs(changed.face_watts[f]-affine.face_watts[f]));
    require(delta>1e-10,"V48: tensor operator accidentally frozen");
    // A second-moment elliptic action is linear for fixed E and optical
    // coefficients even when its Krylov direction contains negative values.
    std::vector<double> v(nc),absorption(nc);
    for(std::size_t c=0;c<nc;++c) {
        v[c]=std::sin(.3+.47*double(c));
        absorption[c]=.07+.01*std::abs(g.cell_centres[c].x);
    }
    const auto frozen=smm.elliptic_action(v,E,sigma,absorption);
    const double h=1e-5;
    auto shifted=J;
    for(std::size_t c=0;c<nc;++c)shifted[c]+=h*v[c];
    const auto shifted_flux=smm.evaluate(shifted,E,sigma,boundary);
    double derivative_error=0.,derivative_scale=1.;
    for(std::size_t c=0;c<nc;++c) {
        const double raw=(shifted_flux.divergence_watts[c]-changed.divergence_watts[c])/h
                         +PI4*absorption[c]*g.cell_volumes[c]*v[c];
        derivative_error=std::max(derivative_error,std::abs(raw-frozen[c]));
        derivative_scale=std::max(derivative_scale,std::abs(frozen[c]));
    }
    require(derivative_error/derivative_scale<5e-8,
            "V48: elliptic Jacobian is not a true frozen-E direction");
    // This case is NOT a proof of coercivity/M-matrix on all polyhedra.
    // Negative tests: fail closed, no input/solution mutation.
    bool refused=false;
    auto invalidE=E; invalidE[0]={1.,0.,0.,0.,0.,0.};
    try {(void)smm.evaluate(J,invalidE,sigma,boundary);}
    catch(const std::invalid_argument&){refused=true;}
    require(refused,"V48: non-SPD second moment accepted");
    refused=false;
    auto invalidSigma=sigma; invalidSigma[0]=0.;
    try {(void)smm.evaluate(J,E,invalidSigma,boundary);}
    catch(const std::invalid_argument&){refused=true;}
    require(refused,"V48: zero extinction accepted");
    refused=false;
    auto invalidBoundary=boundary;
    for(std::size_t f=0;f<nf;++f)
        if(mesh.ownership().neighbour(f)<0) {invalidBoundary[f]=std::numeric_limits<double>::quiet_NaN();break;}
    try {(void)smm.evaluate(J,E,sigma,invalidBoundary);}
    catch(const std::invalid_argument&){refused=true;}
    require(refused,"V48: missing wall-oracle current accepted");
    std::cout<<"V48_SMM_MMS cells="<<nc<<" max_E_gradient_error="<<max_error/max_scale
             <<" affine_tensor_error="<<affine_error/affine_scale
             <<" elliptic_direction="<<derivative_error/derivative_scale
             <<" conservative="<<std::abs(sum-exterior)
             <<" closure_face="<<reconcile.max_face_identity_relative<<'\n';
}
void actual_dom() {
    auto mesh=make_grid(3,3,3,.03);
    auto geo=build_fv_geometry(mesh);
    auto dirs=make_dom_tensor_gauss_legendre_v39(4,8);
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& name=mesh.boundary().patch(p).name;
        walls[name]={name=="bottom"?0.0:name=="west"?.4:.8,
                     name=="west"?500.:300.};
    }
    using Cell=DomMomentSchurSystemV44::CellField;
    Cell temp(mesh.n_cells(),"temperature","K",1);
    for(std::size_t c=0;c<mesh.n_cells();++c)
        temp(c)=350.+25.*geo.cell_centres[c].x;
    RadiationTransportControls ctrl;
    ctrl.absorption=.1;ctrl.scattering=14.9;
    ctrl.tolerance=1e-9;ctrl.linear_tolerance=1e-14;
    ctrl.linear_max_iterations=2000;
    DomMomentSchurSystemV44 t(mesh,geo,temp,dirs,walls,ctrl);
    DomHoloFaceClosureV47 holo(t,mesh,geo,temp,dirs,walls,ctrl);
    std::vector<double> x(t.dimension());
    for(std::size_t c=0;c<mesh.n_cells();++c)x[c]=blackbody_intensity(temp(c));
    for(const auto& wall:t.wall_faces())
        if(wall.reflect_index<t.n_reflecting_faces())
            x[mesh.n_cells()+wall.reflect_index]=wall.ib;
    const auto snap=holo.evaluate(x);
    std::vector<DomSmmSecondMomentV48::Tensor> E(mesh.n_cells());
    std::vector<double> J(mesh.n_cells()),sigma(mesh.n_cells(),15.);
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        J[c]=snap.projected[c];
        for(std::size_t k=0;k<6;++k)E[c][k]=snap.eddington_symmetric[6*c+k];
    }
    std::vector<double> boundary(mesh.n_faces(),0.);
    for(std::size_t f=0;f<mesh.n_faces();++f)
        if(mesh.ownership().neighbour(f)<0)boundary[f]=snap.face_current[f];
    DomSmmSecondMomentV48 smm(mesh,geo);
    const auto lo=smm.evaluate(J,E,sigma,boundary);
    const auto corrected=smm.reconcile(lo,snap.face_current);
    require(corrected.max_face_identity_relative<1e-12 &&
            corrected.max_cell_identity_relative<1e-12,
            "V48: HO DOM currents not recovered on actual gray-wall transport");
    std::cout<<"V48_SMM_DOM cells="<<mesh.n_cells()
             <<" face_error="<<corrected.max_face_identity_relative
             <<" cell_error="<<corrected.max_cell_identity_relative<<'\n';
}
}
int main(int argc,char** argv) {
    try {
        const bool small=argc==2 && std::string(argv[1])=="--sanitizer-small";
        if(argc!=1 && !small)throw std::invalid_argument("V48: unknown arguments");
        manufactured(3,3,3,.03);
        actual_dom();
        if(!small){manufactured(4,4,6,.012);manufactured(8,8,10,.003);}
        std::cout<<"V48_SMM_SECOND_MOMENT PASS\n";
    } catch(const std::exception& e){
        std::cerr<<"V48_SMM_SECOND_MOMENT FAIL: "<<e.what()<<'\n';return 1;
    }
}
