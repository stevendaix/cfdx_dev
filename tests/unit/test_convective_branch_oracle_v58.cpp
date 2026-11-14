#include "cfdx/physics/convective_branch_oracle_v58.h"
#include "bfs152_fixture_v58.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::core;
using namespace cfdx::physics;
static void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
int main(int argc,char** argv){
    check(argc==2,"V58 BFS152 fixture path required");
    const auto mesh=load_bfs152_fixture_v58(argv[1]);
    const auto geo=build_fv_geometry(mesh);
    const auto nc=mesh.n_cells(),nf=mesh.n_faces();
    CellScalarV56 u(nc,"u","",1),du(nc,"du","",1);
    FaceScalarV56 phi(nf,"phi","",1),dphi(nf,"dphi","",1);
    u.fill(2.);du.fill(0.);phi.fill(0.);dphi.fill(0.);
    std::size_t internal=nf;
    for(std::size_t f=0;f<nf;++f)
        if(mesh.ownership().neighbour(f)>=0){internal=f;break;}
    check(internal<nf,"V58 BFS152 has no internal face");
    const auto neighbour=mesh.ownership().neighbour(internal);
    u(static_cast<std::size_t>(neighbour))=5.;
    dphi(internal)=1.;
    ScalarBoundaryConditions bc;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
        const auto& patch=mesh.boundary().patch(p);
        if(patch.type==PatchType::EMPTY)continue;
        bc[patch.name]={ScalarBoundaryType::ZERO_GRADIENT,0.,0.};
    }
    const auto kink=audit_convective_branches_v58(mesh,geo,u,du,phi,dphi,
        ConvectionScheme::UPWIND,bc);
    check(kink.effective_kinks==1 && kink.faces[internal].effective_kink &&
          !kink.bilateral_certifiable && kink.faces[internal].sided_relative>0.5,
          "V58 missed real BFS152 owner/neighbour upwind kink");
    check(std::abs(kink.faces[internal].right-2.)<1e-8 &&
          std::abs(kink.faces[internal].left-5.)<1e-8,
          "V58 left/right slopes disagree with upwind production");
    // A sign transition is harmless if both upwind states coincide.
    u(static_cast<std::size_t>(neighbour))=2.;
    const auto inert=audit_convective_branches_v58(mesh,geo,u,du,phi,dphi,
        ConvectionScheme::UPWIND,bc);
    check(inert.effective_kinks==0 && inert.bilateral_certifiable &&
          inert.inert_branch_flags>0,"V58 false positive from inert sign switch");
    phi(internal)=3.;
    const auto smooth=audit_convective_branches_v58(mesh,geo,u,du,phi,dphi,
        ConvectionScheme::UPWIND,bc);
    check(smooth.effective_kinks==0 && smooth.bilateral_certifiable,
          "V58 smooth face flagged nonlinear");
    bool invalid=false;
    try{(void)audit_convective_branches_v58(mesh,geo,u,du,phi,dphi,
        ConvectionScheme::UPWIND,bc,0.0);}catch(const std::invalid_argument&){invalid=true;}
    check(invalid,"V58 accepted zero difference step");
    // The complete limited reconstruction must also be assessed. A constant
    // field must not produce a physical kink for any of the TVD/MUSCL schemes.
    const ConvectionScheme high[]={ConvectionScheme::TVD,
        ConvectionScheme::TVD_VANLEER,ConvectionScheme::TVD_VANALBADA,
        ConvectionScheme::MUSCL_CELL_LIMITED};
    for(const auto scheme:high){
        const auto a=audit_convective_branches_v58(mesh,geo,u,du,phi,dphi,
            scheme,bc);
        check(a.bilateral_certifiable,"V58 false TVD/MUSCL kink at constant state");
    }
    // SOU triple tie: upper, lower and extrapolated face value coincide.
    // The physical B-derivative is positively homogeneous but NON-ADDITIVE:
    // neither a centered directional difference nor an arbitrary symmetric
    // selection can be silently used as a linear Newton operator at this tie.
    const auto sou_direction=[](double d_up,double d_down,double d_inc){
        const double t=1e-6;
        const auto eval=[&](double sign){
            return limited_high_face_v56(LimitedFaceKindV56::SOU,
                face_jet_v56(sign*t*d_up),face_jet_v56(sign*t*d_down),
                face_jet_v56(sign*t*d_inc),face_jet_v56(0.0),
                face_jet_v56(0.0),face_jet_v56(0.0)).value;
        };
        return (eval(1.0)-eval(-1.0))/(2.0*t);
    };
    const double first=sou_direction(0.,1.,1.);
    const double second=sou_direction(0.,1.,-1.);
    const double summed=sou_direction(0.,2.,0.);
    check(std::abs(first-1.0)<1e-12 && std::abs(second)<1e-12 &&
          std::abs(summed)<1e-12 &&
          std::abs((first+second)-summed)>0.5,
          "V58 failed to detect genuinely non-additive SOU triple tie");
    // Actual BFS152 geometry with an initially uniform streamwise velocity
    // and no-slip wall data. This is the branch configuration on which the
    // SOU limiter can clamp a nonzero reconstructed wall gradient to two
    // exactly tied neighbouring cell values.
    u.fill(44.2); phi.fill(0.0042774108); dphi.fill(0.0);
    for(std::size_t c=0;c<nc;++c)du(c)=0.02*std::sin(0.17*static_cast<double>(c));
    for(auto& item:bc)item.second={ScalarBoundaryType::FIXED_VALUE,0.,0.};
    auto fv=audit_scalar_convective_faces_v56(mesh,geo,u,du,phi,dphi,
        ConvectionScheme::SECOND_ORDER_UPWIND,bc);
    const auto original=fv;
    const auto repaired=apply_symmetric_sou_ties_v58(mesh,geo,u,du,phi,dphi,
        bc,1.0,fv);
    check(repaired>0,"V58 no generalized clipped SOU faces on BFS topology");
    auto xp=u,xm=u;
    const double h=1e-6;
    for(std::size_t c=0;c<nc;++c){xp(c)+=h*du(c);xm(c)-=h*du(c);}
    const auto fp=reconstruct_scalar_transport_flux(mesh,geo,phi,xp,0.,bc,true,
        ConvectionScheme::SECOND_ORDER_UPWIND);
    const auto fm=reconstruct_scalar_transport_flux(mesh,geo,phi,xm,0.,bc,true,
        ConvectionScheme::SECOND_ORDER_UPWIND);
    double original_error=0.,generalized_error=0.;
    std::size_t measured=0;
    for(std::size_t f=0;f<nf;++f){
        check(fv[f].flux==original[f].flux,
              "V58 generalized selection changed the primal face flux");
        if(fv[f].tangent==original[f].tangent)continue;
        const double central=(fp(f)-fm(f))/(2*h);
        original_error=std::max(original_error,
            std::abs(original[f].tangent-central));
        generalized_error=std::max(generalized_error,
            std::abs(fv[f].tangent-central));
        ++measured;
    }
    std::cerr<<"V58_TEST_GENERALIZED measured="<<measured<<" repaired="<<repaired
             <<" original_error="<<original_error
             <<" generalized_error="<<generalized_error<<"\n";
    check(measured>0 && generalized_error<original_error &&
          generalized_error<1e-6,
          "V58 generalized SOU action not closer to centered physical flux");
    std::cout<<"V58_BFS152_BRANCH_PASS faces="<<nf
             <<" kink_face="<<internal
             <<" right="<<kink.faces[internal].right
             <<" left="<<kink.faces[internal].left
             <<" inert="<<inert.inert_branch_flags
             <<" generalized_faces="<<repaired
             <<" original_gap="<<original_error
             <<" generalized_gap="<<generalized_error<<"\n";
}
