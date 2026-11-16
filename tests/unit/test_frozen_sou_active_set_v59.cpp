#include "cfdx/physics/frozen_sou_active_set_v59.h"
#include "bfs152_fixture_v58.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::core;
using namespace cfdx::physics;
static void require(bool b,const char* why){if(!b)throw std::runtime_error(why);}
int main(int argc,char** argv) {
    require(argc==2,"V59 needs the real BFS152 topology fixture");
    const auto mesh=load_bfs152_fixture_v58(argv[1]);
    const auto geo=build_fv_geometry(mesh);
    const auto nc=mesh.n_cells(),nf=mesh.n_faces();
    CellScalarV56 u(nc,"base_u","m/s",1),a(nc,"a","m/s",1),b(nc,"b","m/s",1),sum(nc,"sum","m/s",1);
    FaceScalarV56 phi(nf,"phi","kg/s",1),da(nf,"da","kg/s",1),db(nf,"db","kg/s",1),ds(nf,"ds","kg/s",1);
    u.fill(44.2); phi.fill(0.0042774108);
    for(std::size_t c=0;c<nc;++c) {
        a(c)=0.05*std::sin(0.19*static_cast<double>(c));
        b(c)=0.07*std::cos(0.31*static_cast<double>(c));
        sum(c)=a(c)+b(c);
    }
    for(std::size_t f=0;f<nf;++f) {
        da(f)=0.00001*std::sin(0.07*static_cast<double>(f));
        db(f)=0.00003*std::cos(0.13*static_cast<double>(f));
        ds(f)=da(f)+db(f);
    }
    ScalarBoundaryConditions bc;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& patch=mesh.boundary().patch(p);
        if(patch.type==PatchType::EMPTY)continue;
        bc[patch.name]={ScalarBoundaryType::FIXED_VALUE,0.,0.};
    }
    const auto frozen=freeze_sou_active_set_v59(mesh,geo,u,phi,bc);
    require(frozen.triple_ties>0,"V59 lost the uniform-U BFS152 ties");
    const auto af=apply_frozen_sou_v59(frozen,mesh,geo,u,a,phi,da,bc);
    const auto bf=apply_frozen_sou_v59(frozen,mesh,geo,u,b,phi,db,bc);
    const auto sf=apply_frozen_sou_v59(frozen,mesh,geo,u,sum,phi,ds,bc);
    const auto primal=reconstruct_scalar_transport_flux(mesh,geo,phi,u,0.,bc,true,
        ConvectionScheme::SECOND_ORDER_UPWIND);
    double max_gap=0.,max_primal=0.;
    std::size_t internal=0,nonunique=0;
    for(std::size_t f=0;f<nf;++f) {
        max_gap=std::max(max_gap,std::abs(sf[f].tangent-af[f].tangent-bf[f].tangent));
        max_primal=std::max(max_primal,std::abs(sf[f].flux-primal(f)));
        if(mesh.ownership().neighbour(f)>=0) {
            ++internal;
            nonunique+=sf[f].nonsmooth?1u:0u;
        }
    }
    require(max_gap<1e-11,"V59 Krylov face action not additive");
    require(max_primal<1e-10,"V59 changed the nonlinear primal face flux");
    require(nonunique==frozen.nonunique_internal_faces,
        "V59 nondifferentiable faces falsely disappeared");
    // Equal and opposite internal face scatter -- independent of the Jacobian
    // linearity test and of artificial bounded-convection diagonal terms.
    const auto ar=scatter_bounded_convection_v56(mesh,geo,u,a,phi,da,af,false);
    const auto br=scatter_bounded_convection_v56(mesh,geo,u,b,phi,db,bf,false);
    const auto sr=scatter_bounded_convection_v56(mesh,geo,u,sum,phi,ds,sf,false);
    double max_row_gap=0.;
    for(std::size_t c=0;c<nc;++c)
        max_row_gap=std::max(max_row_gap,
            std::abs(sr.tangent[c]-ar.tangent[c]-br.tangent[c]));
    require(max_row_gap<1e-10,"V59 conservative cell scatter not additive");
    // The same selection must not be silently reused after either the cell
    // field or Rhie-Chow face-flux base has changed.
    bool rejected=false;
    auto mutated=u;mutated(0)+=1e-7;
    try {(void)apply_frozen_sou_v59(frozen,mesh,geo,mutated,a,phi,da,bc);}
    catch(const std::invalid_argument&) {rejected=true;}
    require(rejected,"V59 reused stale U active set");
    rejected=false;
    auto changed_phi=phi;changed_phi(0)+=1e-7;
    try {(void)apply_frozen_sou_v59(frozen,mesh,geo,u,a,changed_phi,da,bc);}
    catch(const std::invalid_argument&) {rejected=true;}
    require(rejected,"V59 reused stale Rhie-Chow flux active set");
    // For strictly smooth faces, the frozen derivative MUST also match the
    // independently recomputed production primal: linearity alone is not
    // sufficient evidence of a correct Jacobian.
    CellScalarV56 smooth_u=u;
    for(std::size_t c=0;c<nc;++c)
        smooth_u(c)=1.0+0.01*static_cast<double>(c)+
                    0.1*std::sin(0.37*static_cast<double>(c));
    const auto regular=freeze_sou_active_set_v59(mesh,geo,smooth_u,phi,bc);
    const auto regular_faces=apply_frozen_sou_v59(regular,mesh,geo,
        smooth_u,a,phi,da,bc);
    auto plus_u=smooth_u,minus_u=smooth_u;
    auto plus_flux=phi,minus_flux=phi;
    constexpr double h=1e-6;
    for(std::size_t c=0;c<nc;++c) {
        plus_u(c)+=h*a(c);
        minus_u(c)-=h*a(c);
    }
    for(std::size_t f=0;f<nf;++f) {
        plus_flux(f)+=h*da(f);
        minus_flux(f)-=h*da(f);
    }
    const auto plus_primal=reconstruct_scalar_transport_flux(
        mesh,geo,plus_flux,plus_u,0.,bc,true,
        ConvectionScheme::SECOND_ORDER_UPWIND);
    const auto minus_primal=reconstruct_scalar_transport_flux(
        mesh,geo,minus_flux,minus_u,0.,bc,true,
        ConvectionScheme::SECOND_ORDER_UPWIND);
    std::size_t smooth_faces=0;
    double maximum_smooth_fidelity=0.;
    for(std::size_t f=0;f<nf;++f) {
        if(mesh.ownership().neighbour(f)<0 || regular.faces[f].nonunique)
            continue;
        const double independent=(plus_primal(f)-minus_primal(f))/(2.*h);
        maximum_smooth_fidelity=std::max(maximum_smooth_fidelity,
            std::abs(independent-regular_faces[f].tangent));
        ++smooth_faces;
    }
    require(smooth_faces>0,"V59 smooth face validation had no sample");
    require(maximum_smooth_fidelity<1e-7,
        "V59 frozen SOU failed on genuinely regular production faces");
    std::cout<<"V59_FROZEN_SOU_PASS cells="<<nc<<" internal="<<internal
        <<" triple_ties="<<frozen.triple_ties
        <<" unresolved_physical_ties="<<frozen.nonunique_internal_faces
        <<" max_additivity="<<max_gap
        <<" max_scatter_additivity="<<max_row_gap
        <<" max_primal_error="<<max_primal
        <<" smooth_faces="<<smooth_faces
        <<" max_smooth_fidelity="<<maximum_smooth_fidelity
        <<" newton_qualified=0\n";
}
