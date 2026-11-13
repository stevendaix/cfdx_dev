// Uses the actual 152-cell BFS mesh from the V54 HDF5 evidence, with explicitly
// SYNTHETIC scalar/flux fields. It is a face-operator test, not a BFS solution.
#include "cfdx/physics/convective_face_audit_v56.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
static void check(bool c,const char* why){if(!c)throw std::runtime_error(why);}
static Mesh load_fixture(const char* file) {
    std::ifstream in(file);
    check(bool(in),"V56 BFS fixture missing");
    std::string title;std::getline(in,title);
    check(title.find("real BFS152 mesh topology")!=std::string::npos,"wrong fixture provenance");
    std::size_t np=0,nf=0,nc=0,npatches=0;
    in>>np>>nf>>nc>>npatches;
    check(np==374 && nf==634 && nc==152 && npatches==8,"wrong BFS152 fixture sizes");
    Mesh m;m.points().resize(np);
    for(std::size_t i=0;i<np;++i){double x,y,z;in>>x>>y>>z;m.points().set(i,x,y,z);}
    m.ownership().resize(nf);
    for(std::size_t f=0;f<nf;++f){
        std::size_t count;in>>count;std::vector<FaceIndex> vertex(count);
        for(auto& v:vertex)in>>v;
        m.faces().push_face(vertex);
        std::size_t owner;std::int64_t neigh;in>>owner>>neigh;
        m.ownership().set_owner(f,owner);
        m.ownership().set_neighbour(f,neigh);
    }
    for(std::size_t c=0;c<nc;++c){
        std::size_t count;in>>count;std::vector<FaceIndex> ids(count);
        for(auto& f:ids)in>>f;
        m.cells().push_cell(ids);
    }
    for(std::size_t p=0;p<npatches;++p){
        std::string name;bool empty;std::size_t count;in>>name>>empty>>count;
        Patch patch;patch.name=name;patch.type=empty?PatchType::EMPTY:PatchType::WALL;
        for(std::size_t j=0;j<count;++j){FaceIndex f;in>>f;patch.face_ids.push_back(f);}
        m.boundary().add_patch(patch);
    }
    check(bool(in),"truncated BFS fixture");
    return m;
}
int main(int argc,char** argv) {
    check(argc==2,"usage: test_bfs152_face_tangent_v56 PATH_TO_FIXTURE");
    auto mesh=load_fixture(argv[1]);
    const auto geo=build_fv_geometry(mesh);
    CellScalarV56 field(mesh.n_cells(),"synthetic","",1),dfield(mesh.n_cells(),"dsynthetic","",1);
    FaceScalarV56 flux(mesh.n_faces(),"synthetic_phi","",1),dflux(mesh.n_faces(),"dphi","",1);
    for(std::size_t c=0;c<mesh.n_cells();++c){
        const auto v=geo.cell_centres[c];
        field(c)=2.0+0.7*std::sin(50*v.x)+0.3*std::cos(130*v.y)+0.07*std::sin(double(c)*0.1);
        dfield(c)=0.1*std::cos(35*v.x)+0.07*std::sin(90*v.y)+0.02*std::sin(double(c)*0.1);
    }
    for(std::size_t f=0;f<mesh.n_faces();++f){
        flux(f)=0.0005*std::sin(static_cast<double>(f)*0.25)+0.0002*std::cos(static_cast<double>(f)*0.11);
        dflux(f)=0.00003*std::sin(static_cast<double>(f)*0.12)+0.00002*std::cos(static_cast<double>(f)*0.19);
    }
    ScalarBoundaryConditions bc;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
        const auto& patch=mesh.boundary().patch(p);
        if(patch.type==PatchType::EMPTY)continue;
        bc[patch.name]={ScalarBoundaryType::FIXED_VALUE,1.7+0.17*static_cast<double>(p),0.0};
    }
    const ConvectionScheme schemes[]={ConvectionScheme::UPWIND,
        ConvectionScheme::SECOND_ORDER_UPWIND,ConvectionScheme::TVD,
        ConvectionScheme::TVD_VANLEER,ConvectionScheme::TVD_VANALBADA,
        ConvectionScheme::MUSCL_CELL_LIMITED};
    std::size_t verified=0,nonsmooth=0,checked_faces=0;
    double max_primal_error=0.,max_tangent_error=0.;
    double max_cell_parity=0.,max_cell_tangent_error=0.,max_additivity=0.;
    CellScalarV56 direction2(mesh.n_cells(),"independent_direction","",1);
    FaceScalarV56 dflux2(mesh.n_faces(),"independent_flux_direction","",1);
    for(std::size_t c=0;c<mesh.n_cells();++c)
        direction2(c)=0.04*std::sin(0.27*static_cast<double>(c))+
                      0.09*std::cos(0.093*static_cast<double>(c));
    for(std::size_t f=0;f<mesh.n_faces();++f)
        dflux2(f)=0.00002*std::cos(0.13*static_cast<double>(f));
    auto direction_sum=dfield; auto flux_sum=dflux;
    for(std::size_t c=0;c<mesh.n_cells();++c)direction_sum(c)+=direction2(c);
    for(std::size_t f=0;f<mesh.n_faces();++f)flux_sum(f)+=dflux2(f);
    for(auto scheme:schemes){
        auto audited=audit_scalar_convective_faces_v56(mesh,geo,field,dfield,flux,dflux,scheme,bc);
        auto audited2=audit_scalar_convective_faces_v56(mesh,geo,field,direction2,flux,dflux2,scheme,bc);
        auto audited_sum=audit_scalar_convective_faces_v56(mesh,geo,field,direction_sum,flux,flux_sum,scheme,bc);
        auto original=reconstruct_scalar_transport_flux(mesh,geo,flux,field,0.0,bc,true,scheme);
        const double h=1.e-6;
        auto fp=field,fm=field; auto plus=flux,minus=flux;
        for(std::size_t c=0;c<mesh.n_cells();++c){fp(c)+=h*dfield(c);fm(c)-=h*dfield(c);}
        for(std::size_t f=0;f<mesh.n_faces();++f){plus(f)+=h*dflux(f);minus(f)-=h*dflux(f);}
        auto yp=reconstruct_scalar_transport_flux(mesh,geo,plus,fp,0.0,bc,true,scheme);
        auto ym=reconstruct_scalar_transport_flux(mesh,geo,minus,fm,0.0,bc,true,scheme);
        for(std::size_t f=0;f<mesh.n_faces();++f){
            max_primal_error=std::max(max_primal_error,std::abs(audited[f].flux-original(f)));
            if(audited[f].nonsmooth || audited2[f].nonsmooth ||
               audited_sum[f].nonsmooth){++nonsmooth;continue;}
            max_additivity=std::max(max_additivity,std::abs(
               audited[f].tangent+audited2[f].tangent-audited_sum[f].tangent));
            const double fd=(yp(f)-ym(f))/(2*h);
            const double err=std::abs(audited[f].tangent-fd);
            max_tangent_error=std::max(max_tangent_error,err);
            check(err<=2.e-7*std::max(1.0,std::abs(fd)),"V56 tangent differs from real-mesh production FD");
            ++verified;
        }
        // Compare the complete bounded-convection equation, including
        // -div(phi)*psi, against the original CSR reassembly at each candidate.
        auto residual=scatter_bounded_convection_v56(mesh,geo,field,dfield,
            flux,dflux,audited,true);
        const auto nc=mesh.n_cells();
        CellScalarV56 zeros(nc,"source_zero","",1);
        zeros.fill(0.0);
        std::vector<double> extra(nc,1.0); // avoid zero diagonals; subtract below
        auto eq=assemble_scalar_equation(mesh,geo,flux,0.0,zeros,zeros,bc,
            true,nullptr,&extra,nullptr,nullptr,scheme,&field);
        auto ep=assemble_scalar_equation(mesh,geo,plus,0.0,zeros,zeros,bc,
            true,nullptr,&extra,nullptr,nullptr,scheme,&fp);
        auto em=assemble_scalar_equation(mesh,geo,minus,0.0,zeros,zeros,bc,
            true,nullptr,&extra,nullptr,nullptr,scheme,&fm);
        const auto mat_res=[](const ScalarEquation& a,const CellScalarV56& u,
                              std::size_t row) {
            double r=-a.rhs(row);
            for(auto k=a.matrix.row_offsets_data()[row];
                k<a.matrix.row_offsets_data()[row+1];++k)
                r+=a.matrix.values_data()[k]*u(a.matrix.columns_data()[k]);
            return r;
        };
        for(std::size_t c=0;c<nc;++c){
            const double original_res=mat_res(eq,field,c)-extra[c]*field(c);
            const double derivative=(mat_res(ep,fp,c)-extra[c]*fp(c)
                -(mat_res(em,fm,c)-extra[c]*fm(c)))/(2*h);
            max_cell_parity=std::max(max_cell_parity,
                std::abs(original_res-residual.value[c]));
            max_cell_tangent_error=std::max(max_cell_tangent_error,
                std::abs(derivative-residual.tangent[c]));
        }
        ++checked_faces;
    }
    check(max_primal_error<5.e-12,"V56 primal not equal to production face flux");
    check(verified>1000,"V56 insufficient smooth faces verified");
    check(max_cell_parity<1.e-10,"V56 bounded assembly primal parity failed");
    check(max_cell_tangent_error<2.e-6,"V56 bounded assembly tangent parity failed");
    check(max_additivity<1.e-11,"V56 active-set tangent is not additive");
    // On the REAL BFS topology, force a nominal zero-flux internal face.
    // Neither a one-sided branch derivative nor a smoothed surrogate may be
    // certified as an ordinary two-sided derivative at this kink.
    auto zero_crossing_flux=flux;
    std::size_t crossing=mesh.n_faces();
    for(std::size_t f=0;f<mesh.n_faces();++f)
        if(mesh.ownership().neighbour(f)>=0 && dflux(f)!=0.0){crossing=f;break;}
    check(crossing<mesh.n_faces(),"BFS fixture has no internal branch face");
    zero_crossing_flux(crossing)=0.0;
    auto at_kink=audit_scalar_convective_faces_v56(mesh,geo,field,dfield,
        zero_crossing_flux,dflux,ConvectionScheme::TVD,bc);
    check(at_kink[crossing].nonsmooth,"V56 ignored zero-flux branch change");
    std::cout<<std::setprecision(12)
        <<"V56_BFS152_FACE_AUDIT: schemes="<<checked_faces
        <<" faces="<<mesh.n_faces()<<" verified="<<verified
        <<" nonsmooth="<<nonsmooth
        <<" max_primal_error="<<max_primal_error
        <<" max_tangent_error="<<max_tangent_error
        <<" max_cell_parity="<<max_cell_parity
        <<" max_cell_tangent_error="<<max_cell_tangent_error
        <<" max_additivity="<<max_additivity
        <<" kink_face_id="<<crossing<<" kink_detected=1"<<"\n";
}
