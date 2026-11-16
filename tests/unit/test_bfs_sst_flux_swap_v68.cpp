#include "cfdx/physics/bfs_sst_flux_swap_v68.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
using namespace cfdx::core;
static void require(bool yes,const char* reason){if(!yes)throw std::runtime_error(reason);}
template<class F>static void rejects(F&& f){bool threw=false;try{f();}catch(const std::exception&){threw=true;}require(threw,"V68 must reject falsified evidence");}
static Mesh cube(){
    Mesh mesh;mesh.points().resize(8);
    constexpr double pts[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                                {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i)mesh.points().set(i,pts[i][0],pts[i][1],pts[i][2]);
    const std::vector<std::vector<FaceIndex>> faces={
        {0,3,2,1},{4,5,6,7},{0,1,5,4},{3,7,6,2},{0,4,7,3},{1,2,6,5}};
    for(const auto& f:faces)mesh.faces().push_face(f);
    mesh.ownership().resize(faces.size());
    for(std::size_t f=0;f<faces.size();++f){mesh.ownership().set_owner(f,0);mesh.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);}
    mesh.cells().push_cell({0,1,2,3,4,5});
    Patch patch;patch.name="wall";patch.type=PatchType::WALL;
    patch.face_ids={0,1,2,3,4,5};mesh.boundary().add_patch(patch);
    return mesh;
}
static Mesh two_cubes(){
    Mesh mesh;mesh.points().resize(12);
    constexpr double pts[12][3]={{0,0,0},{0,1,0},{0,1,1},{0,0,1},
                                  {1,0,0},{1,1,0},{1,1,1},{1,0,1},
                                  {2,0,0},{2,1,0},{2,1,1},{2,0,1}};
    for(std::size_t i=0;i<12;++i)mesh.points().set(i,pts[i][0],pts[i][1],pts[i][2]);
    const std::vector<std::vector<FaceIndex>> faces={
        {4,5,6,7}, // shared face; owner 0, neighbour 1
        {0,3,2,1},{0,4,7,3},{1,2,6,5},{0,1,5,4},{3,7,6,2},
        {8,9,10,11},{4,8,11,7},{5,6,10,9},{4,5,9,8},{7,11,10,6}};
    for(const auto& face:faces)mesh.faces().push_face(face);
    mesh.ownership().resize(faces.size());
    for(std::size_t f=0;f<faces.size();++f){
        mesh.ownership().set_owner(f,f<=5?0:1);
        mesh.ownership().set_neighbour(f,f==0?1:FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0,1,2,3,4,5});
    mesh.cells().push_cell({0,6,7,8,9,10});
    Patch patch;patch.name="wall";patch.type=PatchType::WALL;
    patch.face_ids={1,2,3,4,5,6,7,8,9,10};mesh.boundary().add_patch(patch);
    return mesh;
}
int main(){
    auto mesh=cube();auto geom=build_fv_geometry(mesh);
    Field<double,Location::FACE> carrier(6,"carrier","kg/s",1),physical(6,"physical","kg/s",1);
    carrier.fill(0.);physical.fill(0.);physical(0)=-0.1;physical(1)=0.2;
    Field<double,Location::CELL> k(1,"k","",1),omega(1,"omega","",1),
        strain(1,"strain","",1),distance(1,"distance","",1);
    k(0)=.4;omega(0)=8.;strain(0)=2.;distance(0)=.1;
    TurbulenceTransportControls c;c.model=TurbulenceModel::SST;
    ScalarBoundaryConditions kb,wb;
    kb["wall"]={ScalarBoundaryType::FIXED_VALUE,.01,0.};
    wb["wall"]={ScalarBoundaryType::FIXED_VALUE,2.,0.};
    auto eval=[&](const auto& flux){return sst_coupled_residual_v22(mesh,geom,flux,k,omega,strain,distance,c,kb,wb,nullptr);};
    const auto a=eval(carrier),b=eval(physical);
    const auto audit=analyze_bfs_sst_flux_swap_v68(mesh,geom,carrier,physical,k,omega,kb,wb,nullptr,nullptr,a,b,0);
    require(audit.cells==1 && audit.faces==6 && audit.omega_delta_linf>0.,"V68 physical swap not detected");
    require(audit.max_face_identity_gap<1e-12,"V68 face identity");
    require(std::abs((b.values[1]-a.values[1])-audit.face_omega_delta[0])<1e-11,"V68 focus face contribution");
    require(std::abs(audit.face_omega_delta[0]-0.6)<1e-12,"V68 inflow Dirichlet upwind sign");
    const auto same=analyze_bfs_sst_flux_swap_v68(mesh,geom,carrier,carrier,k,omega,kb,wb,nullptr,nullptr,a,a,0);
    require(same.omega_delta_linf==0. && same.k_delta_linf==0.,"V68 identity flux");
    rejects([&]{auto fake=b;fake.values[1]+=1.;analyze_bfs_sst_flux_swap_v68(mesh,geom,carrier,physical,k,omega,kb,wb,nullptr,nullptr,a,fake,0);});
    rejects([&]{auto fake=physical;fake(0)=INFINITY;analyze_bfs_sst_flux_swap_v68(mesh,geom,carrier,fake,k,omega,kb,wb,nullptr,nullptr,a,b,0);});
    rejects([&]{analyze_bfs_sst_flux_swap_v68(mesh,geom,carrier,physical,k,omega,kb,wb,nullptr,nullptr,a,b,8);});
    {
        // An internal face whose flux reverses direction is the decisive
        // conservative upwind test. No modified production SST sources.
        const auto m=two_cubes();const auto g=build_fv_geometry(m);
        Field<double,Location::FACE> fc(m.n_faces(),"fc","kg/s",1),fp(m.n_faces(),"fp","kg/s",1);
        fc.fill(0.);fp.fill(0.);fc(0)=-.1;fp(0)=.15;
        Field<double,Location::CELL> kk(2,"kk","",1),ww(2,"ww","",1),
            ss(2,"ss","",1),dd(2,"dd","",1);
        kk(0)=.4;kk(1)=.8;ww(0)=8.;ww(1)=16.;
        ss.fill(2.);dd.fill(.1);
        const auto rr=[&](const auto& phi){return sst_coupled_residual_v22(m,g,phi,kk,ww,ss,dd,c,kb,wb,nullptr);};
        const auto before=rr(fc),after=rr(fp);
        const auto contrast=analyze_bfs_sst_flux_swap_v68(m,g,fc,fp,kk,ww,kb,wb,nullptr,nullptr,before,after,0);
        require(contrast.omega_delta_linf>0. && contrast.max_face_identity_gap<1e-12,
                "V68 internal upwind face swap FV identity");
        require(std::abs(contrast.face_omega_delta[0]-.8)<1e-12 &&
                std::abs(contrast.face_omega_delta[1]-1.2)<1e-12,
                "V68 owner/neighbour contributions for reversed flux");
    }
    std::cout<<"V68_SST_FLUX_SWAP_PASS identity="<<audit.max_face_identity_gap
             <<" omega_delta="<<audit.omega_delta_linf<<"\n";
}