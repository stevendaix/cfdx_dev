#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/low_storage_time_integration.h"
#include "verification_metrics.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::verification;

namespace {
struct Mesh1D { Mesh mesh; };

Mesh1D make_channel(std::size_t n, double H)
{
    if(n<2 || !(H>0)) throw std::invalid_argument("invalid 1D mesh");
    Mesh m; m.points().resize(8*n);
    const double dy=H/n;
    for(std::size_t i=0;i<n;++i){
        const double y0=i*dy,y1=(i+1)*dy; const std::size_t b=8*i;
        const double p[8][3]={{0,y0,0},{1,y0,0},{1,y1,0},{0,y1,0},
                              {0,y0,1},{1,y0,1},{1,y1,1},{0,y1,1}};
        for(std::size_t j=0;j<8;++j)m.points().set(b+j,p[j][0],p[j][1],p[j][2]);
    }
    std::vector<std::size_t> lo,hi,x0,x1,z0,z1,internal;
    auto face=[&](std::initializer_list<std::size_t> v){auto id=m.n_faces();m.faces().push_face(std::vector<FaceIndex>(v));return id;};
    lo.push_back(face({0,1,5,4}));
    hi.push_back(face({8*(n-1)+3,8*(n-1)+7,8*(n-1)+6,8*(n-1)+2}));
    for(std::size_t i=0;i<n;++i){auto b=8*i;x0.push_back(face({b,b+4,b+7,b+3}));x1.push_back(face({b+1,b+2,b+6,b+5}));z0.push_back(face({b,b+3,b+2,b+1}));z1.push_back(face({b+4,b+5,b+6,b+7}));}
    for(std::size_t i=0;i+1<n;++i){auto b=8*i;internal.push_back(face({b+3,b+7,b+6,b+2}));}
    std::vector<std::vector<std::size_t>> cf(n);
    for(std::size_t i=0;i<n;++i) cf[i]={i?internal[i-1]:lo[0],i+1<n?internal[i]:hi[0],x0[i],x1[i],z0[i],z1[i]};
    m.ownership().resize(m.n_faces());

    // Boundary faces have a single owner. Internal faces have a fixed owner
    // equal to the lower cell index and a neighbour equal to owner+1. The
    // previous all-cells scan could overwrite an internal face owner while
    // visiting its neighbour, corrupting the mesh connectivity and causing
    // the MMS executable to access invalid cell topology.
    for(std::size_t c=0;c<n;++c) {
        m.ownership().set_owner(lo[0],0);
        m.ownership().set_owner(hi[0],n-1);
        m.ownership().set_owner(x0[c],c);
        m.ownership().set_owner(x1[c],c);
        m.ownership().set_owner(z0[c],c);
        m.ownership().set_owner(z1[c],c);
    }
    m.ownership().set_neighbour(lo[0],FaceOwnership::BOUNDARY);
    m.ownership().set_neighbour(hi[0],FaceOwnership::BOUNDARY);
    for(std::size_t c=0;c<n;++c) {
        m.ownership().set_neighbour(x0[c],FaceOwnership::BOUNDARY);
        m.ownership().set_neighbour(x1[c],FaceOwnership::BOUNDARY);
        m.ownership().set_neighbour(z0[c],FaceOwnership::BOUNDARY);
        m.ownership().set_neighbour(z1[c],FaceOwnership::BOUNDARY);
    }
    for(std::size_t i=0;i+1<n;++i) {
        m.ownership().set_owner(internal[i],i);
        m.ownership().set_neighbour(internal[i],static_cast<int>(i+1));
    }

    for(auto& q:cf)m.cells().push_cell(q);
    auto patch=[&](const char* name,const std::vector<std::size_t>& ids){Patch p;p.name=name;p.type=PatchType::WALL;p.face_ids=ids;m.boundary().add_patch(p);};
    patch("inlet",lo);patch("outlet",hi);patch("x0",x0);patch("x1",x1);patch("z0",z0);patch("z1",z1);
    return {std::move(m)};
}

struct Result { std::vector<double> y,u,v; };

Result solve_conv_diff(std::size_t n,double H,double D,double V)
{
    auto p=make_channel(n,H);auto g=build_fv_geometry(p.mesh);
    Field<double,Location::FACE> flux(p.mesh.n_faces(),"phi","kg/s",1);flux.fill(0);
    // Positive y mass flux on every y-face. Boundary signs are handled by FV ownership.
    for(std::size_t f=0;f<p.mesh.n_faces();++f) {
        const auto Sf=g.face_area_vectors[f];
        if(std::abs(Sf.y)>0.5) flux(f)=V*Sf.y;
    }
    Field<double,Location::CELL> su(n,"source","1",1),sp(n,"sp","1",1),field(n,"field","1",1);
    const double a=V/D;
    const double pi=std::acos(-1.0);
    for(std::size_t c=0;c<n;++c){const double y=g.cell_centres[c].y;field(c)=std::sin(pi*y/H);su(c)=D*pi*pi/(H*H)*std::sin(pi*y/H)+V*pi/H*std::cos(pi*y/H);}
    ScalarBoundaryConditions bc;
    bc["inlet"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    bc["outlet"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    bc["x0"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};bc["x1"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    bc["z0"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};bc["z1"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    auto eq=assemble_scalar_equation(p.mesh,g,flux,D,su,sp,bc,true,nullptr,nullptr,nullptr,nullptr,ConvectionScheme::UPWIND,&field);
    Vector sol(n,0);auto r=solve_scalar_equation(eq,sol,{5000,1e-12,1});
    if(r.status!=SolverStatus::CONVERGED)throw std::runtime_error("convection-diffusion solve failed");
    Result out;out.y.resize(n);out.u.resize(n);out.v.resize(n);
    for(std::size_t i=0;i<n;++i){out.y[i]=g.cell_centres[i].y;out.u[i]=sol(i);out.v[i]=std::sin(pi*out.y[i]/H);}
    return out;
}

void temporal_convergence()
{
    std::vector<double> errors;
    for(int level=0;level<4;++level){
        const double dt=0.1/std::pow(2.0,level); const int steps=static_cast<int>(1.0/dt);
        Field<double,Location::CELL> u(1,"u");u(0)=1.0;
        for(int i=0;i<steps;++i)low_storage_rk2_step(u,dt,[](const auto& x,auto& r){r(0)=-x(0);});
        const double e=std::abs(u(0)-std::exp(-1.0));errors.push_back(e);
        std::cout<<"RK2_TEMPORAL dt="<<dt<<" error="<<e;
        if(level)std::cout<<" order="<<std::log(errors[level-1]/e)/std::log(2.0);
        std::cout<<"\n";
    }
    if(errors.back()>=errors.front()/20.0)throw std::runtime_error("RK2 temporal refinement insufficient");
    std::cout<<"TEMPORAL_RK2_CONVERGENCE: PASS\n";
}

}

int main()
{
    try{
        const double H=1,D=0.01,V=0.1;
        std::vector<double> errs;
        for(std::size_t n:{16u,32u,64u,128u}){
            auto r=solve_conv_diff(n,H,D,V);
            auto e=error_norms(r.u,r.v,std::vector<double>(n,1.0/n));
            errs.push_back(e.l2);
            std::cout<<"CONVECTION_DIFFUSION_MMS N="<<n<<" L2="<<e.l2<<" Linf="<<e.linf;
            if(errs.size()>1)std::cout<<" order="<<std::log(errs[errs.size()-2]/errs.back())/std::log(2.0);
            std::cout<<"\n";
        }
        if(!(errs.back()<errs.front()))throw std::runtime_error("convection-diffusion refinement did not reduce error");
        temporal_convergence();
        std::cout<<"DISCRETIZATION_VERIFICATION: PASS\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"DISCRETIZATION_VERIFICATION: FAIL: "<<e.what()<<"\n";return 1;}
}
