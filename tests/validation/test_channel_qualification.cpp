#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/io/hdf5/hdf5_reader.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

Mesh make_channel_mesh(std::size_t nx, std::size_t ny, double length)
{
    if (nx < 4 || ny < 4 || !(length > 1.0))
        throw std::invalid_argument("invalid channel mesh");

    Mesh mesh;
    const std::size_t plane = (nx + 1) * (ny + 1);
    mesh.points().resize(2 * plane);
    const auto id=[nx](std::size_t i,std::size_t j,std::size_t k){
        return (j*(nx+1)+i)*2+k;
    };
    for(std::size_t j=0;j<=ny;++j)
        for(std::size_t i=0;i<=nx;++i) {
            const double x=length*static_cast<double>(i)/static_cast<double>(nx);
            const double y=static_cast<double>(j)/static_cast<double>(ny);
            mesh.points().set(id(i,j,0),x,y,0.0);
            mesh.points().set(id(i,j,1),x,y,1.0);
        }

    std::map<std::vector<std::size_t>,std::size_t> faces;
    std::vector<std::vector<std::size_t>> cell_faces(nx*ny);
    auto add=[&](std::initializer_list<std::size_t> vv,std::size_t c){
        std::vector<std::size_t> key(vv);
        std::sort(key.begin(),key.end());
        const auto it=faces.find(key);
        if(it!=faces.end()) {
            mesh.ownership().set_neighbour(it->second,static_cast<int>(c));
            return it->second;
        }
        const std::size_t f=mesh.faces().n_faces();
        mesh.faces().push_face(vv);
        mesh.ownership().resize(mesh.faces().n_faces());
        mesh.ownership().set_owner(f,c);
        mesh.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
        faces.emplace(std::move(key),f);
        return f;
    };

    for(std::size_t j=0;j<ny;++j)
        for(std::size_t i=0;i<nx;++i) {
            const std::size_t c=j*nx+i;
            const auto a=id(i,j,0),b=id(i+1,j,0),cc=id(i+1,j+1,0),d=id(i,j+1,0);
            const auto e=id(i,j,1),f=id(i+1,j,1),g=id(i+1,j+1,1),h=id(i,j+1,1);
            cell_faces[c]={
                add({a,d,cc,b},c),add({e,f,g,h},c),add({a,b,f,e},c),
                add({d,h,g,cc},c),add({a,e,h,d},c),add({b,cc,g,f},c)};
        }
    for(const auto& cf:cell_faces) mesh.cells().push_cell(cf);

    Patch inlet{"inlet",PatchType::INLET,{}};
    Patch outlet{"outlet",PatchType::OUTLET,{}};
    Patch bottom{"bottom",PatchType::WALL,{}};
    Patch top{"top",PatchType::WALL,{}};
    Patch front{"front",PatchType::EMPTY,{}};
    Patch back{"back",PatchType::EMPTY,{}};

    for(std::size_t f=0;f<mesh.n_faces();++f) {
        if(mesh.ownership().neighbour(f)>=0) continue;
        const auto& v=mesh.faces().vertices();
        const auto b=v.begin()+static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto e=b+static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));
        double x=0,y=0,z=0;
        for(auto p=b;p!=e;++p){x+=mesh.points().x(*p);y+=mesh.points().y(*p);z+=mesh.points().z(*p);}
        const double n=static_cast<double>(mesh.faces().face_size(f));
        x/=n;y/=n;z/=n;
        constexpr double tol=1e-12;
        if(std::abs(x)<tol) inlet.face_ids.push_back(f);
        else if(std::abs(x-length)<tol) outlet.face_ids.push_back(f);
        else if(std::abs(y)<tol) bottom.face_ids.push_back(f);
        else if(std::abs(y-1.0)<tol) top.face_ids.push_back(f);
        else if(std::abs(z)<tol) front.face_ids.push_back(f);
        else if(std::abs(z-1.0)<tol) back.face_ids.push_back(f);
        else throw std::runtime_error("unclassified channel boundary");
    }
    mesh.boundary().add_patch(inlet);mesh.boundary().add_patch(outlet);
    mesh.boundary().add_patch(bottom);mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(front);mesh.boundary().add_patch(back);
    return mesh;
}

struct Result {
    Field<double,Location::CELL> U;
    IncompressibleSolveResult solve;
    GeometryCache geometry;
};

Result solve_case(const std::string& mesh_path, std::size_t n)
{
    constexpr double G=1.0;
    constexpr double nu=0.1;
    Mesh mesh;
    if(!cfdx::io::read_mesh_hdf5(mesh_path,mesh))
        throw std::runtime_error("failed to read Python-generated channel mesh: "+mesh_path);
    const auto topo=mesh.topo_validate();
    if(!topo.ok) throw std::runtime_error("channel topology invalid: "+mesh_path);

    Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    Field<double,Location::CELL> p(mesh.n_cells(),"p","Pa",1);
    U.fill(0.0); p.fill(0.0);

    VelocityBoundaryConditions ubc;
    ubc["inlet"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ubc["outlet"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ubc["bottom"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["top"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["front"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ubc["back"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};

    ScalarBoundaryConditions pbc;
    for(const char* s:{"inlet","outlet","bottom","top","front","back"})
        pbc[s]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

    IncompressibleSolverControls c;
    c.algorithm=PressureVelocityAlgorithm::SIMPLE;
    c.density=1.0;
    c.kinematic_viscosity=nu;
    c.body_force={G,0.0,0.0};
    c.linear_max_iterations=4000;
    c.linear_tolerance=1e-10;
    c.convergence.max_iterations=6000;
    c.convergence.relative_tolerance=1e-9;
    c.convergence.continuity_tolerance=1e-9;
    c.coupling.alpha_u=0.7;
    c.coupling.alpha_p=0.3;
    c.use_bounded_convection=true;
    c.convection_scheme=ConvectionScheme::UPWIND;
    c.pressure_reference_cell=0;
    c.pressure_reference_value=0.0;

    const auto solve=solve_steady_incompressible(mesh,U,p,ubc,pbc,c);
    if(!solve.converged) throw std::runtime_error("channel solver did not converge");

    GeometryCache geometry;
    compute_geometry_cache(mesh,geometry);
    return {std::move(U),solve,geometry};
}

double profile_error(const Result& r,std::size_t n)
{
    double sum2=0.0,maxerr=0.0;std::size_t count=0;
    const double* u=r.U.component_data(0);
    for(std::size_t cell=0;cell<r.geometry.cell_centres.size();++cell) {
        const auto& x=r.geometry.cell_centres[cell];
        if(x.x<3.0) continue;
        const double exact=(1.0/(2.0*0.1))*x.y*(1.0-x.y);
        const double e=std::abs(u[cell]-exact);
        sum2+=e*e;maxerr=std::max(maxerr,e);++count;
    }
    if(count==0) throw std::runtime_error("no downstream channel cells");
    const double l2=std::sqrt(sum2/static_cast<double>(count));
    std::cout<<"CHANNEL n="<<n<<" iterations="<<r.solve.iterations
             <<" L2="<<l2<<" Linf="<<maxerr
             <<" continuity="<<r.solve.history.back().continuity_linf<<"\n";
    return l2;
}

} // namespace

int main(int argc,char** argv)
{
    try {
        const bool quick=argc==2 && std::string(argv[1])=="--quick";
        if(quick && argc != 3) throw std::invalid_argument("quick channel requires one mesh path");
        if(!quick && argc != 5) throw std::invalid_argument("full channel requires three mesh paths");

        if(argc < 3) throw std::invalid_argument("usage: test_channel_qualification [--quick] mesh_n16.h5 [mesh_n32.h5 mesh_n64.h5]");
        const auto r16=solve_case(argv[2],16);
        const double e16=profile_error(r16,16);
        if(quick) {
            if(e16>0.20 || r16.solve.history.back().continuity_linf>1e-8)
                throw std::runtime_error("channel quick gate failed");
            std::cout<<"CHANNEL_QUICK: PASS\n";
            return 0;
        }

        if(argc < 5) throw std::invalid_argument("full channel qualification requires three Python-generated meshes");
        const auto r32=solve_case(argv[3],32);
        const auto r64=solve_case(argv[4],64);
        const double e32=profile_error(r32,32);
        const double e64=profile_error(r64,64);
        if(!(e64<e32 && e32<e16))
            throw std::runtime_error("channel refinement is not monotone");

        const double p=std::log(e32/e64)/std::log(2.0);
        std::cout<<"CHANNEL observed_order="<<p<<"\n";
        if(!(p>0.5) || r32.solve.history.back().continuity_linf>1e-8 ||
           r64.solve.history.back().continuity_linf>1e-8)
            throw std::runtime_error("channel qualification gate failed");

        std::cout<<"CHANNEL_QUALIFICATION: PASS\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"CHANNEL_QUALIFICATION: FAIL: "<<e.what()<<"\n";
        return 1;
    }
}
