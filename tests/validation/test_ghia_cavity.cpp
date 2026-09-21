#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

struct CavityCase { double reynolds; std::size_t nx; std::size_t ny; std::size_t max_iterations; };
struct Sample { double coordinate; double reference; };
struct Comparison { double rms = 0.0; double max_abs = 0.0; double max_relative_to_lid = 0.0; };
struct CavityResult { Field<double, Location::CELL> velocity; IncompressibleSolveResult solve; FvGeometry geometry; };

Mesh make_cavity_mesh(std::size_t nx, std::size_t ny)
{
    if (nx < 4 || ny < 4) throw std::invalid_argument("cavity mesh is too small");
    Mesh mesh;
    mesh.points().resize((nx + 1) * (ny + 1) * 2);
    const auto point_id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };
    for (std::size_t j=0;j<=ny;++j) for (std::size_t i=0;i<=nx;++i) {
        const double x=static_cast<double>(i)/static_cast<double>(nx);
        const double y=static_cast<double>(j)/static_cast<double>(ny);
        mesh.points().set(point_id(i,j,0),x,y,0.0);
        mesh.points().set(point_id(i,j,1),x,y,1.0);
    }

    std::map<std::vector<std::size_t>,std::size_t> face_map;
    std::vector<std::vector<std::size_t>> cell_faces(nx*ny);
    // Ownership is sized after face creation because shared faces are deduplicated.
    auto add_or_get_face=[&](std::initializer_list<std::size_t> vertices,std::size_t cell) {
        std::vector<std::size_t> sorted(vertices);
        std::sort(sorted.begin(),sorted.end());
        const auto found=face_map.find(sorted);
        if(found!=face_map.end()) {
            mesh.ownership().set_neighbour(found->second,static_cast<int>(cell));
            return found->second;
        }
        const std::size_t face=mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        face_map.emplace(std::move(sorted),face);
        mesh.ownership().set_owner(face,cell);
        mesh.ownership().set_neighbour(face,FaceOwnership::BOUNDARY);
        return face;
    };

    for(std::size_t j=0;j<ny;++j) for(std::size_t i=0;i<nx;++i) {
        const std::size_t cell=j*nx+i;
        const std::size_t a=point_id(i,j,0), b=point_id(i+1,j,0);
        const std::size_t c=point_id(i+1,j+1,0), d=point_id(i,j+1,0);
        const std::size_t e=point_id(i,j,1), f=point_id(i+1,j,1);
        const std::size_t g=point_id(i+1,j+1,1), h=point_id(i,j+1,1);
        cell_faces[cell]={
            add_or_get_face({a,d,c,b},cell), add_or_get_face({e,f,g,h},cell),
            add_or_get_face({a,b,f,e},cell), add_or_get_face({d,h,g,c},cell),
            add_or_get_face({a,e,h,d},cell), add_or_get_face({b,c,g,f},cell)};
    }
    for(const auto& faces:cell_faces) mesh.cells().push_cell(faces);

    Patch bottom; bottom.name="bottom"; bottom.type=PatchType::WALL;
    Patch top; top.name="top"; top.type=PatchType::WALL;
    Patch left; left.name="left"; left.type=PatchType::WALL;
    Patch right; right.name="right"; right.type=PatchType::WALL;
    Patch front; front.name="front"; front.type=PatchType::WALL;
    Patch back; back.name="back"; back.type=PatchType::WALL;

    for(std::size_t f=0;f<mesh.n_faces();++f) {
        if(mesh.ownership().neighbour(f)>=0) continue;
        const auto& all_vertices=mesh.faces().vertices();
        const auto begin=all_vertices.begin()+static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto end=begin+static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));
        double x=0,y=0,z=0;
        for(auto it=begin;it!=end;++it) { x+=mesh.points().x(*it); y+=mesh.points().y(*it); z+=mesh.points().z(*it); }
        const double nverts=static_cast<double>(mesh.faces().face_size(f));
        x/=nverts; y/=nverts; z/=nverts;
        constexpr double tol=1e-12;
        if(std::abs(y)<tol) bottom.face_ids.push_back(f);
        else if(std::abs(y-1)<tol) top.face_ids.push_back(f);
        else if(std::abs(x)<tol) left.face_ids.push_back(f);
        else if(std::abs(x-1)<tol) right.face_ids.push_back(f);
        else if(std::abs(z)<tol) front.face_ids.push_back(f);
        else if(std::abs(z-1)<tol) back.face_ids.push_back(f);
        else throw std::runtime_error("cavity: unclassified boundary face");
    }
    mesh.boundary().add_patch(bottom); mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(left); mesh.boundary().add_patch(right);
    mesh.boundary().add_patch(front); mesh.boundary().add_patch(back);
    return mesh;
}

CavityResult solve_cavity(const CavityCase& test)
{
    Mesh mesh=make_cavity_mesh(test.nx,test.ny);
    Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    Field<double,Location::CELL> p(mesh.n_cells(),"p","Pa",1);
    U.fill(0.0); p.fill(0.0);

    VelocityBoundaryConditions ubc;
    ubc["bottom"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["left"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["right"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["front"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["back"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["top"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};

    ScalarBoundaryConditions pbc;
    for(const char* name:{"bottom","top","left","right","front","back"})
        pbc[name]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

    IncompressibleSolverControls controls;
    controls.algorithm=PressureVelocityAlgorithm::SIMPLE;
    controls.density=1.0;
    controls.kinematic_viscosity=1.0/test.reynolds;
    controls.linear_max_iterations=5000;
    controls.linear_tolerance=1e-10;
    controls.pressure_reference_cell=0;
    controls.pressure_reference_value=0.0;
    controls.use_bounded_convection=true;
    controls.convection_scheme=ConvectionScheme::UPWIND;
    controls.coupling.alpha_u=0.5;
    controls.coupling.alpha_p=0.2;
    controls.convergence.max_iterations=test.max_iterations;
    controls.convergence.relative_tolerance=1e-8;
    controls.convergence.continuity_tolerance=1e-8;

    const auto solve=solve_steady_incompressible(mesh,U,p,ubc,pbc,controls);
    if(!solve.converged)
        throw std::runtime_error("Ghia cavity did not converge for Re="+std::to_string(test.reynolds));
    return {std::move(U),solve,build_fv_geometry(mesh)};
}

double interpolate_line(const Field<double,Location::CELL>& U,
                        const FvGeometry& geometry,
                        std::size_t nx, std::size_t ny,
                        bool horizontal_component,
                        double coordinate)
{
    const double centre=0.5;
    const double eps=1e-12;
    if (coordinate<=eps || coordinate>=1.0-eps)
        return horizontal_component
            ? (coordinate>=1.0-eps ? 1.0 : 0.0)
            : 0.0;

    const double* values=U.component_data(horizontal_component?0:1);
    const double qx=horizontal_component ? centre : coordinate;
    const double qy=horizontal_component ? coordinate : centre;

    const double fx=qx*static_cast<double>(nx)-0.5;
    const double fy=qy*static_cast<double>(ny)-0.5;
    const long ix0=std::clamp(static_cast<long>(std::floor(fx)),0L,static_cast<long>(nx)-2L);
    const long iy0=std::clamp(static_cast<long>(std::floor(fy)),0L,static_cast<long>(ny)-2L);
    const std::size_t i0=static_cast<std::size_t>(ix0);
    const std::size_t j0=static_cast<std::size_t>(iy0);
    const double tx=std::clamp(fx-static_cast<double>(ix0),0.0,1.0);
    const double ty=std::clamp(fy-static_cast<double>(iy0),0.0,1.0);

    const auto at=[&](std::size_t i,std::size_t j) {
        return values[j*nx+i];
    };
    const double v00=at(i0,j0);
    const double v10=at(i0+1,j0);
    const double v01=at(i0,j0+1);
    const double v11=at(i0+1,j0+1);
    return (1.0-ty)*((1.0-tx)*v00+tx*v10)
         + ty*((1.0-tx)*v01+tx*v11);
}

Comparison compare(const Field<double,Location::CELL>& U,const FvGeometry& geometry,
                   std::size_t nx,std::size_t ny,bool horizontal_component,const std::vector<Sample>& samples)
{
    double sum2=0.0,max_abs=0.0;
    for(const auto& sample:samples) {
        const double error=std::abs(interpolate_line(U,geometry,nx,ny,horizontal_component,sample.coordinate)-sample.reference);
        sum2+=error*error; max_abs=std::max(max_abs,error);
    }
    return {std::sqrt(sum2/static_cast<double>(samples.size())),max_abs,max_abs};
}

std::vector<Sample> u_reference(double re)
{
    static const double y[]={0.0000,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5000,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1.0000};
    static const double u100[]={0.00000,-0.03717,-0.04192,-0.04775,-0.06434,-0.10150,-0.15662,-0.21090,-0.20581,-0.13641,0.00332,0.23151,0.68717,0.73722,0.78871,0.84123,1.00000};
    static const double u400[]={0.00000,-0.08186,-0.09266,-0.10338,-0.14612,-0.24299,-0.32726,-0.17119,-0.11477,0.02135,0.16256,0.29093,0.55892,0.61756,0.68439,0.75837,1.00000};
    const double* values=re==100.0?u100:u400;
    std::vector<Sample> result; for(std::size_t i=0;i<17;++i) result.push_back({y[i],values[i]}); return result;
}

std::vector<Sample> v_reference(double re)
{
    static const double x[]={0.0000,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5000,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1.0000};
    static const double v100[]={0.00000,0.09233,0.10091,0.10890,0.12317,0.16077,0.17507,0.17527,0.05454,-0.24533,-0.22445,-0.16914,-0.10313,-0.08864,-0.07391,-0.05906,0.00000};
    static const double v400[]={0.00000,0.18360,0.19713,0.20920,0.22965,0.28124,0.30203,0.30174,0.05186,-0.38598,-0.44993,-0.23827,-0.22847,-0.19254,-0.15663,-0.12146,0.00000};
    const double* values=re==100.0?v100:v400;
    std::vector<Sample> result; for(std::size_t i=0;i<17;++i) result.push_back({x[i],values[i]}); return result;
}

struct CaseMetrics { Comparison u; Comparison v; };

CaseMetrics run_case(const CavityCase& test)
{
    const auto result=solve_cavity(test);
    const auto u=compare(result.velocity,result.geometry,test.nx,test.ny,true,u_reference(test.reynolds));
    const auto v=compare(result.velocity,result.geometry,test.nx,test.ny,false,v_reference(test.reynolds));
    std::cout<<"GHIA Re="<<test.reynolds<<" grid="<<test.nx<<"x"<<test.ny
             <<" iterations="<<result.solve.iterations
             <<" continuity="<<result.solve.history.back().continuity_linf
             <<" continuity_norm="<<result.solve.history.back().continuity_normalized
             <<" momentum_eq="<<result.solve.history.back().momentum_equation_residual
             <<" U_RMS="<<u.rms<<" U_max="<<u.max_abs<<" V_RMS="<<v.rms<<" V_max="<<v.max_abs<<"\n";
    const double max_allowed=test.nx>=64?0.10:0.15;
    const double rms_allowed=0.075;
    if(u.max_abs>max_allowed || v.max_abs>max_allowed ||
       u.rms>rms_allowed || v.rms>rms_allowed)
        throw std::runtime_error("Ghia velocity profile mismatch");
    if(result.solve.history.back().continuity_linf>1e-7 ||
       !std::isfinite(result.solve.history.back().continuity_normalized) ||
       !std::isfinite(result.solve.history.back().momentum_equation_residual) ||
       result.solve.history.back().momentum_equation_residual>1e-7)
        throw std::runtime_error("Ghia physical residual gate failed");
    return {u,v};
}

} // namespace

int main()
{
    try {
        const auto r32 = run_case({100.0,32,32,2500});
        const auto r64 = run_case({100.0,64,64,5000});
        const auto r128 = run_case({100.0,128,128,12000});
        run_case({400.0,64,64,9000});

        const double p_v_max = std::log(r64.v.max_abs / r128.v.max_abs) / std::log(2.0);
        const double p_v_rms = std::log(r64.v.rms / r128.v.rms) / std::log(2.0);
        const double p_u_max = std::log(r64.u.max_abs / r128.u.max_abs) / std::log(2.0);
        const double p_u_rms = std::log(r64.u.rms / r128.u.rms) / std::log(2.0);
        std::cout << "GHIA Re=100 observed_order"
                  << " U_RMS=" << p_u_rms
                  << " U_max=" << p_u_max
                  << " V_RMS=" << p_v_rms
                  << " V_max=" << p_v_max << "\n";
        (void)r32;
        if (!(p_u_rms > 0.50) || !(p_v_rms > 0.50) ||
            !(p_u_max > 0.50) || !(p_v_max > 0.50))
            throw std::runtime_error("Ghia Re=100 mesh convergence is insufficient");
        std::cout<<"GHIA_CAVITY_VALIDATION: PASS\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"GHIA_CAVITY_VALIDATION: FAIL: "<<e.what()<<"\n";
        return 1;
    }
}
