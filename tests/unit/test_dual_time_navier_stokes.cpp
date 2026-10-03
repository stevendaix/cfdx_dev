#include "cfdx/physics/dual_time_navier_stokes.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

Mesh make_channel_mesh(std::size_t nx, std::size_t ny)
{
    Mesh mesh;
    const std::size_t plane = (nx + 1) * (ny + 1);
    mesh.points().resize(2 * plane);
    const auto id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };
    for (std::size_t j = 0; j <= ny; ++j)
        for (std::size_t i = 0; i <= nx; ++i) {
            const double x = static_cast<double>(i) / nx;
            const double y = static_cast<double>(j) / ny;
            mesh.points().set(id(i,j,0), x, y, 0.0);
            mesh.points().set(id(i,j,1), x, y, 1.0);
        }

    std::map<std::vector<std::size_t>, std::size_t> face_map;
    std::vector<std::vector<std::size_t>> cell_faces(nx * ny);
    auto add_face = [&](std::initializer_list<std::size_t> vertices, std::size_t cell) {
        std::vector<std::size_t> key(vertices);
        std::sort(key.begin(), key.end());
        auto it = face_map.find(key);
        if (it != face_map.end()) {
            mesh.ownership().set_neighbour(it->second, static_cast<int>(cell));
            return it->second;
        }
        const auto f = mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        mesh.ownership().set_owner(f, cell);
        mesh.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        face_map.emplace(std::move(key), f);
        return f;
    };

    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i) {
            const auto c = j * nx + i;
            const auto a=id(i,j,0), b=id(i+1,j,0), c0=id(i+1,j+1,0), d=id(i,j+1,0);
            const auto e=id(i,j,1), f=id(i+1,j,1), g=id(i+1,j+1,1), h=id(i,j+1,1);
            cell_faces[c] = {
                add_face({a,d,c0,b},c), add_face({e,f,g,h},c),
                add_face({a,b,f,e},c), add_face({d,h,g,c0},c),
                add_face({a,e,h,d},c), add_face({b,c0,g,f},c)};
        }
    for (const auto& faces : cell_faces) mesh.cells().push_cell(faces);

    Patch inlet{"inlet",PatchType::INLET,{}}, outlet{"outlet",PatchType::OUTLET,{}};
    Patch bottom{"bottom",PatchType::WALL,{}}, top{"top",PatchType::WALL,{}};
    Patch front{"front",PatchType::EMPTY,{}}, back{"back",PatchType::EMPTY,{}};
    for (std::size_t f=0; f<mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f)>=0) continue;
        const auto& vs=mesh.faces().vertices();
        auto it=vs.begin()+static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        auto end=it+static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));
        double x=0,y=0,z=0; std::size_t n=0;
        for(;it!=end;++it){x+=mesh.points().x(*it);y+=mesh.points().y(*it);z+=mesh.points().z(*it);++n;}
        x/=n;y/=n;z/=n;
        if(std::abs(x)<1e-12) inlet.face_ids.push_back(f);
        else if(std::abs(x-1)<1e-12) outlet.face_ids.push_back(f);
        else if(std::abs(y)<1e-12) bottom.face_ids.push_back(f);
        else if(std::abs(y-1)<1e-12) top.face_ids.push_back(f);
        else if(std::abs(z)<1e-12) front.face_ids.push_back(f);
        else if(std::abs(z-1)<1e-12) back.face_ids.push_back(f);
        else throw std::runtime_error("unclassified channel face");
    }
    mesh.boundary().add_patch(inlet); mesh.boundary().add_patch(outlet);
    mesh.boundary().add_patch(bottom); mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(front); mesh.boundary().add_patch(back);
    return mesh;
}

void require(bool ok, const std::string& msg) {
    if (!ok) throw std::runtime_error(msg);
}

IncompressibleSolverControls base_controls()
{
    IncompressibleSolverControls c;
    c.algorithm = PressureVelocityAlgorithm::COUPLED;
    c.coupling.alpha_u = 0.7;
    c.coupling.alpha_p = 0.3;
    c.coupling.coupled_max_iterations = 300;
    c.coupling.coupled_linear_tolerance = 1e-9;
    c.convergence.max_iterations = 300;
    c.convergence.relative_tolerance = 1e-7;
    c.convergence.continuity_tolerance = 1e-7;
    c.linear_max_iterations = 300;
    c.linear_tolerance = 1e-9;
    c.density = 1.0;
    c.kinematic_viscosity = 0.1;
    c.use_bounded_convection = true;
    c.convection_scheme = ConvectionScheme::UPWIND;
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = 0.0;
    c.coupled_linear_solver.krylov = KrylovModel::FGMRES;
    c.coupled_linear_solver.preconditioner = PreconditionerModel::CoupledBlockSchur;
    return c;
}

} // namespace

int main()
{
    try {
        const Mesh mesh = make_channel_mesh(4, 4);
        Field<double, Location::CELL> U(mesh.n_cells(),"U","m/s",3);
        Field<double, Location::CELL> p(mesh.n_cells(),"p","Pa",1);
        U.fill(0.0); p.fill(0.0);

        VelocityBoundaryConditions ubc;
        ubc["inlet"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
        ubc["outlet"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
        ubc["bottom"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        ubc["top"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};
        ubc["front"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
        ubc["back"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};

        ScalarBoundaryConditions pbc;
        for (const auto& name : {"inlet","outlet","bottom","top","front","back"})
            pbc[name]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        DualTimeNavierStokesLifecycle lifecycle(U,p);
        lifecycle.configure_problem(mesh,ubc,pbc);

        auto solver = base_controls();
        DualTimeNavierStokesControls tc;
        tc.dt_min=1e-4; tc.dt_max=0.2; tc.dt_shrink=0.5;
        tc.dt_growth=1.5;
        tc.temporal_absolute_tolerance=1e-8;
        tc.temporal_relative_tolerance=1e-5;
        tc.max_retries=12;

        const auto first=lifecycle.advance(0.02,solver,tc);
        require(first.retries==0,"first step must use BE bootstrap without retry");
        require(lifecycle.step()==1 && lifecycle.history_valid(),"first accepted step must commit history");
        const auto checkpoint=lifecycle.checkpoint();

        const auto second=lifecycle.advance(0.08,solver,tc);
        require(second.retries > 0,"tight temporal gate must exercise physical retry");
        require(second.dt_accepted < 0.08,"retry must reduce the accepted physical dt");
        require(second.dt_accepted <= 0.08,"accepted dt must not exceed requested dt");
        require(second.temporal_error <= 1.0,"accepted temporal error gate must pass");
        require(second.dt_proposed >= tc.dt_min && second.dt_proposed <= tc.dt_max,
                "adaptive dt proposal must remain bounded");
        require(lifecycle.step()==2,"accepted step counter must advance exactly once");
        require(lifecycle.physical_time()>checkpoint.physical_time,"physical time must advance only on acceptance");

        DualTimeNavierStokesLifecycle restarted(U,p);
        restarted.configure_problem(mesh,ubc,pbc);
        restarted.restore(checkpoint);
        DualTimeNavierStokesLifecycle reference(U,p);
        reference.configure_problem(mesh,ubc,pbc);
        reference.restore(checkpoint);

        const auto a=restarted.advance(second.dt_accepted,solver,tc);
        const auto b=reference.advance(second.dt_accepted,solver,tc);
        require(std::abs(restarted.physical_time()-reference.physical_time())<1e-14,
                "restart physical time must be deterministic");
        double max_du=0.0, max_dp=0.0;
        for(std::size_t c=0;c<U.size();++c) {
            for(std::size_t d=0;d<3;++d)
                max_du=std::max(max_du,std::abs(restarted.U().component_data(d)[c]-reference.U().component_data(d)[c]));
            max_dp=std::max(max_dp,std::abs(restarted.p(c)-reference.p(c)));
        }
        require(max_du<1e-12 && max_dp<1e-12,"restart continuation must reproduce the accepted solution");
        require(a.high_order.converged && b.high_order.converged,"restart solves must converge");

        std::cout<<"production dual-time Navier-Stokes PASS"
                 <<" retries="<<second.retries
                 <<" dt="<<second.dt_accepted
                 <<" error="<<second.temporal_error
                 <<" max_du="<<max_du
                 <<" max_dp="<<max_dp<<"\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"production dual-time Navier-Stokes FAIL: "<<e.what()<<"\n";
        return 1;
    }
}
