#include "cfdx/physics/finite_volume_transport.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

Mesh make_channel(std::size_t n)
{
    Mesh m;
    m.points().resize(8*n);
    for (std::size_t i=0;i<n;++i) {
        const double x0=static_cast<double>(i)/n, x1=static_cast<double>(i+1)/n;
        const std::size_t b=8*i;
        const double p[8][3]={{x0,0,0},{x1,0,0},{x1,1,0},{x0,1,0},
                              {x0,0,1},{x1,0,1},{x1,1,1},{x0,1,1}};
        for (std::size_t q=0;q<8;++q) m.points().set(b+q,p[q][0],p[q][1],p[q][2]);
    }
    std::vector<std::size_t> left(n),right(n),internal,walls;
    std::vector<std::vector<std::size_t>> cf(n);
    auto face=[&](std::initializer_list<std::size_t> v) {
        const auto id=m.n_faces(); m.faces().push_face(std::vector<FaceIndex>(v)); return id;
    };
    for (std::size_t i=0;i<n;++i) {
        const std::size_t b=8*i;
        left[i]=face({b,b+4,b+7,b+3});
        right[i]=face({b+1,b+2,b+6,b+5});
        walls.push_back(face({b,b+1,b+5,b+4}));
        walls.push_back(face({b+3,b+7,b+6,b+2}));
        walls.push_back(face({b,b+3,b+2,b+1}));
        walls.push_back(face({b+4,b+5,b+6,b+7}));
    }
    // Re-use the existing right/left face IDs by identifying the shared face.
    // The cell connectivity is therefore a genuine face-based FV topology.
    for (std::size_t i=0;i+1<n;++i) {
        const std::size_t f=right[i];
        left[i+1]=f;
        internal.push_back(f);
    }
    m.ownership().resize(m.n_faces());
    for (std::size_t i=0;i<n;++i) {
        cf[i]={left[i],right[i],6*i+2,6*i+3,6*i+4,6*i+5};
        for (const auto f:cf[i]) m.ownership().set_owner(f,i);
    }
    for (std::size_t i=0;i+1<n;++i) {
        m.ownership().set_owner(internal[i],i);
        m.ownership().set_neighbour(internal[i],static_cast<std::int64_t>(i+1));
    }
    m.ownership().set_neighbour(left[0],FaceOwnership::BOUNDARY);
    m.ownership().set_neighbour(right[n-1],FaceOwnership::BOUNDARY);
    for (const auto f:walls) m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    for (const auto& faces:cf) m.cells().push_cell(faces);

    Patch pin; pin.name="inlet"; pin.type=PatchType::WALL; pin.face_ids={left[0]};
    Patch pout; pout.name="outlet"; pout.type=PatchType::WALL; pout.face_ids={right[n-1]};
    Patch pw; pw.name="walls"; pw.type=PatchType::WALL; pw.face_ids=walls;
    m.boundary().add_patch(pin); m.boundary().add_patch(pout); m.boundary().add_patch(pw);
    return m;
}

struct Result { double min_value, max_value, outlet_flux; };

Result solve_case(std::size_t n, ConvectionScheme scheme)
{
    auto mesh=make_channel(n);
    const auto geometry=build_fv_geometry(mesh);
    Field<double,Location::FACE> flux(mesh.n_faces(),"F","m3/s",1);
    flux.fill(0.0);
    for (std::size_t f=0;f<mesh.n_faces();++f)
        flux(f)=geometry.face_area_vectors[f].x;

    Field<double,Location::CELL> su(mesh.n_cells(),"Su","1",1);
    Field<double,Location::CELL> sp(mesh.n_cells(),"Sp","1",1);
    Field<double,Location::CELL> phi(mesh.n_cells(),"phi","1",1);
    phi.fill(0.0);

    ScalarBoundaryConditions bc;
    bc["inlet"]={ScalarBoundaryType::FIXED_VALUE,1.0,0.0};
    bc["outlet"]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
    bc["walls"]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

    Vector sol(mesh.n_cells(),0.0);
    for (std::size_t iter=0;iter<30;++iter) {
        for (std::size_t c=0;c<mesh.n_cells();++c) phi(c)=sol(c);
        const auto eq=assemble_scalar_equation(
            mesh,geometry,flux,1.0e-4,su,sp,bc,true,nullptr,nullptr,nullptr,nullptr,
            scheme,&phi);
        const auto r=solve_scalar_equation(eq,sol,{2000,1e-12,1.0});
        if (r.status!=SolverStatus::CONVERGED)
            throw std::runtime_error("steep-gradient transport linear solve did not converge");
    }

    Result out{sol(0),sol(0),0.0};
    for (std::size_t c=0;c<mesh.n_cells();++c) {
        out.min_value=std::min(out.min_value,sol(c));
        out.max_value=std::max(out.max_value,sol(c));
    }
    out.outlet_flux=flux(2)*sol(n-1);
    return out;
}

} // namespace

int main()
{
    try {
        for (const auto& [scheme,name]:{
            std::pair{ConvectionScheme::UPWIND,"upwind"},
            std::pair{ConvectionScheme::SECOND_ORDER_UPWIND,"sou"},
            std::pair{ConvectionScheme::TVD,"tvd"},
            std::pair{ConvectionScheme::QUICK_BOUNDED,"quick_bounded"},
            std::pair{ConvectionScheme::CENTRAL,"central"},
            std::pair{ConvectionScheme::BLENDED,"blended"},
            std::pair{ConvectionScheme::QUICK,"quick"}}) {
            const auto r=solve_case(64,scheme);
            std::cout<<"N4_SOLVER_STEEP scheme="<<name
                     <<" min="<<r.min_value<<" max="<<r.max_value
                     <<" outlet_flux="<<r.outlet_flux<<"\n";
            if (!std::isfinite(r.min_value)||!std::isfinite(r.max_value)||
                !std::isfinite(r.outlet_flux))
                throw std::runtime_error(std::string(name)+" produced non-finite solver-level result");
            if (scheme==ConvectionScheme::UPWIND ||
                scheme==ConvectionScheme::SECOND_ORDER_UPWIND ||
                scheme==ConvectionScheme::TVD ||
                scheme==ConvectionScheme::QUICK_BOUNDED) {
                if (r.min_value < -1e-10 || r.max_value > 1.0+1e-10)
                    throw std::runtime_error(std::string(name)+" violates solver-level boundedness");
            }
        }
        std::cout<<"CONVECTION_SOLVER_QUALIFICATION: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr<<"CONVECTION_SOLVER_QUALIFICATION: FAIL: "<<e.what()<<"\n";
        return 1;
    }
}
