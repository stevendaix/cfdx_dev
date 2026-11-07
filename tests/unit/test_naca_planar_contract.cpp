#include "cfdx/physics/naca_planar_contract.h"
#include "cfdx/physics/forces.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cfdx::core;
using namespace cfdx::physics;

static Mesh one_cell_extrusion() {
    Mesh mesh;
    mesh.points().resize(8);
    for (std::size_t i=0; i<8; ++i) {
        const double x = (i%4 == 1 || i%4 == 2) ? 1.0 : 0.0;
        const double y = (i%4 >= 2) ? 1.0 : 0.0;
        const double z = i >= 4 ? 2.0 : 0.0;
        mesh.points().set(i, x, y, z);
    }
    for (auto verts : {std::initializer_list<std::size_t>{0,3,2,1},
                       {4,5,6,7}, {0,1,5,4}, {1,2,6,5},
                       {2,3,7,6}, {3,0,4,7}})
        mesh.faces().push_face(verts);
    mesh.ownership().resize(6);
    for (std::size_t f=0; f<6; ++f) {
        mesh.ownership().set_owner(f, 0);
        mesh.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0,1,2,3,4,5});
    mesh.boundary().add_patch(Patch{"front",PatchType::EMPTY,{0},{}});
    mesh.boundary().add_patch(Patch{"back",PatchType::EMPTY,{1},{}});
    mesh.boundary().add_patch(Patch{"side",PatchType::WALL,{2,3,4,5},{}});
    assert(mesh.topo_validate().ok);
    return mesh;
}

template <typename F> static void reject(F test) {
    bool threw=false;
    try { test(); }
    catch (const std::exception&) { threw=true; }
    if (!threw) throw std::runtime_error("expected NACA2D fail-closed rejection");
}

static double coefficient_for_span(double span) {
    // A dimensional surface integral scales with the fictitious extrusion.
    // The coefficient must not (for c, rho, U held fixed).
    forces::Reference ref;
    ref.rho_inf=1.0;
    ref.u_inf=1.0;
    ref.area=1.5 * span;
    ref.length=1.5;
    ref.moment_center={0.375,0.0,0.0};
    const forces::SurfaceSample sample{{0.0,0.2,0.0}, {-1.0,0.0,0.0},
                                      span,1.0,{0.2,0.0,0.0}};
    const auto result=forces::integrate({sample},ref);
    if (!(result.cd_pressure > 0) || !(result.cd_viscous > 0) ||
        std::abs(result.cd - (result.cd_pressure+result.cd_viscous)) > 1e-14)
        throw std::runtime_error("pressure/shear split is inconsistent");
    return result.cd;
}

int main() {
    auto mesh=one_cell_extrusion();
    const auto geo=require_naca_planar_geometry(mesh);
    if (geo.span!=2.0 || geo.cells!=1) return 1;
    Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
    U.fill(0); phi.fill(0);
    const auto mass0 = audit_naca_face_mass_budget(mesh, phi, 2.0);
    if (mass0.cell_linf != 0 || mass0.global_net != 0 ||
        mass0.faces_checked != mesh.n_faces()) return 4;
    phi(2) = 0.3;
    const auto mass1 = audit_naca_face_mass_budget(mesh, phi, 2.0);
    if (std::abs(mass1.cell_linf-0.3) > 1e-14 ||
        std::abs(mass1.global_net-0.3) > 1e-14 ||
        std::abs(mass1.local_linf_relative-0.15)>1e-14) return 5;
    phi(2)=0;
    phi(2)=std::numeric_limits<double>::quiet_NaN();
    reject([&]{audit_naca_face_mass_budget(mesh, phi, 2.0);});
    phi(2)=0;
    const auto state=require_naca_planar_state(mesh,U,phi,1,2);
    if (state.checked_empty_faces!=2 || state.max_abs_uz!=0 ||
        state.max_abs_empty_face_flux!=0) return 2;
    U.component_data(2)[0]=1e-3;
    reject([&]{require_naca_planar_state(mesh,U,phi,1,2);});
    U.component_data(2)[0]=0;
    phi(0)=1e-5;
    reject([&]{require_naca_planar_state(mesh,U,phi,1,2);});
    phi(0)=0;
    mesh.boundary().patch(0).type=PatchType::WALL;
    reject([&]{require_naca_planar_geometry(mesh);});
    mesh.boundary().patch(0).type=PatchType::EMPTY;
    mesh.points().set(0,0,0,1);
    reject([&]{require_naca_planar_geometry(mesh);});
    mesh.points().set(0,0,0,0);
    if (std::abs(coefficient_for_span(2)-coefficient_for_span(4))>1e-14)
        return 3;
    std::cout << "NACA_PLANAR_CONTRACT: PASS geometry, Uz, EMPTY flux, "
                 "span-invariant pressure/viscous coefficients\n";
}
