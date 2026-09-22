#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {

Mesh make_channel(std::size_t n, double height, double length)
{
    Mesh m;
    m.points().resize(8 * n);
    const double dy = height / static_cast<double>(n);

    for (std::size_t i = 0; i < n; ++i) {
        const double y0 = dy * static_cast<double>(i);
        const double y1 = dy * static_cast<double>(i + 1);
        const std::size_t b = 8 * i;
        const double p[8][3] = {
            {0.0,y0,0.0}, {length,y0,0.0}, {length,y1,0.0}, {0.0,y1,0.0},
            {0.0,y0,1.0}, {length,y0,1.0}, {length,y1,1.0}, {0.0,y1,1.0}
        };
        for (std::size_t j = 0; j < 8; ++j)
            m.points().set(b + j, p[j][0], p[j][1], p[j][2]);
    }

    std::vector<std::size_t> bottom, top, x0, x1, z0, z1, internal;
    bottom.push_back(m.faces().n_faces());
    m.faces().push_face({0,1,5,4});
    top.push_back(m.faces().n_faces());
    m.faces().push_face({8*(n-1)+3,8*(n-1)+7,8*(n-1)+6,8*(n-1)+2});

    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t b = 8 * i;
        x0.push_back(m.faces().n_faces()); m.faces().push_face({b+0,b+4,b+7,b+3});
        x1.push_back(m.faces().n_faces()); m.faces().push_face({b+1,b+2,b+6,b+5});
        z0.push_back(m.faces().n_faces()); m.faces().push_face({b+0,b+3,b+2,b+1});
        z1.push_back(m.faces().n_faces()); m.faces().push_face({b+4,b+5,b+6,b+7});
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        const std::size_t b = 8 * i;
        internal.push_back(m.faces().n_faces());
        m.faces().push_face({b+3,b+7,b+6,b+2});
    }

    m.ownership().resize(m.n_faces());
    for (std::size_t i = 0; i + 1 < n; ++i) {
        m.ownership().set_owner(internal[i], i);
        m.ownership().set_neighbour(internal[i], static_cast<int>(i + 1));
    }
    for (const auto f : bottom) {
        m.ownership().set_owner(f, 0);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    for (const auto f : top) {
        m.ownership().set_owner(f, n - 1);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    for (std::size_t i = 0; i < n; ++i) {
        m.ownership().set_owner(x0[i], i);
        m.ownership().set_neighbour(x0[i], FaceOwnership::BOUNDARY);
        m.ownership().set_owner(x1[i], i);
        m.ownership().set_neighbour(x1[i], FaceOwnership::BOUNDARY);
        m.ownership().set_owner(z0[i], i);
        m.ownership().set_neighbour(z0[i], FaceOwnership::BOUNDARY);
        m.ownership().set_owner(z1[i], i);
        m.ownership().set_neighbour(z1[i], FaceOwnership::BOUNDARY);
    }
    for (std::size_t i = 0; i < n; ++i) {
        const std::vector<std::size_t> faces = {
            i == 0 ? bottom[0] : internal[i-1],
            i + 1 == n ? top[0] : internal[i],
            x0[i], x1[i], z0[i], z1[i]
        };
        m.cells().push_cell(faces);
    }

    Patch p;
    p.name = "bottom"; p.type = PatchType::WALL; p.face_ids = bottom; m.boundary().add_patch(p);
    p = {}; p.name = "top"; p.type = PatchType::WALL; p.face_ids = top; m.boundary().add_patch(p);
    p = {}; p.name = "x0"; p.type = PatchType::INLET; p.face_ids = x0; m.boundary().add_patch(p);
    p = {}; p.name = "x1"; p.type = PatchType::OUTLET; p.face_ids = x1; m.boundary().add_patch(p);
    p = {}; p.name = "z0"; p.type = PatchType::SYMMETRY; p.face_ids = z0; m.boundary().add_patch(p);
    p = {}; p.name = "z1"; p.type = PatchType::SYMMETRY; p.face_ids = z1; m.boundary().add_patch(p);
    return m;
}

void check_case(std::size_t n)
{
    constexpr double H = 1.0;
    constexpr double L = 1.5;
    constexpr double rho = 1.0;
    constexpr double mu = 1.0;
    constexpr double Uwall = 3.0;
    constexpr double pin = 18.0;
    constexpr double pout = 0.0;

    Mesh mesh = make_channel(n, H, L);
    Field<double,Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    Field<double,Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    U.fill(0.0);
    p.fill(0.0);

    VelocityBoundaryConditions ubc;
    ubc["bottom"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
    ubc["top"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{Uwall,0.0,0.0}};
    ubc["x0"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0.0,0.0,0.0}};
    ubc["x1"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0.0,0.0,0.0}};
    ubc["z0"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0.0,0.0,0.0}};
    ubc["z1"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0.0,0.0,0.0}};

    ScalarBoundaryConditions pbc;
    pbc["bottom"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
    pbc["top"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
    pbc["x0"] = {ScalarBoundaryType::FIXED_VALUE,pin,0.0};
    pbc["x1"] = {ScalarBoundaryType::FIXED_VALUE,pout,0.0};
    pbc["z0"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
    pbc["z1"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

    IncompressibleSolverControls c;
    c.algorithm = PressureVelocityAlgorithm::SIMPLE;
    c.convergence.max_iterations = 500;
    c.convergence.relative_tolerance = 1e-8;
    c.convergence.continuity_tolerance = 1e-10;
    c.linear_max_iterations = 5000;
    c.linear_tolerance = 1e-8;
    c.density = rho;
    c.kinematic_viscosity = mu / rho;
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = pin - (pin-pout)/(2.0*static_cast<double>(n));

    const auto result = solve_steady_incompressible(mesh,U,p,ubc,pbc,c);
    if (!result.converged)
        throw std::runtime_error("VMFL004 coupled solver did not converge");


    const auto geometry = build_fv_geometry(mesh);
    double max_error = 0.0;
    double mean_u = 0.0;
    double volume = 0.0;
    for (std::size_t i = 0; i < mesh.n_cells(); ++i) {
        const double y = geometry.cell_centres[i].y;
        const double exact = Uwall*y/H + (-12.0)/(2.0*mu)*y*(y-H);
        max_error = std::max(max_error, std::abs(U(i,0)-exact));
        mean_u += U(i,0) * geometry.cell_volumes[i];
        volume += geometry.cell_volumes[i];
    }
    mean_u /= volume;

    const double relative_max = max_error / Uwall;
    const double mean_exact = 2.5;
    const double relative_mean = std::abs(mean_u-mean_exact)/mean_exact;
    const auto& last = result.history.back();

    std::cout << "VMFL004 N=" << n
              << " max_abs_error=" << max_error
              << " max_rel_error=" << relative_max
              << " mean_velocity=" << mean_u
              << " mean_rel_error=" << relative_mean
              << " continuity_linf=" << last.continuity_linf
              << " iterations=" << result.iterations << "\n";

    if (relative_max > 5e-2 || relative_mean > 5e-2)
        throw std::runtime_error("VMFL004 velocity profile mismatch");
    if (last.continuity_linf > 1e-8)
        throw std::runtime_error("VMFL004 continuity error too large");
}

} // namespace

int main()
{
    try {
        check_case(32);
        check_case(64);
        std::cout << "VMFL004_VALIDATION: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "VMFL004_VALIDATION: FAIL: " << e.what() << "\n";
        return 1;
    }
}
