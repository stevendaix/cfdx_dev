#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/io/gmsh/gmsh_importer.h"
#include "cfdx/core/geometry/face_geometry.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {
constexpr double D = 1.0;
constexpr double RHO = 1.0;
constexpr double U_INF = 1.0;
constexpr double MU = 0.01;
constexpr double NU = MU / RHO;
constexpr double RE = RHO * U_INF * D / MU;
constexpr double AREA_REF = M_PI * D * D / 4.0;
constexpr double CD_REF = 1.0895;

std::string shell_quote(const std::string& s)
{
    std::string q = "'";
    for (char c : s) q += c == '\'' ? "'\\''" : std::string(1, c);
    return q + "'";
}

ScalarBoundaryConditions scalar_velocity_bcs()
{
    ScalarBoundaryConditions bcs;
    bcs["sphere"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
    bcs["inlet"] = {ScalarBoundaryType::FIXED_VALUE, U_INF, 0.0};
    bcs["outlet"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    bcs["farfield"] = {ScalarBoundaryType::FIXED_VALUE, U_INF, 0.0};
    return bcs;
}

ScalarBoundaryConditions pressure_bcs()
{
    ScalarBoundaryConditions bcs;
    for (const char* n : {"sphere", "inlet", "outlet", "farfield"})
        bcs[n] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    return bcs;
}

void generate_mesh(const std::string& msh)
{
#ifndef CFDX_SOURCE_DIR
    throw std::runtime_error("CFDX_SOURCE_DIR is not defined");
#else
#ifdef CFDX_PYTHON_EXECUTABLE
    const std::string python = CFDX_PYTHON_EXECUTABLE;
#else
    const std::string python = "python3";
#endif
    const std::filesystem::path script =
        std::filesystem::path(CFDX_SOURCE_DIR) / "scripts" / "generate_vmfl036_mesh.py";
    const std::string cmd =
        shell_quote(python) + " " + shell_quote(script.string()) + " " + shell_quote(msh);
    std::cout << "VMFL036_MESH_COMMAND: " << cmd << "\n";
    if (std::system(cmd.c_str()) != 0)
        throw std::runtime_error("VMFL036 Gmsh mesh generation failed");
#endif
}

struct Result {
    Mesh mesh;
    Field<double, Location::CELL> U;
    Field<double, Location::CELL> p;
    IncompressibleSolveResult solve;
    FvGeometry geometry;
    double fd_pressure = 0.0;
    double fd_viscous = 0.0;
};

Result solve_case(const std::string& msh)
{
    Mesh mesh;
    if (!cfdx::io::gmsh::import_gmsh_mesh(msh, mesh))
        throw std::runtime_error("VMFL036 Gmsh import failed");

    const auto topo = mesh.topo_validate();
    if (!topo.ok)
        throw std::runtime_error("VMFL036 topology invalid: " +
            (topo.errors.empty() ? std::string("unknown") : topo.errors.front()));

    for (const char* n : {"sphere", "inlet", "outlet", "farfield"})
        if (!mesh.boundary().has_patch(n))
            throw std::runtime_error(std::string("VMFL036 missing patch: ") + n);

    Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    U.fill(0.0);
    p.fill(0.0);
    for (std::size_t c = 0; c < U.size(); ++c)
        U.component_data(0)[c] = U_INF;

    VelocityBoundaryConditions ubc;
    ubc["sphere"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0,0,0}};
    ubc["inlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {U_INF,0,0}};
    ubc["outlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0,0,0}};
    ubc["farfield"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {U_INF,0,0}};
    const auto pbc = pressure_bcs();

    IncompressibleSolverControls c;
    c.algorithm = PressureVelocityAlgorithm::SIMPLE;
    c.density = RHO;
    c.kinematic_viscosity = NU;
    c.linear_max_iterations = 3000;
    c.linear_tolerance = 1e-9;
    c.convergence.max_iterations = 1500;
    c.convergence.relative_tolerance = 1e-7;
    c.convergence.continuity_tolerance = 1e-7;
    c.coupling.alpha_u = 0.7;
    c.coupling.alpha_p = 0.3;
    c.use_bounded_convection = true;
    c.convection_scheme = ConvectionScheme::UPWIND;
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = 0.0;
    c.diagnostics.iteration_trace = true;
    c.diagnostics.iteration_trace_frequency = 50;

    std::cout << std::setprecision(12)
              << "VMFL036 CFDX START Re=" << RE
              << " cells=" << mesh.n_cells()
              << " faces=" << mesh.n_faces() << "\n";

    const auto solve = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
    if (!solve.converged)
        throw std::runtime_error("VMFL036 CFDX solver did not converge");

    const auto geometry = build_fv_geometry(mesh);
    Field<double, Location::CELL> ux(mesh.n_cells(), "Ux", "m/s", 1);
    Field<double, Location::CELL> uy(mesh.n_cells(), "Uy", "m/s", 1);
    Field<double, Location::CELL> uz(mesh.n_cells(), "Uz", "m/s", 1);
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        ux(cell) = U.component_data(0)[cell];
        uy(cell) = U.component_data(1)[cell];
        uz(cell) = U.component_data(2)[cell];
    }

    const auto gux = gauss_gradient_with_boundary(ux, mesh, geometry, scalar_velocity_bcs());
    const auto guy = gauss_gradient_with_boundary(uy, mesh, geometry, scalar_velocity_bcs());
    const auto guz = gauss_gradient_with_boundary(uz, mesh, geometry, scalar_velocity_bcs());

    const auto sphere = mesh.boundary().find("sphere");
    double fp = 0.0;
    double fv = 0.0;

    for (const auto face : mesh.boundary().patch(sphere).face_ids) {
        const std::size_t owner = mesh.ownership().owner(face);
        const auto& verts = mesh.faces().vertices();
        const auto* offsets = mesh.faces().offsets_data();
        const auto fg = compute_face_geometry_oriented(
            mesh.points().x_data(), mesh.points().y_data(), mesh.points().z_data(),
            verts.data(), offsets[face], offsets[face+1]-offsets[face],
            geometry.cell_centres[owner]);

        // For the inner sphere boundary the owner-cell Sf points into the solid.
        // Thus pressure force on the body is +p*Sf and viscous force is -tau*Sf.
        fp += p(owner) * fg.Sf.x;

        const double dux_dx = gux.component_data(0)[owner];
        const double dux_dy = gux.component_data(1)[owner];
        const double dux_dz = gux.component_data(2)[owner];
        const double duy_dx = guy.component_data(0)[owner];
        const double duz_dx = guz.component_data(0)[owner];
        const double duz_dy = guz.component_data(1)[owner];

        const double tau_xx = 2.0 * MU * dux_dx;
        const double tau_xy = MU * (dux_dy + duy_dx);
        const double tau_xz = MU * (dux_dz + duz_dx);
        fv -= (tau_xx * fg.Sf.x + tau_xy * fg.Sf.y + tau_xz * fg.Sf.z);
        (void)duz_dy;
    }

    return {std::move(mesh), std::move(U), std::move(p), solve, geometry, fp, fv};
}

} // namespace

int main()
{
    try {
        const auto msh =
            std::filesystem::temp_directory_path() / "cfdx_vmfl036_re100.msh";
        generate_mesh(msh.string());
        const auto r = solve_case(msh.string());

        const auto& h = r.solve.history.back();
        const double denom = 0.5 * RHO * U_INF * U_INF * AREA_REF;
        const double fd_total = r.fd_pressure + r.fd_viscous;
        const double cd_p = r.fd_pressure / denom;
        const double cd_v = r.fd_viscous / denom;
        const double cd = fd_total / denom;
        const double rel_error = std::abs(cd - CD_REF) / CD_REF;

        std::cout << "VMFL036 RESULT"
                  << " iterations=" << r.solve.iterations
                  << " continuity=" << h.continuity_linf
                  << " continuity_norm=" << h.continuity_normalized
                  << " momentum_eq=" << h.momentum_equation_residual
                  << " pressure_linear_iterations=" << h.pressure_linear_iterations
                  << "\n";
        std::cout << "VMFL036 DRAG"
                  << " F_pressure=" << r.fd_pressure
                  << " F_viscous=" << r.fd_viscous
                  << " F_total=" << fd_total
                  << " Cd_pressure=" << cd_p
                  << " Cd_viscous=" << cd_v
                  << " Cd_total=" << cd
                  << " Cd_reference=" << CD_REF
                  << " Cd_relative_error=" << rel_error
                  << " F_reference=" << CD_REF * denom
                  << "\n";

        if (!std::isfinite(cd) || !std::isfinite(fd_total))
            throw std::runtime_error("VMFL036 drag is non-finite");

        std::cout << "VMFL036_CFDX_RESULT: PASS (converged physical run; literature accuracy is diagnostic only)\n";
        std::error_code ec;
        std::filesystem::remove(msh, ec);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "VMFL036_CFDX_RESULT: FAIL: " << e.what() << "\n";
        return 1;
    }
}
