#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

Mesh make_channel_mesh(std::size_t nx, std::size_t ny)
{
    if (nx < 2 || ny < 4) throw std::invalid_argument("channel mesh too small");

    Mesh mesh;
    const std::size_t plane = (nx + 1) * (ny + 1);
    mesh.points().resize(2 * plane);
    const auto id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };

    for (std::size_t j = 0; j <= ny; ++j) {
        for (std::size_t i = 0; i <= nx; ++i) {
            const double x = static_cast<double>(i) / static_cast<double>(nx);
            const double y = static_cast<double>(j) / static_cast<double>(ny);
            mesh.points().set(id(i, j, 0), x, y, 0.0);
            mesh.points().set(id(i, j, 1), x, y, 1.0);
        }
    }

    std::map<std::vector<std::size_t>, std::size_t> face_map;
    std::vector<std::vector<std::size_t>> cell_faces(nx * ny);

    auto add_face = [&](std::initializer_list<std::size_t> vertices,
                        std::size_t cell) {
        std::vector<std::size_t> key(vertices);
        std::sort(key.begin(), key.end());
        const auto found = face_map.find(key);
        if (found != face_map.end()) {
            mesh.ownership().set_neighbour(found->second, static_cast<int>(cell));
            return found->second;
        }

        const std::size_t f = mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        mesh.ownership().set_owner(f, cell);
        mesh.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        face_map.emplace(std::move(key), f);
        return f;
    };

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t c = j * nx + i;
            const auto a = id(i, j, 0);
            const auto b = id(i + 1, j, 0);
            const auto c0 = id(i + 1, j + 1, 0);
            const auto d = id(i, j + 1, 0);
            const auto e = id(i, j, 1);
            const auto f = id(i + 1, j, 1);
            const auto g = id(i + 1, j + 1, 1);
            const auto h = id(i, j + 1, 1);
            cell_faces[c] = {
                add_face({a, d, c0, b}, c),
                add_face({e, f, g, h}, c),
                add_face({a, b, f, e}, c),
                add_face({d, h, g, c0}, c),
                add_face({a, e, h, d}, c),
                add_face({b, c0, g, f}, c)};
        }
    }

    for (const auto& faces : cell_faces) mesh.cells().push_cell(faces);

    Patch inlet{"inlet", PatchType::INLET, {}};
    Patch outlet{"outlet", PatchType::OUTLET, {}};
    Patch bottom{"bottom", PatchType::WALL, {}};
    Patch top{"top", PatchType::WALL, {}};
    Patch front{"front", PatchType::EMPTY, {}};
    Patch back{"back", PatchType::EMPTY, {}};

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) >= 0) continue;

        const auto& vertices = mesh.faces().vertices();
        const auto begin = vertices.begin() +
            static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto end = begin +
            static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));

        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        for (auto it = begin; it != end; ++it) {
            x += mesh.points().x(*it);
            y += mesh.points().y(*it);
            z += mesh.points().z(*it);
        }
        const double n = static_cast<double>(mesh.faces().face_size(f));
        x /= n; y /= n; z /= n;

        constexpr double tol = 1e-12;
        if (std::abs(x) < tol) inlet.face_ids.push_back(f);
        else if (std::abs(x - 1.0) < tol) outlet.face_ids.push_back(f);
        else if (std::abs(y) < tol) bottom.face_ids.push_back(f);
        else if (std::abs(y - 1.0) < tol) top.face_ids.push_back(f);
        else if (std::abs(z) < tol) front.face_ids.push_back(f);
        else if (std::abs(z - 1.0) < tol) back.face_ids.push_back(f);
        else throw std::runtime_error("channel boundary face is unclassified");
    }

    mesh.boundary().add_patch(inlet);
    mesh.boundary().add_patch(outlet);
    mesh.boundary().add_patch(bottom);
    mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(front);
    mesh.boundary().add_patch(back);
    return mesh;
}

struct RunResult {
    Field<double, Location::CELL> U;
    Field<double, Location::CELL> p;
    IncompressibleSolveResult solve;
    FvGeometry geometry;
};

RunResult run_couette_channel(
    PressureVelocityAlgorithm algorithm,
    ConvectionScheme scheme,
    bool bounded,
    std::size_t nx = 8,
    std::size_t ny = 16)
{
    Mesh mesh = make_channel_mesh(nx, ny);
    Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    U.fill(0.0);
    p.fill(0.5);

    VelocityBoundaryConditions ubc;
    ubc["inlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    ubc["outlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    ubc["bottom"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0, 0, 0}};
    ubc["top"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {1, 0, 0}};
    ubc["front"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    ubc["back"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};

    ScalarBoundaryConditions pbc;
    pbc["inlet"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    pbc["outlet"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    pbc["bottom"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    pbc["top"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    pbc["front"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    pbc["back"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

    IncompressibleSolverControls c;
    c.algorithm = algorithm;
    c.coupling.alpha_u = 0.7;
    c.coupling.alpha_p = 0.3;
    c.coupling.n_pressure_correctors =
        algorithm == PressureVelocityAlgorithm::PISO ||
        algorithm == PressureVelocityAlgorithm::PIMPLE ? 2 : 1;
    c.coupling.n_outer_correctors =
        algorithm == PressureVelocityAlgorithm::PIMPLE ? 2 : 1;
    c.convergence.max_iterations = 250;
    c.convergence.relative_tolerance = 1e-8;
    c.convergence.continuity_tolerance = 1e-8;
    c.linear_max_iterations = 2000;
    c.linear_tolerance = 1e-10;
    c.density = 1.0;
    c.kinematic_viscosity = 0.1;
    c.body_force = {0.0, 0.0, 0.0};
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = 0.0;
    c.use_bounded_convection = bounded;
    c.convection_scheme = scheme;

    const auto solve = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
    if (!solve.converged || solve.history.empty())
        throw std::runtime_error("Couette channel did not converge");

    const auto& h = solve.history.back();
    if (!(h.continuity_linf < 1e-7) ||
        !(h.continuity_normalized < 1e-7) ||
        !(h.momentum_equation_residual_relative < 1e-7))
        throw std::runtime_error("physical convergence gate failed");

    return {std::move(U), std::move(p), solve, build_fv_geometry(mesh)};
}

double profile_error(const RunResult& r, std::size_t nx, std::size_t ny)
{
    double l2 = 0.0;
    double linf = 0.0;
    std::size_t count = 0;

    for (std::size_t c = 0; c < r.U.size(); ++c) {
        const std::size_t i = c % nx;
        const std::size_t j = c / nx;
        if (i >= nx || j >= ny) continue;
        const double y = (static_cast<double>(j) + 0.5) / static_cast<double>(ny);
        const double exact = y;
        const double error = r.U.component_data(0)[c] - exact;
        l2 += error * error;
        linf = std::max(linf, std::abs(error));
        ++count;
    }

    return std::sqrt(l2 / static_cast<double>(count)) + linf;
}

void assert_close(double value, double reference, double tolerance, const char* label)
{
    if (!std::isfinite(value) || std::abs(value - reference) > tolerance)
        throw std::runtime_error(
            std::string(label) + " mismatch: value=" + std::to_string(value) +
            " reference=" + std::to_string(reference));
}

void run_pure_neumann_gauge()
{
    Mesh mesh = make_channel_mesh(4, 8);
    Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    U.fill(0.0);
    p.fill(37.0);

    VelocityBoundaryConditions ubc;
    for (const char* n : {"inlet", "outlet", "front", "back"})
        ubc[n] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    ubc["bottom"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0, 0, 0}};
    ubc["top"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0, 0, 0}};

    ScalarBoundaryConditions pbc;
    for (const char* n : {"inlet", "outlet", "bottom", "top", "front", "back"})
        pbc[n] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

    IncompressibleSolverControls c;
    c.algorithm = PressureVelocityAlgorithm::SIMPLE;
    c.convergence.max_iterations = 5;
    c.convergence.relative_tolerance = 1e-10;
    c.convergence.continuity_tolerance = 1e-10;
    c.linear_max_iterations = 100;
    c.linear_tolerance = 1e-12;
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = 0.0;

    const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
    if (!result.converged)
        throw std::runtime_error("pure-Neumann gauge case did not converge");
    assert_close(p(0), 0.0, 1e-12, "pressure reference");
}

} // namespace

int main()
{
    try {
        std::cout << "PHASE9: Couette analytical verification\n";

        const auto simple = run_couette_channel(
            PressureVelocityAlgorithm::SIMPLE, ConvectionScheme::UPWIND, true);
        const auto simplec = run_couette_channel(
            PressureVelocityAlgorithm::SIMPLEC, ConvectionScheme::UPWIND, true);
        const auto piso = run_couette_channel(
            PressureVelocityAlgorithm::PISO, ConvectionScheme::UPWIND, true);
        const auto pimple = run_couette_channel(
            PressureVelocityAlgorithm::PIMPLE, ConvectionScheme::UPWIND, true);

        const double e_simple = profile_error(simple, 8, 16);
        const double e_simplec = profile_error(simplec, 8, 16);
        const double e_piso = profile_error(piso, 8, 16);
        const double e_pimple = profile_error(pimple, 8, 16);
        if (pimple.solve.iterations < 2)
            throw std::runtime_error("PIMPLE n_outer_correctors was not honored");

        std::cout << "Couette profile error: SIMPLE=" << e_simple
                  << " SIMPLEC=" << e_simplec
                  << " PISO=" << e_piso
                  << " PIMPLE=" << e_pimple << "\n";

        if (!(e_simple < 0.20 && e_simplec < 0.20 &&
              e_piso < 0.20 && e_pimple < 0.20))
            throw std::runtime_error("Couette profile validation failed");

        const auto second_order = run_couette_channel(
            PressureVelocityAlgorithm::SIMPLE,
            ConvectionScheme::SECOND_ORDER_UPWIND, true);
        const auto unbounded = run_couette_channel(
            PressureVelocityAlgorithm::SIMPLE,
            ConvectionScheme::UPWIND, false);

        const double e_second = profile_error(second_order, 8, 16);
        const double e_unbounded = profile_error(unbounded, 8, 16);
        std::cout << "Convection/boundedness error: SOU=" << e_second
                  << " unbounded-upwind=" << e_unbounded << "\n";
        if (!(e_second < 0.20 && e_unbounded < 0.20))
            throw std::runtime_error("convection/bounded-deferred-correction gate failed");

        run_pure_neumann_gauge();

        // The analytic profile is u(y)=f_x/(2 nu)*y*(1-y).
        // With f_x=1 and nu=0.1, u_max=1.25.
        double max_u = 0.0;
        for (std::size_t c = 0; c < simple.U.size(); ++c)
            max_u = std::max(max_u, simple.U.component_data(0)[c]);
        if (!(max_u > 0.9 && max_u < 1.5))
            throw std::runtime_error("Couette peak velocity is outside the physical gate");

        std::cout << "PHASE9_ACCEPTANCE: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "PHASE9_ACCEPTANCE: FAIL: " << e.what() << "\n";
        return 1;
    }
}
