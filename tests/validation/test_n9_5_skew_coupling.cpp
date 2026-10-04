#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

Mesh make_skew_channel_mesh(std::size_t nx, std::size_t ny)
{
    if (nx < 4 || ny < 4)
        throw std::invalid_argument("skew channel mesh too small");

    Mesh mesh;
    constexpr double pi = 3.141592653589793238462643383279502884;
    const std::size_t plane = (nx + 1) * (ny + 1);
    mesh.points().resize(2 * plane);

    const auto id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };

    // Preserve all boundary vertices exactly. The interior displacement is a
    // smooth manufactured mesh deformation, so the physical patch locations
    // and boundary conditions remain those of the planar Couette problem. The
    // deformation is deliberately moderate: this qualification targets
    // pressure-velocity equivalence on a skew/non-orthogonal mesh, not a
    // high-distortion mesh-robustness limit (covered separately by N10).
    for (std::size_t j = 0; j <= ny; ++j) {
        for (std::size_t i = 0; i <= nx; ++i) {
            const double x0 = static_cast<double>(i) / static_cast<double>(nx);
            const double y0 = static_cast<double>(j) / static_cast<double>(ny);
            const bool boundary = i == 0 || i == nx || j == 0 || j == ny;
            const double sx = boundary
                ? 0.0
                : 0.03 * std::sin(pi * x0) * std::sin(pi * y0);
            const double sy = boundary
                ? 0.0
                : 0.015 * std::sin(2.0 * pi * x0) * std::sin(pi * y0);
            for (std::size_t k = 0; k < 2; ++k)
                mesh.points().set(id(i, j, k), x0 + sx, y0 + sy, static_cast<double>(k));
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
            const std::size_t cell = j * nx + i;
            const auto a = id(i, j, 0);
            const auto b = id(i + 1, j, 0);
            const auto c = id(i + 1, j + 1, 0);
            const auto d = id(i, j + 1, 0);
            const auto e = id(i, j, 1);
            const auto f = id(i + 1, j, 1);
            const auto g = id(i + 1, j + 1, 1);
            const auto h = id(i, j + 1, 1);
            cell_faces[cell] = {
                add_face({a, d, c, b}, cell),
                add_face({e, f, g, h}, cell),
                add_face({a, b, f, e}, cell),
                add_face({d, h, g, c}, cell),
                add_face({a, e, h, d}, cell),
                add_face({b, c, g, f}, cell)};
        }
    }

    for (const auto& faces : cell_faces)
        mesh.cells().push_cell(faces);

    Patch inlet{"inlet", PatchType::INLET, {}};
    Patch outlet{"outlet", PatchType::OUTLET, {}};
    Patch bottom{"bottom", PatchType::WALL, {}};
    Patch top{"top", PatchType::WALL, {}};
    Patch front{"front", PatchType::EMPTY, {}};
    Patch back{"back", PatchType::EMPTY, {}};

    for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
        if (mesh.ownership().neighbour(face) >= 0)
            continue;

        const auto& vertices = mesh.faces().vertices();
        const auto begin = vertices.begin() +
            static_cast<std::ptrdiff_t>(mesh.faces().face_offset(face));
        const auto end = begin +
            static_cast<std::ptrdiff_t>(mesh.faces().face_size(face));

        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        for (auto it = begin; it != end; ++it) {
            x += mesh.points().x(*it);
            y += mesh.points().y(*it);
            z += mesh.points().z(*it);
        }
        const double n = static_cast<double>(mesh.faces().face_size(face));
        x /= n;
        y /= n;
        z /= n;

        constexpr double tol = 1e-12;
        if (std::abs(x) < tol)
            inlet.face_ids.push_back(face);
        else if (std::abs(x - 1.0) < tol)
            outlet.face_ids.push_back(face);
        else if (std::abs(y) < tol)
            bottom.face_ids.push_back(face);
        else if (std::abs(y - 1.0) < tol)
            top.face_ids.push_back(face);
        else if (std::abs(z) < tol)
            front.face_ids.push_back(face);
        else if (std::abs(z - 1.0) < tol)
            back.face_ids.push_back(face);
        else
            throw std::runtime_error("skew channel boundary face is unclassified");
    }

    mesh.boundary().add_patch(inlet);
    mesh.boundary().add_patch(outlet);
    mesh.boundary().add_patch(bottom);
    mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(front);
    mesh.boundary().add_patch(back);
    return mesh;
}

struct Run {
    Field<double, Location::CELL> U;
    Field<double, Location::CELL> p;
    IncompressibleSolveResult solve;
};

Run run_case(PressureVelocityAlgorithm algorithm,
             KrylovModel krylov = KrylovModel::Auto,
             PreconditionerModel preconditioner = PreconditionerModel::Auto)
{
    Mesh mesh = make_skew_channel_mesh(8, 16);
    const auto topo = mesh.topo_validate();
    if (!topo.ok)
        throw std::runtime_error(
            "skew channel topology invalid: " +
            (topo.errors.empty() ? std::string("unknown") : topo.errors.front()));

    Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    U.fill(0.0);
    p.fill(0.0);

    VelocityBoundaryConditions ubc;
    ubc["inlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    ubc["outlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    ubc["bottom"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0, 0, 0}};
    ubc["top"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {1, 0, 0}};
    ubc["front"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    ubc["back"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};

    ScalarBoundaryConditions pbc;
    for (const char* name : {"inlet", "outlet", "bottom", "top", "front", "back"})
        pbc[name] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

    IncompressibleSolverControls c;
    c.algorithm = algorithm;
    c.coupling.alpha_u = 0.7;
    c.coupling.alpha_p = 0.3;
    c.coupling.n_pressure_correctors =
        algorithm == PressureVelocityAlgorithm::PISO ? 2 : 1;
    c.coupling.coupled_max_iterations = 2000;
    c.coupling.coupled_linear_tolerance = 1e-12;
    c.convergence.max_iterations = 3000;
    c.convergence.relative_tolerance = 1e-10;
    c.convergence.continuity_tolerance = 1e-10;
    c.linear_max_iterations = 2000;
    c.linear_tolerance = 1e-12;
    c.density = 1.0;
    c.kinematic_viscosity = 0.1;
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = 0.0;
    c.use_bounded_convection = true;
    c.convection_scheme = ConvectionScheme::UPWIND;
    c.coupled_linear_solver.krylov = krylov;
    c.coupled_linear_solver.preconditioner = preconditioner;

    c.iteration_output_callback =
        [](std::size_t iter, double, const Mesh&,
           const Field<double, Location::CELL>&,
           const Field<double, Location::CELL>&) {
            return true;
        };

    const auto result = solve_steady_incompressible(
        mesh, U, p, ubc, pbc, c);
    return {std::move(U), std::move(p), result};
}

void print_iteration_diagnostics(const char* name, const Run& run)
{
    std::cout << name;
    if (run.solve.history.empty()) {
        std::cout << " history=empty\\n";
        return;
    }
    const auto& h = run.solve.history.back();
    std::cout << " iterations=" << run.solve.iterations
              << " converged=" << run.solve.converged
              << " status=" << static_cast<int>(run.solve.convergence_status)
              << " reason=\\\"" << run.solve.convergence_reason << "\\\""
              << " continuity=" << h.continuity_linf
              << " continuity_norm=" << h.continuity_normalized
              << " corrected_flux_continuity=" << h.corrected_flux_continuity_linf
              << " reconstructed_velocity_continuity=" << h.reconstructed_velocity_continuity_linf
              << " flux_velocity_mismatch=" << h.flux_velocity_mismatch_linf
              << " mass_local_linf=" << h.mass_local_linf
              << " mass_normalized=" << h.mass_normalized_imbalance
              << " momentum_eq_rel=" << h.momentum_equation_residual_relative
              << " momentum_no_pressure=" << h.momentum_residual_no_pressure
              << " momentum_pressure=" << h.momentum_pressure_contribution
              << " pressure_grad_linf=" << h.pressure_gradient_linf
              << " velocity_change=" << h.velocity_change_inf
              << " pressure_change=" << h.pressure_change_inf
              << "\\n";
}

void require_converged(const char* name, const Run& run)
{
    if (!run.solve.converged || run.solve.history.empty()) {
        std::ostringstream os;
        os << name << ": solver did not converge"
           << " iterations=" << run.solve.iterations
           << " status=" << static_cast<int>(run.solve.convergence_status)
           << " reason=" << run.solve.convergence_reason;
        if (!run.solve.history.empty()) {
            const auto& h = run.solve.history.back();
            os << " continuity=" << h.continuity_linf
               << " continuity_norm=" << h.continuity_normalized
               << " momentum_eq_rel=" << h.momentum_equation_residual_relative
               << " velocity_change=" << h.velocity_change_inf
               << " pressure_change=" << h.pressure_change_inf;
        }
        throw std::runtime_error(os.str());
    }

    const auto& h = run.solve.history.back();
    if (!(h.continuity_linf < 1e-7))
        throw std::runtime_error(std::string(name) + ": continuity gate failed");
    if (!(h.continuity_normalized < 1e-7))
        throw std::runtime_error(std::string(name) + ": normalized continuity gate failed");
    if (!(h.momentum_equation_residual_relative < 1e-7))
        throw std::runtime_error(std::string(name) + ": momentum true-residual gate failed");

    for (std::size_t c = 0; c < run.U.size(); ++c) {
        for (std::size_t component = 0; component < 3; ++component) {
            if (!std::isfinite(run.U.component_data(component)[c]))
                throw std::runtime_error(std::string(name) + ": non-finite velocity");
        }
        if (!std::isfinite(run.p(c)))
            throw std::runtime_error(std::string(name) + ": non-finite pressure");
    }
}

} // namespace

int main()
{
    try {
        std::cout << "N9.5_SKEW_COUPLING: 8x16 interior-deformed Couette\n";

        const Run simple = run_case(PressureVelocityAlgorithm::SIMPLE);
        const Run block = run_case(
            PressureVelocityAlgorithm::COUPLED,
            KrylovModel::FGMRES,
            PreconditionerModel::CoupledBlockSchur);
        const Run mgr = run_case(
            PressureVelocityAlgorithm::COUPLED,
            KrylovModel::FGMRES,
            PreconditionerModel::MGR);

        // Always print the terminal diagnostics before applying acceptance gates.
        // Keep all gates unchanged: the output is diagnostic evidence, not a
        // relaxation of the N9.5 qualification contract.
        print_iteration_diagnostics("SIMPLE", simple);
        print_iteration_diagnostics("COUPLED/BlockSchur", block);
        print_iteration_diagnostics("COUPLED/MGR", mgr);

        require_converged("SIMPLE", simple);
        require_converged("COUPLED/BlockSchur", block);
        require_converged("COUPLED/MGR", mgr);

        if (!block.solve.coupled_linear_plan_resolved ||
            block.solve.coupled_linear_plan.krylov != KrylovModel::FGMRES ||
            block.solve.coupled_linear_plan.preconditioner !=
                PreconditionerModel::CoupledBlockSchur)
            throw std::runtime_error("BlockSchur request was not the resolved production path");

        if (!mgr.solve.coupled_linear_plan_resolved ||
            mgr.solve.coupled_linear_plan.krylov != KrylovModel::FGMRES ||
            mgr.solve.coupled_linear_plan.preconditioner != PreconditionerModel::MGR)
            throw std::runtime_error("MGR request was not the resolved production path");

        if (simple.solve.authoritative_mass_flux.size() !=
            block.solve.authoritative_mass_flux.size() ||
            simple.solve.authoritative_mass_flux.size() !=
            mgr.solve.authoritative_mass_flux.size())
            throw std::runtime_error("authoritative flux size mismatch");

        double max_simple_block_flux = 0.0;
        double max_simple_mgr_flux = 0.0;
        double scale = 1.0;
        double max_du_block = 0.0;
        double max_du_mgr = 0.0;

        for (std::size_t face = 0;
             face < simple.solve.authoritative_mass_flux.size(); ++face) {
            const double a = simple.solve.authoritative_mass_flux(face);
            const double b = block.solve.authoritative_mass_flux(face);
            const double d = mgr.solve.authoritative_mass_flux(face);
            if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(d))
                throw std::runtime_error("non-finite authoritative mass flux");

            max_simple_block_flux = std::max(max_simple_block_flux, std::abs(b - a));
            max_simple_mgr_flux = std::max(max_simple_mgr_flux, std::abs(d - a));
            scale = std::max({scale, std::abs(a), std::abs(b), std::abs(d)});
        }

        for (std::size_t cell = 0; cell < simple.U.size(); ++cell) {
            for (std::size_t component = 0; component < 3; ++component) {
                max_du_block = std::max(
                    max_du_block,
                    std::abs(block.U.component_data(component)[cell] -
                             simple.U.component_data(component)[cell]));
                max_du_mgr = std::max(
                    max_du_mgr,
                    std::abs(mgr.U.component_data(component)[cell] -
                             simple.U.component_data(component)[cell]));
            }
        }

        const double block_flux_equivalence = max_simple_block_flux / scale;
        const double mgr_flux_equivalence = max_simple_mgr_flux / scale;

        const auto& hs = simple.solve.history.back();
        const auto& hb = block.solve.history.back();
        const auto& hm = mgr.solve.history.back();

        std::cout << "SIMPLE iterations=" << simple.solve.iterations
                  << " continuity=" << hs.continuity_linf
                  << " momentum_eq_rel=" << hs.momentum_equation_residual_relative << "\n";
        std::cout << "BLOCK iterations=" << block.solve.iterations
                  << " continuity=" << hb.continuity_linf
                  << " momentum_eq_rel=" << hb.momentum_equation_residual_relative
                  << " flux_equivalence=" << block_flux_equivalence
                  << " max_abs_dU=" << max_du_block << "\n";
        std::cout << "MGR iterations=" << mgr.solve.iterations
                  << " continuity=" << hm.continuity_linf
                  << " momentum_eq_rel=" << hm.momentum_equation_residual_relative
                  << " flux_equivalence=" << mgr_flux_equivalence
                  << " max_abs_dU=" << max_du_mgr << "\n";

        constexpr double flux_gate = 1e-8;
        constexpr double velocity_gate = 1e-5;
        if (!(block_flux_equivalence < flux_gate))
            throw std::runtime_error("BlockSchur/SIMPLE skew flux equivalence failed");
        if (!(mgr_flux_equivalence < flux_gate))
            throw std::runtime_error("MGR/SIMPLE skew flux equivalence failed");
        if (!(max_du_block < velocity_gate))
            throw std::runtime_error("BlockSchur/SIMPLE skew velocity equivalence failed");
        if (!(max_du_mgr < velocity_gate))
            throw std::runtime_error("MGR/SIMPLE skew velocity equivalence failed");

        std::cout << "N9.5_SKEW_COUPLING: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "N9.5_SKEW_COUPLING: FAIL: " << e.what() << "\n";
        return 1;
    }
}
