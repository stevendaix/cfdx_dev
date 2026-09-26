#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
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
    const auto topo = mesh.topo_validate();
    if (!topo.ok) {
        throw std::runtime_error(
            "Couette channel mesh topology invalid: " +
            (topo.errors.empty() ? std::string("unknown error") : topo.errors.front()));
    }
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
    c.coupling.n_fractional_steps =
        algorithm == PressureVelocityAlgorithm::FRACTIONAL_STEP ? 2 : 1;
    c.coupling.coupled_max_iterations = 2000;
    c.coupling.coupled_linear_tolerance = 1e-10;
    c.coupling.n_outer_correctors =
        algorithm == PressureVelocityAlgorithm::PIMPLE ? 2 : 1;
    c.convergence.max_iterations = 3000;
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
    c.diagnostics.coupled_matrix_summary = true;
    c.diagnostics.freeze_state_probe = true;
    c.diagnostics.debug_cell = 33;

    const auto diagnostic_geometry = build_fv_geometry(mesh);
    c.iteration_output_callback =
        [algorithm, scheme, bounded, diagnostic_geometry, ubc, pbc](
            std::size_t iter,
            double,
            const Mesh& callback_mesh,
            const Field<double, Location::CELL>& callback_U,
            const Field<double, Location::CELL>& callback_p) {
            const auto grad_p = gauss_gradient_with_boundary(
                callback_p, callback_mesh, diagnostic_geometry, pbc);
            double gradp_linf = 0.0;
            double gradp_l2_sum = 0.0;
            double u_linf = 0.0;
            for (std::size_t cell = 0; cell < callback_mesh.n_cells(); ++cell) {
                const double gx = grad_p.component_data(0)[cell];
                const double gy = grad_p.component_data(1)[cell];
                const double gz = grad_p.component_data(2)[cell];
                const double gp2 = gx*gx + gy*gy + gz*gz;
                gradp_linf = std::max(gradp_linf, std::sqrt(gp2));
                gradp_l2_sum += gp2;
                for (std::size_t d = 0; d < 3; ++d)
                    u_linf = std::max(u_linf,
                        std::abs(callback_U.component_data(d)[cell]));
            }
            auto reconstructed_flux = make_mass_flux(
                callback_mesh, diagnostic_geometry, callback_U, 1.0, ubc);
            double reconstructed_cont_linf = 0.0;
            const auto* faces = callback_mesh.cells().faces_data();
            const auto* offsets = callback_mesh.cells().offsets_data();
            for (std::size_t cell = 0; cell < callback_mesh.n_cells(); ++cell) {
                double div = 0.0;
                for (Offset k = offsets[cell]; k < offsets[cell + 1]; ++k) {
                    const std::size_t face = faces[k];
                    div += callback_mesh.ownership().owner(face) == cell
                        ? reconstructed_flux(face) : -reconstructed_flux(face);
                }
                reconstructed_cont_linf =
                    std::max(reconstructed_cont_linf, std::abs(div));
            }
            std::cout << "TRACE iter=" << iter
                      << " algorithm=" << static_cast<int>(algorithm)
                      << " scheme=" << static_cast<int>(scheme)
                      << " bounded=" << (bounded ? "true" : "false")
                      << " gradp_linf=" << gradp_linf
                      << " gradp_l2=" << std::sqrt(
                             gradp_l2_sum / std::max<std::size_t>(
                                 callback_mesh.n_cells(), 1))
                      << " U_linf=" << u_linf
                      << " reconstructed_continuity=" << reconstructed_cont_linf
                      << "\n";
            return true;
        };

    const auto solve = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
    return {std::move(U), std::move(p), solve, build_fv_geometry(mesh)};
}

struct ProfileError {
    double l2 = 0.0;
    double linf = 0.0;
};

ProfileError profile_error(const RunResult& r, std::size_t nx, std::size_t ny)
{
    double l2 = 0.0;
    double linf = 0.0;
    std::size_t count = 0;

    for (std::size_t c = 0; c < r.U.size(); ++c) {
        const std::size_t i = c % nx;
        const std::size_t j = c / nx;
        if (i >= nx || j >= ny) continue;
        const double y = (static_cast<double>(j) + 0.5) / static_cast<double>(ny);
        const double error = r.U.component_data(0)[c] - y;
        l2 += error * error;
        linf = std::max(linf, std::abs(error));
        ++count;
    }

    if (count == 0) throw std::runtime_error("empty Couette profile");
    return {std::sqrt(l2 / static_cast<double>(count)), linf};
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
    c.pressure_linear_solver.null_space = NullSpaceModel::Constant;

    const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
    if (!result.converged)
        throw std::runtime_error("pure-Neumann gauge case did not converge");
    assert_close(p(0), 0.0, 1e-12, "pressure reference");
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const bool quick = argc == 2 && std::string(argv[1]) == "--quick";
        if (argc > 1 && !quick)
            throw std::invalid_argument("usage: test_phase9_acceptance [--quick]");

        std::cout << (quick
            ? "PHASE9: Couette SIMPLE/PISO/COUPLED smoke verification\n"
            : "PHASE9: Couette full pressure-velocity-system verification\n");

        struct Case {
            const char* name;
            PressureVelocityAlgorithm algorithm;
            ConvectionScheme scheme;
            bool bounded;
        };

        // The smoke set spans segregated, multi-corrector and monolithic
        // coupling. The full campaign exercises every exposed algorithm.
        // Physical gates remain identical in both modes.
        std::vector<Case> algorithm_cases = {
            {"SIMPLE/upwind/bounded", PressureVelocityAlgorithm::SIMPLE,
             ConvectionScheme::UPWIND, true},
            {"PISO/upwind/bounded", PressureVelocityAlgorithm::PISO,
             ConvectionScheme::UPWIND, true},
            {"COUPLED/upwind/bounded", PressureVelocityAlgorithm::COUPLED,
             ConvectionScheme::UPWIND, true},
        };
        if (!quick) {
            algorithm_cases.insert(algorithm_cases.end(), {
                {"SIMPLEC/upwind/bounded", PressureVelocityAlgorithm::SIMPLEC,
                 ConvectionScheme::UPWIND, true},
                {"PIMPLE/upwind/bounded", PressureVelocityAlgorithm::PIMPLE,
                 ConvectionScheme::UPWIND, true},
                {"FRACTIONAL_STEP/upwind/bounded", PressureVelocityAlgorithm::FRACTIONAL_STEP,
                 ConvectionScheme::UPWIND, true},
            });
        }

        // Couette is an affine velocity field.  On the cell-centred 16-cell
        // mesh the exact discrete maximum is (16 - 0.5) / 16 = 0.96875.
        // These gates must therefore be tight enough to detect a biased
        // discretisation or an incompletely converged segregated solve.
        // Algorithm invariance is a physical-equivalence gate, not merely a
        // smoke test: the coupling algorithm must agree with SIMPLE to 1e-5
        // in the cell-centred velocity field.
        constexpr double profile_l2_tolerance = 1.0e-6;
        constexpr double profile_linf_tolerance = 1.0e-6;
        constexpr double couette_umax_exact = (16.0 - 0.5) / 16.0;
        constexpr double transverse_velocity_tolerance = 1.0e-7;
        constexpr double pressure_uniformity_tolerance = 1.0e-7;
        constexpr double boundary_velocity_tolerance = 1.0e-8;

        std::vector<RunResult> results;
        std::vector<std::string> successful_models;
        std::vector<std::string> failed_models;
        results.reserve(algorithm_cases.size());

        // Diagnostic history intentionally records several distinct operators:
        //   continuity              -> authoritative conservative flux balance
        //   continuity_norm         -> scale-normalized conservation measure
        //   corrected_flux_continuity -> post-pressure-correction flux balance
        //   reconstructed_velocity_continuity -> flux reconstructed from cell U
        //   flux_velocity_mismatch  -> diagnostic only (Rhie-Chow is pressure-dependent)
        //   mom_internal/boundary   -> localisation of a possible residual floor
        // These are evidence channels, not interchangeable acceptance metrics.
        // In particular, Rhie-Chow flux and arithmetic interpolation of U are
        // different operators by construction, so their mismatch is not required
        // to vanish even when the conservative continuity equation converges.

        auto print_history = [](const char* name, const RunResult& result) {
            std::cout << "MODEL_BEGIN " << name << "\n";
            std::cout << "ITERATION_HISTORY_BEGIN " << name << "\n";
            for (const auto& ih : result.solve.history) {
                std::cout << "ITER " << ih.iteration
                          << " continuity=" << ih.continuity_linf
                          << " continuity_norm=" << ih.continuity_normalized
                          << " momentum=" << ih.momentum_residual
                          << " momentum_eq_rel=" << ih.momentum_equation_residual_relative
                          << " pressure=" << ih.pressure_residual
                          << " dU=" << ih.velocity_change_inf
                          << " dp=" << ih.pressure_change_inf
                          << " corrected_flux_continuity=" << ih.corrected_flux_continuity_linf
                          << " reconstructed_velocity_continuity=" << ih.reconstructed_velocity_continuity_linf
                          << " flux_velocity_mismatch=" << ih.flux_velocity_mismatch_linf
                          << " mom_x=" << ih.momentum_equation_residual_components[0]
                          << " mom_y=" << ih.momentum_equation_residual_components[1]
                          << " mom_z=" << ih.momentum_equation_residual_components[2]
                          << " mom_internal=" << ih.momentum_equation_residual_internal
                          << " mom_boundary=" << ih.momentum_equation_residual_boundary
                          << " mom_cell=" << ih.momentum_residual_cell
                          << " mom_patch=" << ih.momentum_residual_patch
                          << " mom_no_pressure=" << ih.momentum_residual_no_pressure
                          << " mom_pressure=" << ih.momentum_pressure_contribution
                          << " gradp_linf=" << ih.pressure_gradient_linf
                          << " gradp_l2=" << ih.pressure_gradient_l2
                          << " momentum_linear_iterations=" << ih.momentum_linear_iterations
                          << " pressure_linear_iterations=" << ih.pressure_linear_iterations
                          << "\n";
            }
            std::cout << "ITERATION_HISTORY_END " << name << "\n";
        };

        for (const auto& test : algorithm_cases) {
            std::cout << "MODEL_CONFIG algorithm=" << test.name
                      << " nx=8 ny=16 bounded=" << (test.bounded ? "true" : "false")
                      << " alpha_u=0.7 alpha_p=0.3"
                      << " pressure_correctors="
                      << (test.algorithm == PressureVelocityAlgorithm::PISO ||
                          test.algorithm == PressureVelocityAlgorithm::PIMPLE ? 2 : 1)
                      << " fractional_steps="
                      << (test.algorithm == PressureVelocityAlgorithm::FRACTIONAL_STEP ? 2 : 1)
                      << "\n";

            try {
                auto result = run_couette_channel(
                    test.algorithm, test.scheme, test.bounded, 8, 16);

                print_history(test.name, result);

                const auto error = profile_error(result, 8, 16);
                double max_abs_uy = 0.0;
                double max_abs_uz = 0.0;
                double max_u = -std::numeric_limits<double>::infinity();
                double min_u = std::numeric_limits<double>::infinity();
                double p_min = std::numeric_limits<double>::infinity();
                double p_max = -std::numeric_limits<double>::infinity();

                for (std::size_t c = 0; c < result.U.size(); ++c) {
                    const double ux = result.U.component_data(0)[c];
                    const double uy = result.U.component_data(1)[c];
                    const double uz = result.U.component_data(2)[c];
                    const double pc = result.p(c);
                    if (!std::isfinite(ux) || !std::isfinite(uy) ||
                        !std::isfinite(uz) || !std::isfinite(pc))
                        throw std::runtime_error("non-finite solution field");
                    max_abs_uy = std::max(max_abs_uy, std::abs(uy));
                    max_abs_uz = std::max(max_abs_uz, std::abs(uz));
                    max_u = std::max(max_u, ux);
                    min_u = std::min(min_u, ux);
                    p_min = std::min(p_min, pc);
                    p_max = std::max(p_max, pc);
                }

                std::vector<std::string> gates;
                if (!result.solve.converged)
                    gates.push_back("solver_not_converged");
                if (result.solve.history.empty())
                    gates.push_back("empty_iteration_history");
                else {
                    const auto& h = result.solve.history.back();
                    if (!(h.continuity_linf < 1e-7))
                        gates.push_back("continuity_linf");
                    if (!(h.continuity_normalized < 1e-7))
                        gates.push_back("continuity_normalized");
                    if (!(h.momentum_equation_residual_relative < 1e-7))
                        gates.push_back("momentum_equation_residual_relative");
                    if (!(h.corrected_flux_continuity_linf < 1e-7))
                        gates.push_back("corrected_flux_continuity");
                    if (test.algorithm == PressureVelocityAlgorithm::PIMPLE &&
                        result.solve.iterations < 2)
                        gates.push_back("pimple_outer_correctors");
                }
                if (!(error.l2 < profile_l2_tolerance && error.linf < profile_linf_tolerance))
                    gates.push_back("analytical_profile");
                if (!(max_abs_uy < transverse_velocity_tolerance))
                    gates.push_back("Uy");
                if (!(max_abs_uz < transverse_velocity_tolerance))
                    gates.push_back("Uz");
                if (!(std::isfinite(max_u) &&
                      std::abs(max_u - couette_umax_exact) < profile_linf_tolerance))
                    gates.push_back("Umax");
                if (!(min_u > -boundary_velocity_tolerance))
                    gates.push_back("Umin");
                if (!std::isfinite(p_min) || !std::isfinite(p_max) ||
                    p_max - p_min > pressure_uniformity_tolerance)
                    gates.push_back("pressure_uniformity");

                const auto& pressure_context =
                    result.solve.pressure_linear_context;
                if (test.algorithm != PressureVelocityAlgorithm::COUPLED) {
                    if (pressure_context.full_setups != 1)
                        gates.push_back("pressure_context_full_setup");
                    if (pressure_context.solves <= 1)
                        gates.push_back("pressure_context_not_reused");
                    if (pressure_context.full_setups +
                            pressure_context.numeric_updates +
                            pressure_context.unchanged_reuses !=
                        pressure_context.solves)
                        gates.push_back("pressure_context_accounting");
                }

                const auto& h = result.solve.history.empty()
                    ? IncompressibleIteration{}
                    : result.solve.history.back();

                std::cout << "MODEL_RESULT " << test.name
                          << " solver_converged=" << (result.solve.converged ? "true" : "false")
                          << " iterations=" << result.solve.iterations
                          << " history_size=" << result.solve.history.size()
                          << " profile_L2=" << error.l2
                          << " profile_Linf=" << error.linf
                          << " Umax=" << max_u
                          << " Umin=" << min_u
                          << " |Uy|max=" << max_abs_uy
                          << " |Uz|max=" << max_abs_uz
                          << " dp_range=" << (p_max - p_min)
                          << " continuity=" << h.continuity_linf
                          << " continuity_norm=" << h.continuity_normalized
                          << " momentum_eq_rel=" << h.momentum_equation_residual_relative
                          << " corrected_flux_continuity=" << h.corrected_flux_continuity_linf
                          << " reconstructed_velocity_continuity=" << h.reconstructed_velocity_continuity_linf
                          << " flux_velocity_mismatch=" << h.flux_velocity_mismatch_linf
                          << " mom_x=" << h.momentum_equation_residual_components[0]
                          << " mom_y=" << h.momentum_equation_residual_components[1]
                          << " mom_z=" << h.momentum_equation_residual_components[2]
                          << " mom_internal=" << h.momentum_equation_residual_internal
                          << " mom_boundary=" << h.momentum_equation_residual_boundary
                          << " mom_cell=" << h.momentum_residual_cell
                          << " mom_patch=" << h.momentum_residual_patch
                          << " mom_no_pressure=" << h.momentum_residual_no_pressure
                          << " mom_pressure=" << h.momentum_pressure_contribution
                          << " gradp_linf=" << h.pressure_gradient_linf
                          << " gradp_l2=" << h.pressure_gradient_l2
                          << " pressure_full_setups="
                          << pressure_context.full_setups
                          << " pressure_numeric_updates="
                          << pressure_context.numeric_updates
                          << " pressure_unchanged_reuses="
                          << pressure_context.unchanged_reuses
                          << " pressure_context_solves="
                          << pressure_context.solves
                          << " gates_failed=" << gates.size()
                          << "\n";

                if (gates.empty()) {
                    successful_models.push_back(test.name);
                    results.push_back(std::move(result));
                } else {
                    failed_models.push_back(test.name);
                    std::cout << "MODEL_FAILURES " << test.name;
                    for (const auto& gate : gates) std::cout << " " << gate;
                    std::cout << "\n";
                }
            } catch (const std::exception& e) {
                failed_models.push_back(test.name);
                std::cout << "MODEL_RESULT " << test.name
                          << " execution_exception=" << e.what() << "\n";
                std::cout << "MODEL_FAILURES " << test.name
                          << " execution_exception\n";
            }
        }

        std::cout << "MODEL_SUMMARY successful=" << successful_models.size()
                  << " failed=" << failed_models.size() << "\n";
        if (!failed_models.empty()) {
            std::cout << "FAILED_MODELS";
            for (const auto& name : failed_models) std::cout << " " << name;
            std::cout << "\n";
        }

        // Algorithm invariance is evaluated only across models that actually
        // produced a valid result. A failure in one model must not prevent the
        // remaining models from running and exposing their diagnostics.
        if (!results.empty()) {
            for (std::size_t k = 1; k < results.size(); ++k) {
                double max_du = 0.0;
                for (std::size_t c = 0; c < results[k].U.size(); ++c)
                    max_du = std::max(
                        max_du,
                        std::abs(results[k].U.component_data(0)[c] -
                                 results.front().U.component_data(0)[c]));
                std::cout << "ALGORITHM_INVARIANCE model=" << successful_models[k]
                          << " vs=" << successful_models.front()
                          << " max_abs_dU=" << max_du << "\n";
                if (!(max_du < 1.0e-5))
                    failed_models.push_back(
                        successful_models[k] + ":algorithm_invariance");
            }
        }

        if (!quick) {
            // Exercise the convection/flux assembly independently of the
            // pressure-velocity algorithm. These cases use the same SIMPLE
            // pressure correction but cover the alternate convection controls.
            const auto second_order = run_couette_channel(
                PressureVelocityAlgorithm::SIMPLE,
                ConvectionScheme::SECOND_ORDER_UPWIND, true);
            const auto unbounded = run_couette_channel(
                PressureVelocityAlgorithm::SIMPLE,
                ConvectionScheme::UPWIND, false);

            for (const auto& pair : {
                     std::pair<const char*, const RunResult*>{
                         "SIMPLE/SOU/bounded", &second_order},
                     {"SIMPLE/upwind/unbounded", &unbounded}}) {
                const auto error = profile_error(*pair.second, 8, 16);
                if (!(error.l2 < profile_l2_tolerance &&
                      error.linf < profile_linf_tolerance))
                    throw std::runtime_error(
                        std::string(pair.first) + ": convection gate failed");
                const auto& h = pair.second->solve.history.back();
                if (!(h.continuity_linf < 1e-7) ||
                    !(h.momentum_equation_residual_relative < 1e-7))
                    throw std::runtime_error(
                        std::string(pair.first) + ": conservation gate failed");
                std::cout << pair.first
                          << ": profile L2/Linf=" << error.l2 << "/" << error.linf
                          << " continuity=" << h.continuity_linf << "\n";
            }

            // A pure-Neumann pressure field has a gauge freedom. Starting from
            // a non-zero uniform pressure must converge to the same physical
            // state without introducing a local pressure jump.
            run_pure_neumann_gauge();
        }

        if (!failed_models.empty()) {
            std::cout << "PHASE9_ACCEPTANCE: FAIL\n";
            return 1;
        }
        std::cout << (quick ? "COUETTE_QUICK: PASS\n" : "PHASE9_ACCEPTANCE: PASS\n");
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "PHASE9_ACCEPTANCE: FAIL: " << e.what() << "\n";
        return 1;
    }
}
