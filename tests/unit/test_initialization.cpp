#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/io/restart/dat_restart.h"
#include "common/test_harness.h"

#include <cmath>
#include <filesystem>
#include <string>

using namespace cfdx::physics;
using namespace cfdx::testing;
using namespace cfdx::core;

namespace {

Mesh make_single_cell_cube()
{
    Mesh m;
    m.points().resize(8);
    const double coords[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i = 0; i < 8; ++i)
        m.points().set(i, coords[i][0], coords[i][1], coords[i][2]);
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for (std::size_t f = 0; f < 6; ++f) {
        m.ownership().set_owner(f, 0);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,1,2,3,4,5};
    m.boundary().add_patch(wall);
    return m;
}

VelocityBoundaryConditions solid_walls()
{
    VelocityBoundaryConditions ubc;
    ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
    return ubc;
}

ScalarBoundaryConditions zero_pressure_gradient()
{
    ScalarBoundaryConditions pbc;
    pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    return pbc;
}

IncompressibleSolverControls one_iteration_controls()
{
    IncompressibleSolverControls controls;
    controls.convergence.max_iterations = 1;
    controls.convergence.continuity_tolerance = 1e-12;
    controls.linear_tolerance = 1e-12;
    return controls;
}

}  // namespace

int main()
{
    // The initial state is part of the numerical method: it sets the reference
    // residual that every later relative residual is measured against. Each
    // strategy is therefore declared and applied explicitly.

    run_case("initialization_modes_parse_and_reject_unknown", [] {
        InitializationMode mode{};
        EXPECT_TRUE(parse_initialization_mode("provided", mode));
        EXPECT_TRUE(mode == InitializationMode::Provided);
        EXPECT_TRUE(parse_initialization_mode("uniform", mode));
        EXPECT_TRUE(mode == InitializationMode::Uniform);
        EXPECT_TRUE(parse_initialization_mode("restart", mode));
        EXPECT_TRUE(mode == InitializationMode::Restart);

        // An unrecognised name must be rejected rather than silently defaulting.
        InitializationMode sentinel = InitializationMode::Restart;
        EXPECT_FALSE(parse_initialization_mode("automatic", sentinel));
        EXPECT_FALSE(parse_initialization_mode("", sentinel));
        EXPECT_FALSE(parse_initialization_mode("UNIFORM", sentinel));
        EXPECT_TRUE(sentinel == InitializationMode::Restart);

        EXPECT_TRUE(std::string(to_string(InitializationMode::Provided)) == "provided");
        EXPECT_TRUE(std::string(to_string(InitializationMode::Uniform)) == "uniform");
        EXPECT_TRUE(std::string(to_string(InitializationMode::Restart)) == "restart");
    });

    run_case("initialization_controls_validation", [] {
        InitializationControls provided;
        validate_initialization_controls(provided);

        InitializationControls uniform;
        uniform.mode = InitializationMode::Uniform;
        uniform.uniform_velocity = Vec3{1.0, -2.0, 3.0};
        uniform.uniform_pressure = 12.5;
        validate_initialization_controls(uniform);

        InitializationControls restart;
        restart.mode = InitializationMode::Restart;
        restart.restart_path = "checkpoint.dat";
        validate_initialization_controls(restart);

        // A non-finite uniform state is a configuration error, not a value to
        // propagate into the first residual.
        InitializationControls nonfinite = uniform;
        nonfinite.uniform_pressure = std::nan("");
        EXPECT_THROW(validate_initialization_controls(nonfinite), std::invalid_argument);
        nonfinite = uniform;
        nonfinite.uniform_velocity.x = std::nan("");
        EXPECT_THROW(validate_initialization_controls(nonfinite), std::invalid_argument);
    });

    run_case("uniform_initialization_overwrites_caller_state", [] {
        const Mesh m = make_single_cell_cube();
        auto ubc = solid_walls();
        auto pbc = zero_pressure_gradient();

        // Start from a deliberately wrong state: if Uniform did not overwrite
        // it, the solve would start from the caller's values.
        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p(1, "p", "Pa", 1);
        U.set(0, 9.0, 9.0, 9.0);
        p(0) = 99.0;

        auto controls = one_iteration_controls();
        controls.initialization.mode = InitializationMode::Uniform;
        controls.initialization.uniform_velocity = Vec3{0.25, -0.5, 0.75};
        controls.initialization.uniform_pressure = 7.0;

        cfdx::io::DatRestartFields fields{};
        apply_initialization(
            controls.initialization, m, U, p, fields, std::string{});

        EXPECT_NEAR(U(0, 0), 0.25, 0.0);
        EXPECT_NEAR(U(0, 1), -0.5, 0.0);
        EXPECT_NEAR(U(0, 2), 0.75, 0.0);
        EXPECT_NEAR(p(0), 7.0, 0.0);
    });

    run_case("provided_initialization_preserves_caller_state", [] {
        const Mesh m = make_single_cell_cube();
        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p(1, "p", "Pa", 1);
        U.set(0, 0.125, 0.25, -0.5);
        p(0) = 3.5;

        InitializationControls provided;  // default mode
        EXPECT_TRUE(provided.mode == InitializationMode::Provided);
        cfdx::io::DatRestartFields fields{};
        apply_initialization(provided, m, U, p, fields, std::string{});

        EXPECT_NEAR(U(0, 0), 0.125, 0.0);
        EXPECT_NEAR(U(0, 1), 0.25, 0.0);
        EXPECT_NEAR(U(0, 2), -0.5, 0.0);
        EXPECT_NEAR(p(0), 3.5, 0.0);
    });

    run_case("restart_initialization_loads_checkpoint", [] {
        const Mesh m = make_single_cell_cube();
        const auto dat = std::filesystem::temp_directory_path() /
                         "cfdx_init_strategy_restart.dat";

        Field<double, Location::CELL> seed_u(1, "U", "m/s", 3);
        Field<double, Location::CELL> seed_p(1, "p", "Pa", 1);
        seed_u.set(0, 0.3, -0.2, 0.1);
        seed_p(0) = 41.0;
        cfdx::io::write_dat_restart(dat.string(), m, seed_u, seed_p, 5, 1.25);

        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p(1, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);

        InitializationControls restart;
        restart.mode = InitializationMode::Restart;
        restart.restart_path = dat.string();
        cfdx::io::DatRestartFields fields{};
        apply_initialization(restart, m, U, p, fields, std::string{});

        EXPECT_NEAR(U(0, 0), 0.3, 1e-14);
        EXPECT_NEAR(U(0, 1), -0.2, 1e-14);
        EXPECT_NEAR(U(0, 2), 0.1, 1e-14);
        EXPECT_NEAR(p(0), 41.0, 1e-14);

        std::filesystem::remove(dat);
    });

    run_case("restart_without_checkpoint_fails_loudly", [] {
        const Mesh m = make_single_cell_cube();
        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p(1, "p", "Pa", 1);

        InitializationControls restart;
        restart.mode = InitializationMode::Restart;  // declared, but no path
        cfdx::io::DatRestartFields fields{};
        // No silent fallback to the caller's state: the capability is missing.
        EXPECT_THROW(apply_initialization(
            restart, m, U, p, fields, std::string{}), std::invalid_argument);
    });

    run_case("restart_rejects_missing_checkpoint_file", [] {
        const Mesh m = make_single_cell_cube();
        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p(1, "p", "Pa", 1);

        InitializationControls restart;
        restart.mode = InitializationMode::Restart;
        restart.restart_path =
            (std::filesystem::temp_directory_path() / "cfdx_absent_checkpoint.dat").string();
        cfdx::io::DatRestartFields fields{};
        // A declared checkpoint that cannot be read must fail rather than
        // silently continue from the caller's state.
        EXPECT_THROW(apply_initialization(
            restart, m, U, p, fields, std::string{}), std::runtime_error);
    });

    run_case("declared_uniform_ignores_stray_legacy_restart_path", [] {
        const Mesh m = make_single_cell_cube();
        const auto dat = std::filesystem::temp_directory_path() /
                         "cfdx_init_legacy_path.dat";

        Field<double, Location::CELL> seed_u(1, "U", "m/s", 3);
        Field<double, Location::CELL> seed_p(1, "p", "Pa", 1);
        seed_u.set(0, 5.0, 5.0, 5.0);
        seed_p(0) = 55.0;
        cfdx::io::write_dat_restart(dat.string(), m, seed_u, seed_p, 1, 0.0);

        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p(1, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);

        // A declared Uniform strategy must not be silently overridden by a
        // positional path argument: that would make the initial state depend on
        // an input the numerical contract does not mention.
        InitializationControls uniform;
        uniform.mode = InitializationMode::Uniform;
        uniform.uniform_velocity = Vec3{1.0, 1.0, 1.0};
        uniform.uniform_pressure = 2.0;
        cfdx::io::DatRestartFields fields{};
        apply_initialization(uniform, m, U, p, fields, dat.string());

        EXPECT_NEAR(U(0, 0), 1.0, 0.0);
        EXPECT_NEAR(p(0), 2.0, 0.0);

        std::filesystem::remove(dat);
    });

    run_case("restart_mode_accepts_positional_checkpoint_path", [] {
        const Mesh m = make_single_cell_cube();
        const auto dat = std::filesystem::temp_directory_path() /
                         "cfdx_init_positional.dat";

        Field<double, Location::CELL> seed_u(1, "U", "m/s", 3);
        Field<double, Location::CELL> seed_p(1, "p", "Pa", 1);
        seed_u.set(0, 0.6, 0.0, -0.6);
        seed_p(0) = 17.0;
        cfdx::io::write_dat_restart(dat.string(), m, seed_u, seed_p, 2, 0.5);

        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p(1, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);

        // Declared Restart with the path supplied positionally: existing callers
        // that only pass the positional argument must keep working.
        InitializationControls restart;
        restart.mode = InitializationMode::Restart;
        cfdx::io::DatRestartFields fields{};
        apply_initialization(restart, m, U, p, fields, dat.string());

        EXPECT_NEAR(U(0, 0), 0.6, 1e-14);
        EXPECT_NEAR(U(0, 2), -0.6, 1e-14);
        EXPECT_NEAR(p(0), 17.0, 1e-14);

        std::filesystem::remove(dat);
    });

    run_case("declared_strategy_changes_the_residual_reference", [] {
        // The initial state fixes result.reference_momentum_residual, so two
        // different declared strategies are numerically different methods even
        // when both converge. This is why the strategy must be declared.
        const Mesh m = make_single_cell_cube();
        auto ubc = solid_walls();
        auto pbc = zero_pressure_gradient();

        auto uniform_controls = one_iteration_controls();
        uniform_controls.initialization.mode = InitializationMode::Uniform;
        uniform_controls.initialization.uniform_velocity = Vec3{0.0, 0.0, 0.0};
        uniform_controls.initialization.uniform_pressure = 0.0;

        auto provided_controls = one_iteration_controls();
        provided_controls.initialization.mode = InitializationMode::Provided;

        Field<double, Location::CELL> u_a(1, "U", "m/s", 3), p_a(1, "p", "Pa", 1);
        u_a.fill(0.0);
        p_a.fill(0.0);
        const auto uniform_result = solve_steady_incompressible(
            m, u_a, p_a, ubc, pbc, uniform_controls);

        Field<double, Location::CELL> u_b(1, "U", "m/s", 3), p_b(1, "p", "Pa", 1);
        u_b.set(0, 4.0, 0.0, 0.0);
        p_b(0) = 8.0;

        // Preservation of the caller's state is asserted on the initialization
        // step itself; the solve below legitimately advances the field.
        {
            Field<double, Location::CELL> probe_u = u_b;
            Field<double, Location::CELL> probe_p = p_b;
            cfdx::io::DatRestartFields no_fields{};
            apply_initialization(
                provided_controls.initialization, m, probe_u, probe_p,
                no_fields, std::string{});
            EXPECT_NEAR(probe_u(0, 0), 4.0, 0.0);
            EXPECT_NEAR(probe_p(0), 8.0, 0.0);
        }

        const auto provided_result = solve_steady_incompressible(
            m, u_b, p_b, ubc, pbc, provided_controls);

        // Two declared strategies are numerically different methods: the
        // reference residual that every relative residual is measured against is
        // taken from the first iteration, which depends on the initial state.
        EXPECT_TRUE(provided_result.reference_momentum_residual !=
                    uniform_result.reference_momentum_residual);
    });

    return run_all();
}
