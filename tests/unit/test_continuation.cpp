#include "cfdx/physics/continuation.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

using namespace cfdx::physics;
using namespace cfdx::testing;


static cfdx::core::Mesh make_unit_cube_with_lid()
{
    using namespace cfdx::core;
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for (std::size_t i=0; i<8; ++i)
        m.points().set(i,p[i][0],p[i][1],p[i][2]);

    m.faces().push_face({0,3,2,1}); // bottom
    m.faces().push_face({4,5,6,7}); // top
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for (std::size_t f=0; f<6; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});

    Patch bottom{"bottom",PatchType::WALL,{0}};
    Patch top{"top",PatchType::WALL,{1}};
    Patch sides{"sides",PatchType::WALL,{2,3,4,5}};
    m.boundary().add_patch(bottom);
    m.boundary().add_patch(top);
    m.boundary().add_patch(sides);
    return m;
}

int main()
{
    run_case("continuation_controls_validation", [] {
        ContinuationControls c;
        validate_continuation_controls(c);
        EXPECT_NEAR(continuation_next_target(0.0, 0.25), 0.25, 1e-12);
        EXPECT_NEAR(continuation_next_target(0.9, 0.25), 1.0, 1e-12);
    });

    run_case("continuation_step_policy", [] {
        ContinuationControls c;
        c.initial_step = 0.25;
        c.maximum_step = 0.5;
        c.minimum_step = 0.01;
        c.step_growth = 1.5;
        c.step_reduction = 0.5;
        EXPECT_NEAR(continuation_step_after_success(0.25, c), 0.375, 1e-12);
        EXPECT_NEAR(continuation_step_after_success(0.5, c), 0.5, 1e-12);
        EXPECT_NEAR(continuation_step_after_failure(0.25, c), 0.125, 1e-12);
    });

    run_case("continuation_invalid_controls", [] {
        ContinuationControls c;
        c.step_reduction = 1.0;
        bool rejected = false;
        try {
            validate_continuation_controls(c);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        EXPECT_TRUE(rejected);
    });


    run_case("continuation_integrates_with_steady_solver", [] {
        const auto mesh = make_unit_cube_with_lid();
        cfdx::core::Field<double,cfdx::core::Location::CELL> U(1,"U","m/s",3);
        cfdx::core::Field<double,cfdx::core::Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(0.0);

        VelocityBoundaryConditions ubc;
        ubc["bottom"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ubc["top"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{1.0,0.0,0.0}};
        ubc["sides"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["bottom"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        pbc["top"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        pbc["sides"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls controls;
        controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
        controls.convergence.max_iterations = 50;
        controls.convergence.relative_tolerance = 1e-8;
        controls.convergence.continuity_tolerance = 1e-8;
        controls.linear_tolerance = 1e-10;
        controls.pressure_reference_cell = 0;
        controls.pressure_reference_value = 0.0;

        ContinuationControls continuation;
        continuation.enabled = true;
        continuation.initial_step = 0.5;
        continuation.minimum_step = 0.25;
        continuation.maximum_step = 0.5;
        continuation.max_stage_attempts = 4;
        continuation.max_stages = 8;

        const auto result = solve_steady_incompressible_continuation(
            mesh,U,p,ubc,pbc,controls,continuation);

        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(result.final_parameter,1.0,1e-14);
        EXPECT_TRUE(result.total_attempts >= 2);
        EXPECT_TRUE(result.stages.size() >= 2);
        for (const auto& stage : result.stages) {
            EXPECT_TRUE(stage.converged);
            EXPECT_TRUE(stage.status == cfdx::core::ConvergenceStatus::CONVERGED);
            EXPECT_TRUE(std::isfinite(stage.parameter));
            EXPECT_TRUE(std::isfinite(stage.step));
        }
    });

    run_case("continuation_step_arithmetic_is_bounded", [] {
        ContinuationControls c;
        c.minimum_step = 0.05;
        c.maximum_step = 0.5;
        c.step_growth = 1.5;
        c.step_reduction = 0.5;

        EXPECT_NEAR(continuation_next_target(0.0, 0.25), 0.25, 1e-14);
        EXPECT_NEAR(continuation_next_target(0.9, 0.25), 1.0, 1e-14);

        EXPECT_NEAR(continuation_step_after_success(0.25, c), 0.375, 1e-14);
        EXPECT_NEAR(continuation_step_after_success(0.4, c), 0.5, 1e-14);

        EXPECT_NEAR(continuation_step_after_failure(0.25, c), 0.125, 1e-14);
        EXPECT_NEAR(continuation_step_after_failure(0.05, c), 0.025, 1e-14);
    });

    run_case("continuation_rolls_back_every_rejected_stage_exactly", [] {
        auto mesh = make_unit_cube_with_lid();
        cfdx::core::Field<double, cfdx::core::Location::CELL> U(1,"U","m/s",3);
        cfdx::core::Field<double, cfdx::core::Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["bottom"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ubc["top"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{1.0,0.0,0.0}};
        ubc["sides"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["bottom"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        pbc["top"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        pbc["sides"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        IncompressibleSolverControls controls;
        controls.body_force = cfdx::core::Vec3{0.0, 0.0, -1.0};
        controls.coupling.alpha_u = 0.7;
        controls.coupling.alpha_p = 0.3;
        controls.convergence.relative_tolerance = 1e-10;
        controls.convergence.continuity_tolerance = 1e-8;
        controls.linear_tolerance = 1e-10;
        controls.pressure_reference_cell = 0;
        controls.pressure_reference_value = 0.0;
        // One outer corrector cannot satisfy the nonlinear criteria, so every
        // stage attempt is rejected and the continuation must give up at the
        // minimum step without ever leaking a partial iterate.
        controls.convergence.max_iterations = 1;

        const auto u_before = U;
        const auto p_before = p;

        ContinuationControls continuation;
        continuation.enabled = true;
        continuation.initial_step = 0.5;
        continuation.minimum_step = 0.125;
        continuation.maximum_step = 0.5;
        continuation.max_stage_attempts = 8;
        continuation.max_stages = 4;

        const auto result = solve_steady_incompressible_continuation(
            mesh, U, p, ubc, pbc, controls, continuation);

        EXPECT_TRUE(!result.converged);
        EXPECT_TRUE(!result.stages.empty());
        for (const auto& stage : result.stages)
            EXPECT_TRUE(!stage.converged);
        EXPECT_NEAR(result.stages.back().step, continuation.minimum_step, 1e-14);

        // Exact rollback: no rejected attempt may leave a trace in the fields.
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            EXPECT_TRUE(p(c) == p_before(c));
            for (std::size_t d = 0; d < 3; ++d)
                EXPECT_TRUE(U.component_data(d)[c] == u_before.component_data(d)[c]);
        }
    });

    return run_all();
}
