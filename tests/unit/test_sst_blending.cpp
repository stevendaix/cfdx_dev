#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/sst_strain_rate.h"
#include "cfdx/physics/sst_case_contract.h"
#include "cfdx/io/cfdx_io/case_schema.h"
#include <array>
#include "common/test_harness.h"
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <limits>
#include <map>
#include <utility>
#include <vector>

using namespace cfdx::testing;
using namespace cfdx::physics;

static cfdx::core::Mesh make_sst_wall_operator_test_mesh()
{
    using namespace cfdx::core;
    Mesh mesh;
    mesh.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{4,0,0},{4,2,0},{0,2,0},
        {0,0,2},{4,0,2},{4,2,2},{0,2,2}
    };
    for (std::size_t i = 0; i < 8; ++i)
        mesh.points().set(i, p[i][0], p[i][1], p[i][2]);
    const std::vector<std::vector<FaceIndex>> faces = {
        {0,4,7,3}, {1,2,6,5}, {0,1,5,4},
        {3,7,6,2}, {0,3,2,1}, {4,5,6,7}
    };
    for (const auto& face : faces)
        mesh.faces().push_face(face);
    mesh.ownership().resize(faces.size());
    for (std::size_t f = 0; f < faces.size(); ++f) {
        mesh.ownership().set_owner(f, 0);
        mesh.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0,1,2,3,4,5});
    const std::array<std::pair<const char*, PatchType>, 6> patches{{
        {"wall", PatchType::WALL}, {"inlet", PatchType::INLET},
        {"side_low", PatchType::SYMMETRY}, {"side_high", PatchType::SYMMETRY},
        {"bottom", PatchType::SYMMETRY}, {"top", PatchType::SYMMETRY}
    }};
    for (std::size_t f = 0; f < patches.size(); ++f) {
        Patch patch;
        patch.name = patches[f].first;
        patch.type = patches[f].second;
        patch.face_ids = {static_cast<FaceIndex>(f)};
        mesh.boundary().add_patch(patch);
    }
    return mesh;
}



static cfdx::core::Mesh make_sst_transaction_unit_cube()
{
    using namespace cfdx::core;
    Mesh mesh;
    mesh.points().resize(8);
    const double points[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i = 0; i < 8; ++i)
        mesh.points().set(i, points[i][0], points[i][1], points[i][2]);

    const std::vector<std::vector<FaceIndex>> faces = {
        {0,3,2,1}, {4,5,6,7}, {0,1,5,4},
        {3,7,6,2}, {0,4,7,3}, {1,2,6,5}
    };
    for (const auto& face : faces)
        mesh.faces().push_face(face);

    mesh.ownership().resize(faces.size());
    for (std::size_t face = 0; face < faces.size(); ++face) {
        mesh.ownership().set_owner(face, 0);
        mesh.ownership().set_neighbour(face, FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0,1,2,3,4,5});

    Patch x0;
    x0.name = "x0";
    x0.type = PatchType::WALL;
    x0.face_ids = {4};
    mesh.boundary().add_patch(x0);
    Patch x1;
    x1.name = "x1";
    x1.type = PatchType::WALL;
    x1.face_ids = {5};
    mesh.boundary().add_patch(x1);
    Patch walls;
    walls.name = "walls";
    walls.type = PatchType::WALL;
    walls.face_ids = {0,1,2,3};
    mesh.boundary().add_patch(walls);
    return mesh;
}

int main()
{
    run_case("sst_dynamic_blending_diagnostics_accumulate_all_bound_events", [] {
        TurbulenceTransportResult aggregate;
        aggregate.k_floor_events = 2;
        aggregate.second_floor_events = 3;
        aggregate.nonfinite_bound_events = 1;
        aggregate.k_total_floor_correction = 0.25;
        aggregate.second_total_floor_correction = 0.5;
        aggregate.k_max_floor_correction = 0.2;
        aggregate.second_max_floor_correction = 0.4;

        TurbulenceTransportResult step;
        step.k_floor_events = 4;
        step.second_floor_events = 5;
        step.nonfinite_bound_events = 2;
        step.k_total_floor_correction = 0.75;
        step.second_total_floor_correction = 1.5;
        step.k_max_floor_correction = 0.6;
        step.second_max_floor_correction = 1.2;

        accumulate_turbulence_diagnostics(aggregate, step);
        EXPECT_TRUE(aggregate.k_floor_events == 6);
        EXPECT_TRUE(aggregate.second_floor_events == 8);
        EXPECT_TRUE(aggregate.nonfinite_bound_events == 3);
        EXPECT_NEAR(aggregate.k_total_floor_correction, 1.0, 0.0);
        EXPECT_NEAR(aggregate.second_total_floor_correction, 2.0, 0.0);
        EXPECT_NEAR(aggregate.k_max_floor_correction, 0.6, 0.0);
        EXPECT_NEAR(aggregate.second_max_floor_correction, 1.2, 0.0);
    });

    run_case("sst_blending_is_bounded", [] {
        const auto [f1, f2] = compute_sst_blending(0.5, 10.0, 0.01, 1.5e-5, 0.09);
        EXPECT_TRUE(f1 >= 0.0 && f1 <= 1.0);
        EXPECT_TRUE(f2 >= 0.0 && f2 <= 1.0);
    });

    run_case("sst_blending_responds_to_wall_distance", [] {
        const auto near_wall = compute_sst_blending(0.5, 10.0, 1.0e-3, 1.5e-5, 0.09);
        const auto far_wall = compute_sst_blending(0.5, 10.0, 10.0, 1.5e-5, 0.09);
        EXPECT_TRUE(near_wall.first > far_wall.first);
        EXPECT_TRUE(near_wall.second > far_wall.second);
    });

    run_case("sst_production_viscosity_uses_cellwise_sst_closure", [] {
        cfdx::core::Field<double, cfdx::core::Location::CELL> k(2, "k", "m2/s2", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> omega(2, "omega", "1/s", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> strain(2, "strain", "1/s", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> f2(2, "F2", "1", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> nu_t(2, "nu_t", "m2/s", 1);
        k(0) = 0.5; k(1) = 1.0;
        omega(0) = 10.0; omega(1) = 20.0;
        strain(0) = 2.0; strain(1) = 8.0;
        f2(0) = 0.5; f2(1) = 1.0;
        nu_t.fill(7.0);
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.a1 = 0.31;

        refresh_sst_production_eddy_viscosity(k, omega, strain, f2, controls, nu_t);
        EXPECT_NEAR(nu_t(0), 0.31 * 0.5 / std::max(0.31 * 10.0, 2.0 * 0.5), 1e-14);
        EXPECT_NEAR(nu_t(1), 0.31 * 1.0 / std::max(0.31 * 20.0, 8.0), 1e-14);
        EXPECT_TRUE(nu_t(0) != nu_t(1));
    });

    run_case("sst_production_viscosity_rejects_invalid_state_atomically", [] {
        cfdx::core::Field<double, cfdx::core::Location::CELL> k(2), omega(2), strain(2), f2(2), nu_t(2);
        k(0) = 0.5; k(1) = 1.0;
        omega(0) = 10.0; omega(1) = 20.0;
        strain(0) = 2.0; strain(1) = 4.0;
        f2(0) = 0.5; f2(1) = std::numeric_limits<double>::quiet_NaN();
        nu_t(0) = 3.0; nu_t(1) = 4.0;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        EXPECT_THROW(refresh_sst_production_eddy_viscosity(k, omega, strain, f2, controls, nu_t),
                     std::invalid_argument);
        EXPECT_NEAR(nu_t(0), 3.0, 0.0);
        EXPECT_NEAR(nu_t(1), 4.0, 0.0);
    });

    run_case("sst_transaction_preflight_rejects_field_mismatch_without_mutation", [] {
        cfdx::core::Mesh mesh;
        cfdx::core::Field<double, cfdx::core::Location::CELL> k(1), omega(1), strain(1), wall(1), nu_t(1);
        cfdx::core::Field<double, cfdx::core::Location::FACE> flux(0);
        k(0) = 0.4; omega(0) = 8.0; strain(0) = 1.0; wall(0) = 0.01; nu_t(0) = 0.02;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        EXPECT_THROW(advance_sst_transport_transactionally(
            mesh, FvGeometry{}, flux, k, omega, strain, wall, controls, nu_t),
            std::invalid_argument);
        EXPECT_NEAR(k(0), 0.4, 0.0);
        EXPECT_NEAR(omega(0), 8.0, 0.0);
        EXPECT_NEAR(nu_t(0), 0.02, 0.0);
    });


    run_case("sst_transaction_rolls_back_fields_after_nonconverged_transport", [] {
        using cfdx::core::FaceOwnership;
        using cfdx::core::Mesh;
        using cfdx::core::Patch;
        using cfdx::physics::build_fv_geometry;

        Mesh mesh;
        mesh.points().resize(8);
        const double points[8][3] = {
            {0,0,0},{1,0,0},{1,1,0},{0,1,0},
            {0,0,1},{1,0,1},{1,1,1},{0,1,1}
        };
        for (std::size_t i = 0; i < 8; ++i)
            mesh.points().set(i, points[i][0], points[i][1], points[i][2]);
        mesh.faces().push_face({0,3,2,1});
        mesh.faces().push_face({4,5,6,7});
        mesh.faces().push_face({0,1,5,4});
        mesh.faces().push_face({3,7,6,2});
        mesh.faces().push_face({0,4,7,3});
        mesh.faces().push_face({1,2,6,5});
        mesh.ownership().resize(6);
        for (std::size_t f = 0; f < 6; ++f) {
            mesh.ownership().set_owner(f, 0);
            mesh.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        }
        mesh.cells().push_cell({0,1,2,3,4,5});
        Patch wall;
        wall.name = "wall";
        wall.type = cfdx::core::PatchType::WALL;
        wall.face_ids = {0,1,2,3,4,5};
        mesh.boundary().add_patch(wall);
        const auto geometry = build_fv_geometry(mesh);

        cfdx::core::Field<double, cfdx::core::Location::FACE> flux(
            mesh.n_faces(), "phi", "kg/s", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> k(
            mesh.n_cells(), "k", "m2/s2", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> omega(
            mesh.n_cells(), "omega", "1/s", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> strain(
            mesh.n_cells(), "strain", "1/s", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> wall_distance(
            mesh.n_cells(), "wall_distance", "m", 1);
        cfdx::core::Field<double, cfdx::core::Location::CELL> nu_t(
            mesh.n_cells(), "nu_t", "m2/s", 1);
        flux.fill(0.0);
        k.fill(0.4);
        omega.fill(8.0);
        strain.fill(1.0);
        wall_distance.fill(0.1);
        nu_t.fill(0.02);

        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        ScalarBoundaryConditions k_bcs;
        k_bcs["wall"] = {ScalarBoundaryType::FIXED_VALUE, 0.01, 0.0};
        ScalarBoundaryConditions omega_bcs;
        omega_bcs["wall"] = {ScalarBoundaryType::FIXED_VALUE, 2.0, 0.0};

        auto direct_k = k;
        auto direct_omega = omega;
        const auto direct = solve_sst_transport_dynamic_blending(
            mesh, geometry, flux, direct_k, direct_omega, strain, wall_distance,
            controls, k_bcs, omega_bcs, 1, 1e-12);
        EXPECT_TRUE(!direct.converged);
        EXPECT_TRUE(direct_k(0) != k(0) || direct_omega(0) != omega(0));

        const auto step = advance_sst_transport_transactionally(
            mesh, geometry, flux, k, omega, strain, wall_distance, controls, nu_t,
            k_bcs, omega_bcs, 1, 1e-12);
        EXPECT_TRUE(!step.committed);
        EXPECT_TRUE(!step.transport.converged);
        EXPECT_NEAR(k(0), 0.4, 0.0);
        EXPECT_NEAR(omega(0), 8.0, 0.0);
        EXPECT_NEAR(nu_t(0), 0.02, 0.0);

        // The transport solver currently accepts this field's cell count and
        // reads its first component; the production-closure contract rejects
        // the non-scalar dimension after transport has converged. This forces
        // the transactional wrapper's exception rollback path without a
        // production-only fault-injection hook.
        cfdx::core::Field<double, cfdx::core::Location::CELL> malformed_strain(
            mesh.n_cells(), "strain_bad_dimension", "1/s", 2);
        malformed_strain.fill(1.0);
        auto exception_direct_k = k;
        auto exception_direct_omega = omega;
        const auto exception_direct = solve_sst_transport_dynamic_blending(
            mesh, geometry, flux, exception_direct_k, exception_direct_omega,
            malformed_strain, wall_distance, controls, k_bcs, omega_bcs, 100, 1e-4);
        EXPECT_TRUE(exception_direct.converged);
        EXPECT_TRUE(exception_direct_k(0) != k(0) ||
                    exception_direct_omega(0) != omega(0));

        EXPECT_THROW(advance_sst_transport_transactionally(
            mesh, geometry, flux, k, omega, malformed_strain, wall_distance,
            controls, nu_t, k_bcs, omega_bcs, 100, 1e-4), std::invalid_argument);
        EXPECT_NEAR(k(0), 0.4, 0.0);
        EXPECT_NEAR(omega(0), 8.0, 0.0);
        EXPECT_NEAR(nu_t(0), 0.02, 0.0);
    });

    run_case("sst_strain_rate_magnitude_matches_simple_shear", [] {
        // Ux = gamma*y, Uy = Uz = 0 gives |S| = |gamma|.
        const double gamma = 3.25;
        const std::array<double, 9> grad{
            0.0, gamma, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0};
        EXPECT_NEAR(sst_strain_rate_magnitude(grad), gamma, 1e-14);
    });

    run_case("sst_strain_rate_magnitude_matches_extension", [] {
        // Ux = a*x and Uy = -a*y gives |S| = 2*|a|.
        const double a = 1.75;
        const std::array<double, 9> grad{
            a, 0.0, 0.0,
            0.0, -a, 0.0,
            0.0, 0.0, 0.0};
        EXPECT_NEAR(sst_strain_rate_magnitude(grad), 2.0 * a, 1e-14);
    });

    run_case("sst_strain_rate_magnitude_rejects_nonfinite_gradient", [] {
        const std::array<double, 9> grad{
            0.0, 0.0, 0.0,
            0.0, std::numeric_limits<double>::infinity(), 0.0,
            0.0, 0.0, 0.0};
        EXPECT_THROW(sst_strain_rate_magnitude(grad), std::invalid_argument);
    });

    run_case("sst_strain_rate_uses_explicit_fixed_velocity_boundaries", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_wall_operator_test_mesh();
        Field<double, Location::CELL> velocity(mesh.n_cells(), "U", "m/s", 3);
        velocity.fill(0.0);

        const auto zero_gradient = compute_sst_strain_rate(velocity, mesh);
        const std::map<std::string, std::array<double, 3>> fixed_values{
            {"inlet", {1.0, 0.0, 0.0}}
        };
        const auto boundary_aware =
            compute_sst_strain_rate_with_boundary_values(velocity, mesh, fixed_values);

        EXPECT_NEAR(zero_gradient(0), 0.0, 1e-14);
        EXPECT_TRUE(std::isfinite(boundary_aware(0)));
        EXPECT_TRUE(boundary_aware(0) > 0.0);
    });

    run_case("sst_strain_rate_reconstructs_affine_field_with_fixed_boundaries", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_wall_operator_test_mesh();
        Field<double, Location::CELL> velocity(mesh.n_cells(), "U", "m/s", 3);
        velocity.fill(0.0);
        velocity.component_data(0)[0] = 2.0;

        // Ux=x is exactly representable on this orthogonal box. Prescribing
        // the affine field at every boundary face should recover dUx/dx=1.
        const std::map<std::string, std::array<double, 3>> fixed_values{
            {"wall", {0.0, 0.0, 0.0}},
            {"inlet", {4.0, 0.0, 0.0}},
            {"side_low", {2.0, 0.0, 0.0}},
            {"side_high", {2.0, 0.0, 0.0}},
            {"bottom", {2.0, 0.0, 0.0}},
            {"top", {2.0, 0.0, 0.0}}
        };

        const auto strain =
            compute_sst_strain_rate_with_boundary_values(velocity, mesh, fixed_values);
        EXPECT_NEAR(strain(0), std::sqrt(2.0), 1e-13);
    });

    run_case("sst_strain_rate_boundary_values_reject_nonfinite_input", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_wall_operator_test_mesh();
        Field<double, Location::CELL> velocity(mesh.n_cells(), "U", "m/s", 3);
        velocity.fill(0.0);
        const std::map<std::string, std::array<double, 3>> fixed_values{
            {"inlet", {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0}}
        };
        EXPECT_THROW(
            compute_sst_strain_rate_with_boundary_values(velocity, mesh, fixed_values),
            std::invalid_argument);
    });

    run_case("sst_strain_rate_rejects_unknown_fixed_value_patch", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_wall_operator_test_mesh();
        Field<double, Location::CELL> velocity(mesh.n_cells(), "U", "m/s", 3);
        velocity.fill(0.0);
        const std::map<std::string, std::array<double, 3>> fixed_values{
            {"missing_patch", {1.0, 0.0, 0.0}}
        };
        EXPECT_THROW(
            compute_sst_strain_rate_with_boundary_values(velocity, mesh, fixed_values),
            std::invalid_argument);
    });

    run_case("sst_blending_rejects_invalid_distance", [] {
        EXPECT_THROW(compute_sst_blending(0.5, 10.0, 0.0, 1.5e-5, 0.09), std::invalid_argument);
    });


    run_case("sst_inlet_turbulence_state_uses_velocity_intensity_and_length_scale", [] {
        cfdx::io::BoundarySpec inlet;
        inlet.patch_name = "inlet";
        inlet.type = cfdx::io::BCType::INLET;
        inlet.velocity_vector = {3.0, 4.0, 0.0};
        inlet.turbulence_intensity = 0.05;
        inlet.turbulence_length_scale = 0.2;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        const auto state = derive_sst_inlet_turbulence_state(inlet, controls);
        const double expected_k = 1.5 * std::pow(5.0 * 0.05, 2);
        const double expected_omega = std::sqrt(expected_k) /
            (std::pow(controls.beta_star, 0.25) * 0.2);
        EXPECT_NEAR(state.k, expected_k, 1e-14);
        EXPECT_NEAR(state.omega, expected_omega, 1e-14);
    });

    run_case("sst_inlet_turbulence_state_rejects_percentage_intensity", [] {
        cfdx::io::BoundarySpec inlet;
        inlet.type = cfdx::io::BCType::INLET;
        inlet.velocity_magnitude = 10.0;
        inlet.turbulence_intensity = 5.0;
        inlet.turbulence_length_scale = 0.1;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        EXPECT_THROW(derive_sst_inlet_turbulence_state(inlet, controls),
                     std::invalid_argument);
    });

    run_case("sst_case_prerequisites_require_explicit_initial_k_omega", [] {
        cfdx::io::CaseSetup setup;
        setup.turbulence_model = "k_omega_sst";
        setup.has_initial_condition = true;
        setup.initial_condition.scalar_fields["k"] = 0.01;
        setup.initial_condition.scalar_fields["omega"] = 2.0;
        cfdx::io::BoundarySpec inlet;
        inlet.patch_name = "inlet";
        inlet.type = cfdx::io::BCType::INLET;
        inlet.velocity_magnitude = 10.0;
        inlet.turbulence_intensity = 0.05;
        inlet.turbulence_length_scale = 0.1;
        setup.boundary_conditions.push_back(inlet);
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        validate_sst_case_prerequisites(setup, controls);

        setup.initial_condition.scalar_fields.erase("omega");
        EXPECT_THROW(validate_sst_case_prerequisites(setup, controls),
                     std::invalid_argument);
    });

    run_case("sst_case_prerequisites_reject_nonpositive_initial_fields_and_duplicate_patches", [] {
        cfdx::io::CaseSetup setup;
        setup.turbulence_model = "k_omega_sst";
        setup.has_initial_condition = true;
        setup.initial_condition.scalar_fields["k"] = 0.01;
        setup.initial_condition.scalar_fields["omega"] = 2.0;
        cfdx::io::BoundarySpec wall;
        wall.patch_name = "wall";
        wall.type = cfdx::io::BCType::WALL;
        setup.boundary_conditions.push_back(wall);
        setup.boundary_conditions.push_back(wall);
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        EXPECT_THROW(validate_sst_case_prerequisites(setup, controls),
                     std::invalid_argument);
        setup.boundary_conditions.clear();
        setup.initial_condition.scalar_fields["k"] = 0.0;
        EXPECT_THROW(validate_sst_case_prerequisites(setup, controls),
                     std::invalid_argument);
        setup.initial_condition.scalar_fields["k"] = 0.01;
        cfdx::io::BoundarySpec unknown;
        unknown.patch_name = "mystery";
        unknown.type = cfdx::io::BCType::UNKNOWN;
        setup.boundary_conditions.push_back(unknown);
        EXPECT_THROW(validate_sst_case_prerequisites(setup, controls),
                     std::invalid_argument);
    });

    run_case("sst_scalar_boundary_contract_accepts_explicit_wall_inlet_outlet", [] {
        cfdx::io::CaseSetup setup;
        cfdx::io::BoundarySpec inlet;
        inlet.patch_name = "inlet";
        inlet.type = cfdx::io::BCType::INLET;
        inlet.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 0.01},
            {"omega", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 2.0}
        };
        cfdx::io::BoundarySpec wall;
        wall.patch_name = "wall";
        wall.type = cfdx::io::BCType::WALL;
        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::WALL_OMEGA_SST, 0.0}
        };
        cfdx::io::BoundarySpec outlet;
        outlet.patch_name = "outlet";
        outlet.type = cfdx::io::BCType::PRESSURE_OUTLET;
        outlet.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::ZERO_GRADIENT, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::ZERO_GRADIENT, 0.0}
        };
        setup.boundary_conditions = {inlet, wall, outlet};
        validate_sst_scalar_boundary_contract(setup);

        setup.boundary_conditions.push_back(outlet);
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup), std::invalid_argument);
        setup.boundary_conditions.clear();
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup), std::invalid_argument);
    });

    run_case("sst_scalar_boundary_contract_rejects_zero_fixed_k", [] {
        cfdx::io::CaseSetup setup;
        cfdx::io::BoundarySpec inlet;
        inlet.patch_name = "inlet";
        inlet.type = cfdx::io::BCType::INLET;
        inlet.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 2.0}
        };
        setup.boundary_conditions = {inlet};
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup), std::invalid_argument);
    });

    run_case("sst_scalar_boundary_contract_rejects_missing_duplicate_and_incompatible_conditions", [] {
        cfdx::io::CaseSetup setup;
        cfdx::io::BoundarySpec wall;
        wall.patch_name = "wall";
        wall.type = cfdx::io::BCType::WALL;
        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0}
        };
        setup.boundary_conditions = {wall};
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup), std::invalid_argument);

        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0},
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::WALL_OMEGA_SST, 0.0}
        };
        setup.boundary_conditions = {wall};
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup), std::invalid_argument);

        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 0.01},
            {"omega", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 2.0}
        };
        setup.boundary_conditions = {wall};
        EXPECT_THROW(validate_sst_scalar_boundary_contract(setup), std::invalid_argument);
    });


    run_case("sst_production_state_validates_owned_fields", [] {
        SSTProductionState state(2);
        state.k(0) = 0.1; state.k(1) = 0.2;
        state.omega(0) = 10.0; state.omega(1) = 20.0;
        state.nu_t(0) = 0.001; state.nu_t(1) = 0.002;
        state.wall_distance(0) = 0.01; state.wall_distance(1) = 0.02;
        state.strain_rate(0) = 1.0; state.strain_rate(1) = 2.0;

        validate_sst_production_state(state, 2);
        EXPECT_THROW(validate_sst_production_state(state, 3), std::invalid_argument);

        state.wall_distance(1) = std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(validate_sst_production_state(state, 2), std::invalid_argument);
    });

    run_case("sst_production_state_snapshot_restores_all_owned_fields", [] {
        SSTProductionState state(2);
        state.k(0) = 0.1; state.k(1) = 0.2;
        state.omega(0) = 10.0; state.omega(1) = 20.0;
        state.nu_t(0) = 0.001; state.nu_t(1) = 0.002;
        state.wall_distance(0) = 0.01; state.wall_distance(1) = 0.02;
        state.strain_rate(0) = 1.0; state.strain_rate(1) = 2.0;

        const auto snapshot = snapshot_sst_production_state(state);
        state.k(0) = 9.0;
        state.omega(1) = 99.0;
        state.nu_t(0) = 0.9;
        state.wall_distance(1) = 9.0;
        state.strain_rate(0) = 99.0;

        restore_sst_production_state(state, snapshot);
        EXPECT_NEAR(state.k(0), 0.1, 0.0);
        EXPECT_NEAR(state.k(1), 0.2, 0.0);
        EXPECT_NEAR(state.omega(0), 10.0, 0.0);
        EXPECT_NEAR(state.omega(1), 20.0, 0.0);
        EXPECT_NEAR(state.nu_t(0), 0.001, 0.0);
        EXPECT_NEAR(state.nu_t(1), 0.002, 0.0);
        EXPECT_NEAR(state.wall_distance(0), 0.01, 0.0);
        EXPECT_NEAR(state.wall_distance(1), 0.02, 0.0);
        EXPECT_NEAR(state.strain_rate(0), 1.0, 0.0);
        EXPECT_NEAR(state.strain_rate(1), 2.0, 0.0);
    });

    run_case("sst_coupled_transaction_rolls_back_carrier_flux_and_all_turbulence_fields", [] {
        using cfdx::core::Field;
        Field<double, cfdx::core::Location::CELL> U(2, "U", "m/s", 3);
        Field<double, cfdx::core::Location::CELL> p(2, "p", "Pa", 1);
        Field<double, cfdx::core::Location::FACE> phi(3, "phi", "kg/s", 1);
        SSTProductionState state(2);
        for (std::size_t c = 0; c < 2; ++c) {
            U.component_data(0)[c] = 1.0 + c;
            U.component_data(1)[c] = 2.0 + c;
            U.component_data(2)[c] = 3.0 + c;
            p(c) = 10.0 + c;
            state.k(c) = 0.1 + 0.1 * c;
            state.omega(c) = 10.0 + 10.0 * c;
            state.nu_t(c) = 0.001 + 0.001 * c;
            state.wall_distance(c) = 0.01 + 0.01 * c;
            state.strain_rate(c) = 1.0 + c;
        }
        phi(0) = 0.25; phi(1) = -0.5; phi(2) = 0.75;

        SSTCoupledStateTransaction transaction(U, p, phi, state);
        transaction.begin();
        U.component_data(0)[0] = -100.0;
        p(1) = -200.0;
        phi(2) = 300.0;
        state.k(0) = 400.0;
        state.omega(1) = 500.0;
        state.nu_t(0) = 600.0;
        state.wall_distance(1) = 700.0;
        state.strain_rate(0) = 800.0;
        transaction.rollback();

        EXPECT_TRUE(!transaction.active());
        EXPECT_NEAR(U.component_data(0)[0], 1.0, 0.0);
        EXPECT_NEAR(p(1), 11.0, 0.0);
        EXPECT_NEAR(phi(2), 0.75, 0.0);
        EXPECT_NEAR(state.k(0), 0.1, 0.0);
        EXPECT_NEAR(state.omega(1), 20.0, 0.0);
        EXPECT_NEAR(state.nu_t(0), 0.001, 0.0);
        EXPECT_NEAR(state.wall_distance(1), 0.02, 0.0);
        EXPECT_NEAR(state.strain_rate(0), 1.0, 0.0);
    });

    run_case("sst_coupled_transaction_commit_keeps_all_candidate_fields", [] {
        using cfdx::core::Field;
        Field<double, cfdx::core::Location::CELL> U(1, "U", "m/s", 3);
        Field<double, cfdx::core::Location::CELL> p(1, "p", "Pa", 1);
        Field<double, cfdx::core::Location::FACE> phi(1, "phi", "kg/s", 1);
        SSTProductionState state(1);
        U.fill(1.0); p.fill(2.0); phi.fill(3.0);
        state.k(0) = 0.1; state.omega(0) = 10.0;
        state.nu_t(0) = 0.01; state.wall_distance(0) = 0.1;
        state.strain_rate(0) = 1.0;

        SSTCoupledStateTransaction transaction(U, p, phi, state);
        transaction.begin();
        U.component_data(0)[0] = 4.0;
        p(0) = 5.0;
        phi(0) = 6.0;
        state.k(0) = 0.2;
        transaction.commit();

        EXPECT_TRUE(!transaction.active());
        EXPECT_NEAR(U.component_data(0)[0], 4.0, 0.0);
        EXPECT_NEAR(p(0), 5.0, 0.0);
        EXPECT_NEAR(phi(0), 6.0, 0.0);
        EXPECT_NEAR(state.k(0), 0.2, 0.0);
        EXPECT_THROW(transaction.rollback(), std::logic_error);
    });


    run_case("sst_runtime_boundary_adapter_preserves_wall_omega_obligation", [] {
        cfdx::io::CaseSetup setup;
        cfdx::io::BoundarySpec wall;
        wall.patch_name = "wall";
        wall.type = cfdx::io::BCType::WALL;
        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::WALL_OMEGA_SST, 0.0}
        };
        setup.boundary_conditions.push_back(wall);

        const auto runtime = make_sst_runtime_boundary_contract(setup);
        EXPECT_TRUE(runtime.k.count("wall") == 1);
        EXPECT_TRUE(runtime.k.at("wall").type == ScalarBoundaryType::FIXED_VALUE);
        EXPECT_NEAR(runtime.k.at("wall").value, 0.0, 0.0);
        EXPECT_TRUE(runtime.omega.count("wall") == 0);
        EXPECT_TRUE(runtime.wall_omega_sst_patches.size() == 1);
        EXPECT_TRUE(runtime.wall_omega_sst_patches.front() == "wall");
        EXPECT_THROW(require_supported_sst_runtime_boundaries(runtime),
                     std::invalid_argument);
    });

    run_case("sst_runtime_boundary_adapter_maps_fixed_and_zero_gradient", [] {
        cfdx::io::CaseSetup setup;
        cfdx::io::BoundarySpec inlet;
        inlet.patch_name = "inlet";
        inlet.type = cfdx::io::BCType::INLET;
        inlet.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 0.4},
            {"omega", cfdx::io::SSTScalarBoundaryType::FIXED_VALUE, 8.0}
        };
        cfdx::io::BoundarySpec outlet;
        outlet.patch_name = "outlet";
        outlet.type = cfdx::io::BCType::OUTLET;
        outlet.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::ZERO_GRADIENT, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::ZERO_GRADIENT, 0.0}
        };
        setup.boundary_conditions = {inlet, outlet};

        const auto runtime = make_sst_runtime_boundary_contract(setup);
        EXPECT_NEAR(runtime.k.at("inlet").value, 0.4, 0.0);
        EXPECT_NEAR(runtime.omega.at("inlet").value, 8.0, 0.0);
        EXPECT_TRUE(runtime.k.at("outlet").type == ScalarBoundaryType::ZERO_GRADIENT);
        EXPECT_TRUE(runtime.omega.at("outlet").type == ScalarBoundaryType::ZERO_GRADIENT);
        EXPECT_TRUE(runtime.wall_omega_sst_patches.empty());
        bool supported = true;
        try {
            require_supported_sst_runtime_boundaries(runtime);
        } catch (...) {
            supported = false;
        }
        EXPECT_TRUE(supported);
    });

    run_case("sst_wall_omega_face_operator_only_populates_requested_wall_faces", [] {
        const auto mesh = make_sst_wall_operator_test_mesh();
        const auto geometry = build_fv_geometry(mesh);
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.molecular_viscosity = 1.0e-5;
        controls.beta1 = 0.075;
        const auto face_values = make_sst_wall_omega_face_values(
            mesh, geometry, std::vector<std::string>{"wall"}, controls);
        EXPECT_TRUE(face_values.values.count("wall") == 1);
        EXPECT_TRUE(face_values.values.count("inlet") == 0);
        EXPECT_NEAR(face_values.values.at("wall")[0], 2.0e-4, 1.0e-15);
        for (std::size_t face = 1; face < mesh.n_faces(); ++face)
            EXPECT_TRUE(!std::isfinite(face_values.values.at("wall")[face]));
        const auto runtime = SSTRuntimeBoundaryContract{{}, {}, {"wall"}};
        EXPECT_THROW(require_supported_sst_runtime_boundaries(runtime),
                     std::invalid_argument);
        require_supported_sst_runtime_boundaries(runtime, &face_values);
    });

    run_case("sst_wall_omega_value_uses_resolved_wall_formula", [] {
        const double value = sst_wall_omega_value(1.0e-5, 0.075, 2.0e-3);
        EXPECT_NEAR(value, 200.0, 1.0e-10);
        EXPECT_THROW(sst_wall_omega_value(1.0e-5, 0.075, 0.0),
                     std::invalid_argument);
        EXPECT_THROW(sst_wall_omega_value(0.0, 0.075, 2.0e-3),
                     std::invalid_argument);
    });

    run_case("sst_production_state_initializes_from_explicit_case_and_mesh_fields", [] {
        auto mesh = make_sst_wall_operator_test_mesh();
        cfdx::io::CaseSetup setup;
        setup.turbulence_model = "k_omega_sst";
        setup.has_initial_condition = true;
        setup.initial_condition.scalar_fields["k"] = 0.12;
        setup.initial_condition.scalar_fields["omega"] = 18.0;
        cfdx::io::BoundarySpec wall;
        wall.patch_name = "wall";
        wall.type = cfdx::io::BCType::WALL;
        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::WALL_OMEGA_SST, 0.0}
        };
        setup.boundary_conditions = {wall};

        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.molecular_viscosity = 1.0e-5;
        const auto state = initialize_sst_production_state(
            mesh, setup, controls, std::vector<double>{0.01},
            std::vector<double>{2.0});
        validate_sst_production_state(state, mesh.n_cells());
        EXPECT_NEAR(state.k(0), 0.12, 0.0);
        EXPECT_NEAR(state.omega(0), 18.0, 0.0);
        EXPECT_NEAR(state.wall_distance(0), 0.01, 0.0);
        EXPECT_NEAR(state.strain_rate(0), 2.0, 0.0);
        EXPECT_TRUE(std::isfinite(state.nu_t(0)) && state.nu_t(0) >= 0.0);

        EXPECT_THROW(initialize_sst_production_state(
            mesh, setup, controls, std::vector<double>{0.0},
            std::vector<double>{2.0}), std::invalid_argument);
        EXPECT_THROW(initialize_sst_production_state(
            mesh, setup, controls, std::vector<double>{0.01},
            std::vector<double>{std::numeric_limits<double>::quiet_NaN()}),
            std::invalid_argument);
    });

    run_case("sst_production_state_builds_wall_distance_and_strain_from_mesh", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_wall_operator_test_mesh();
        Field<double, Location::CELL> velocity(mesh.n_cells(), "U", "m/s", 3);
        velocity.fill(0.0);

        cfdx::io::CaseSetup setup;
        setup.physics_model = "incompressible";
        setup.turbulence_model = "SST";
        setup.has_initial_condition = true;
        setup.initial_condition.scalar_fields["k"] = 0.12;
        setup.initial_condition.scalar_fields["omega"] = 18.0;
        cfdx::io::BoundarySpec wall;
        wall.patch_name = "wall";
        wall.type = cfdx::io::BCType::WALL;
        wall.sst_scalar_conditions = {
            {"k", cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO, 0.0},
            {"omega", cfdx::io::SSTScalarBoundaryType::WALL_OMEGA_SST, 0.0}
        };
        setup.boundary_conditions = {wall};

        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.molecular_viscosity = 1.0e-5;
        const std::map<std::string, std::array<double, 3>> fixed_values{
            {"inlet", {1.0, 0.0, 0.0}}
        };

        auto state = initialize_sst_production_state_from_mesh(
            mesh, setup, controls, velocity, fixed_values);
        validate_sst_production_state(state, mesh.n_cells());
        EXPECT_NEAR(state.k(0), 0.12, 0.0);
        EXPECT_NEAR(state.omega(0), 18.0, 0.0);
        EXPECT_TRUE(std::isfinite(state.wall_distance(0)) &&
                    state.wall_distance(0) > 0.0);
        EXPECT_TRUE(std::isfinite(state.strain_rate(0)) &&
                    state.strain_rate(0) > 0.0);
        EXPECT_TRUE(std::isfinite(state.nu_t(0)) && state.nu_t(0) >= 0.0);

        const double initial_strain = state.strain_rate(0);
        const std::map<std::string, std::array<double, 3>> updated_values{
            {"inlet", {2.0, 0.0, 0.0}}
        };
        refresh_sst_production_derived_state_from_mesh(
            state, mesh, controls, velocity, updated_values);
        EXPECT_TRUE(state.strain_rate(0) > initial_strain);
        validate_sst_production_state(state, mesh.n_cells());

        const double valid_strain = state.strain_rate(0);
        const std::map<std::string, std::array<double, 3>> invalid_values{
            {"missing_patch", {1.0, 0.0, 0.0}}
        };
        EXPECT_THROW(refresh_sst_production_derived_state_from_mesh(
            state, mesh, controls, velocity, invalid_values), std::invalid_argument);
        EXPECT_NEAR(state.strain_rate(0), valid_strain, 0.0);

        Field<double, Location::CELL> invalid_velocity(mesh.n_cells(), "bad", "m/s", 1);
        EXPECT_THROW(initialize_sst_production_state_from_mesh(
            mesh, setup, controls, invalid_velocity, fixed_values),
            std::invalid_argument);
    });

    run_case("sst_dat_restart_restores_named_scalars_and_rebuilds_derived_state", [] {
        using cfdx::core::Field;
        using cfdx::core::Location;
        const auto mesh = make_sst_wall_operator_test_mesh();
        const auto path = std::filesystem::temp_directory_path() /
            "cfdx_sst_restart_hydration.dat";
        std::filesystem::remove(path);
        std::filesystem::remove(path.string() + ".tmp");

        Field<double, Location::CELL> saved_u(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> saved_p(mesh.n_cells(), "p", "Pa", 1);
        Field<double, Location::CELL> saved_k(mesh.n_cells(), "k", "m2/s2", 1);
        Field<double, Location::CELL> saved_omega(mesh.n_cells(), "omega", "1/s", 1);
        saved_u.component_data(0)[0] = 0.0;
        saved_u.component_data(1)[0] = 0.0;
        saved_u.component_data(2)[0] = 0.0;
        saved_p(0) = 14.0;
        saved_k(0) = 0.42;
        saved_omega(0) = 125.0;
        cfdx::io::DatRestartFields write_fields{};
        write_fields.k = &saved_k;
        write_fields.omega = &saved_omega;
        cfdx::io::write_dat_restart_fields(
            path.string(), mesh, saved_u, saved_p, write_fields, 17, 2.5);

        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.component_data(0)[0] = 3.0;
        U.component_data(1)[0] = -2.0;
        U.component_data(2)[0] = 1.0;
        p(0) = -9.0;
        SSTProductionState state(mesh.n_cells());
        state.k(0) = 9.0;
        state.omega(0) = 9.0;
        state.nu_t(0) = 9.0;
        state.wall_distance(0) = 0.5;

        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.molecular_viscosity = 1.5e-5;
        const std::map<std::string, std::array<double, 3>> fixed_velocity_boundary_values{
            {"inlet", {4.0, 0.0, 0.0}}
        };
        const auto metadata = restore_sst_production_state_from_dat(
            path.string(), mesh, U, p, state, std::vector<double>{0.01}, controls,
            fixed_velocity_boundary_values);

        EXPECT_TRUE(metadata.iteration == 17);
        EXPECT_NEAR(metadata.time, 2.5, 0.0);
        EXPECT_NEAR(U.component_data(0)[0], 0.0, 0.0);
        EXPECT_NEAR(p(0), 14.0, 0.0);
        EXPECT_NEAR(state.k(0), 0.42, 0.0);
        EXPECT_NEAR(state.omega(0), 125.0, 0.0);
        EXPECT_NEAR(state.wall_distance(0), 0.01, 0.0);
        EXPECT_TRUE(std::isfinite(state.strain_rate(0)) && state.strain_rate(0) > 0.0);
        EXPECT_TRUE(std::isfinite(state.nu_t(0)) && state.nu_t(0) >= 0.0);
        validate_sst_production_state(state, mesh.n_cells());
        std::filesystem::remove(path);
    });

    run_case("sst_dat_restart_missing_omega_leaves_carrier_and_state_unchanged", [] {
        using cfdx::core::Field;
        using cfdx::core::Location;
        const auto mesh = make_sst_wall_operator_test_mesh();
        const auto path = std::filesystem::temp_directory_path() /
            "cfdx_sst_restart_missing_omega.dat";
        std::filesystem::remove(path);
        std::filesystem::remove(path.string() + ".tmp");

        Field<double, Location::CELL> saved_u(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> saved_p(mesh.n_cells(), "p", "Pa", 1);
        Field<double, Location::CELL> saved_k(mesh.n_cells(), "k", "m2/s2", 1);
        saved_u.component_data(0)[0] = 1.0;
        saved_u.component_data(1)[0] = 2.0;
        saved_u.component_data(2)[0] = 3.0;
        saved_p(0) = 4.0;
        saved_k(0) = 0.2;
        cfdx::io::DatRestartFields write_fields{};
        write_fields.k = &saved_k;
        cfdx::io::write_dat_restart_fields(
            path.string(), mesh, saved_u, saved_p, write_fields, 4, 1.0);

        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        U.component_data(0)[0] = 8.0;
        U.component_data(1)[0] = 7.0;
        U.component_data(2)[0] = 6.0;
        p(0) = 5.0;
        SSTProductionState state(mesh.n_cells());
        state.k(0) = 0.3;
        state.omega(0) = 6.0;
        state.nu_t(0) = 0.04;
        state.wall_distance(0) = 0.01;
        state.strain_rate(0) = 2.0;
        const auto old_state = state;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;

        EXPECT_THROW(restore_sst_production_state_from_dat(
            path.string(), mesh, U, p, state, std::vector<double>{0.01}, controls),
            std::runtime_error);
        EXPECT_NEAR(U.component_data(0)[0], 8.0, 0.0);
        EXPECT_NEAR(U.component_data(1)[0], 7.0, 0.0);
        EXPECT_NEAR(U.component_data(2)[0], 6.0, 0.0);
        EXPECT_NEAR(p(0), 5.0, 0.0);
        EXPECT_NEAR(state.k(0), old_state.k(0), 0.0);
        EXPECT_NEAR(state.omega(0), old_state.omega(0), 0.0);
        EXPECT_NEAR(state.nu_t(0), old_state.nu_t(0), 0.0);
        EXPECT_NEAR(state.wall_distance(0), old_state.wall_distance(0), 0.0);
        EXPECT_NEAR(state.strain_rate(0), old_state.strain_rate(0), 0.0);
        std::filesystem::remove(path);
    });

    run_case("sst_transaction_rejects_non_scalar_transport_fields", [] {
        using namespace cfdx::core;
        const Mesh mesh = make_sst_transaction_unit_cube();
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::FACE> flux(mesh.n_faces(), "phi", "kg/s", 1);
        Field<double, Location::CELL> k(mesh.n_cells(), "k", "m2/s2", 1);
        Field<double, Location::CELL> omega(mesh.n_cells(), "omega", "1/s", 1);
        Field<double, Location::CELL> strain(mesh.n_cells(), "strain", "1/s", 1);
        Field<double, Location::CELL> wall(mesh.n_cells(), "wall_distance", "m", 1);
        Field<double, Location::CELL> nu_t(mesh.n_cells(), "nu_t", "m2/s", 1);
        Field<double, Location::CELL> malformed_k(mesh.n_cells(), "k_bad", "m2/s2", 2);
        flux.fill(0.0);
        k(0) = 0.4;
        omega(0) = 8.0;
        strain(0) = 1.0;
        wall(0) = 0.01;
        nu_t(0) = 0.02;
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;

        EXPECT_THROW(advance_sst_transport_transactionally(
            mesh, geometry, flux, malformed_k, omega, strain, wall,
            controls, nu_t, {}, {}, 1, 1e-8), std::invalid_argument);

        Field<double, Location::FACE> malformed_flux(
            mesh.n_faces(), "phi_bad", "kg/s", 2);
        EXPECT_THROW(advance_sst_transport_transactionally(
            mesh, geometry, malformed_flux, k, omega, strain, wall,
            controls, nu_t, {}, {}, 1, 1e-8), std::invalid_argument);
    });

    run_case("sst_transaction_rolls_back_k_omega_after_nonconverged_transport", [] {
        using namespace cfdx::core;
        const Mesh mesh = make_sst_transaction_unit_cube();
        const auto geometry = build_fv_geometry(mesh);
        Field<double, Location::FACE> flux(mesh.n_faces(), "phi", "kg/s", 1);
        Field<double, Location::CELL> k(mesh.n_cells(), "k", "m2/s2", 1);
        Field<double, Location::CELL> omega(mesh.n_cells(), "omega", "1/s", 1);
        Field<double, Location::CELL> strain(mesh.n_cells(), "strain", "1/s", 1);
        Field<double, Location::CELL> wall(mesh.n_cells(), "wall_distance", "m", 1);
        Field<double, Location::CELL> nu_t(mesh.n_cells(), "nu_t", "m2/s", 1);
        flux.fill(0.0);
        k(0) = 0.4;
        omega(0) = 8.0;
        strain(0) = 1.0;
        wall(0) = 0.01;
        nu_t(0) = 0.02;

        ScalarBoundaryConditions k_bcs{
            {"x0", {ScalarBoundaryType::FIXED_VALUE, 0.9, 0.0}},
            {"x1", {ScalarBoundaryType::FIXED_VALUE, 0.2, 0.0}},
            {"walls", {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0}}
        };
        ScalarBoundaryConditions omega_bcs{
            {"x0", {ScalarBoundaryType::FIXED_VALUE, 20.0, 0.0}},
            {"x1", {ScalarBoundaryType::FIXED_VALUE, 4.0, 0.0}},
            {"walls", {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0}}
        };
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;

        const auto step = advance_sst_transport_transactionally(
            mesh, geometry, flux, k, omega, strain, wall, controls, nu_t,
            k_bcs, omega_bcs, 1, 1e-30);

        EXPECT_TRUE(!step.committed);
        EXPECT_TRUE(!step.transport.converged);
        EXPECT_TRUE(step.transport.iterations == 1);
        EXPECT_NEAR(k(0), 0.4, 0.0);
        EXPECT_NEAR(omega(0), 8.0, 0.0);
        EXPECT_NEAR(nu_t(0), 0.02, 0.0);
    });


    run_case("sst_production_iterate_transaction_rolls_back_derived_and_carrier_state", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_transaction_unit_cube();
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        Field<double, Location::FACE> phi(mesh.n_faces(), "phi", "kg/s", 1);
        U.fill(0.0);
        U.component_data(0)[0] = 0.5;
        p.fill(2.0);
        phi.fill(0.0);

        SSTProductionState state(mesh.n_cells());
        state.k(0) = 0.4;
        state.omega(0) = 8.0;
        state.nu_t(0) = 0.02;
        state.wall_distance(0) = 0.25;
        state.strain_rate(0) = 3.0;
        const auto original_state = snapshot_sst_production_state(state);
        const auto original_U = U;
        const auto original_p = p;
        const auto original_phi = phi;

        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.molecular_viscosity = 1.0e-5;
        const std::map<std::string, std::array<double, 3>> fixed_values{
            {"x0", {1.0, 0.0, 0.0}}
        };
        const ScalarBoundaryConditions k_bcs{
            {"x0", {ScalarBoundaryType::FIXED_VALUE, 0.9, 0.0}},
            {"x1", {ScalarBoundaryType::FIXED_VALUE, 0.2, 0.0}},
            {"walls", {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0}}
        };
        const ScalarBoundaryConditions omega_bcs{
            {"x0", {ScalarBoundaryType::FIXED_VALUE, 20.0, 0.0}},
            {"x1", {ScalarBoundaryType::FIXED_VALUE, 4.0, 0.0}},
            {"walls", {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0}}
        };

        const auto result = advance_sst_production_iterate_transactionally(
            mesh, U, p, phi, state, controls, fixed_values, {}, k_bcs,
            omega_bcs, 1, 1e-30);

        EXPECT_TRUE(!result.committed);
        EXPECT_TRUE(!result.transport.converged);
        EXPECT_TRUE(result.reason.find("fully rolled back") != std::string::npos);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            EXPECT_NEAR(U.component_data(0)[c], original_U.component_data(0)[c], 0.0);
            EXPECT_NEAR(U.component_data(1)[c], original_U.component_data(1)[c], 0.0);
            EXPECT_NEAR(U.component_data(2)[c], original_U.component_data(2)[c], 0.0);
            EXPECT_NEAR(p(c), original_p(c), 0.0);
            EXPECT_NEAR(state.k(c), original_state.k(c), 0.0);
            EXPECT_NEAR(state.omega(c), original_state.omega(c), 0.0);
            EXPECT_NEAR(state.nu_t(c), original_state.nu_t(c), 0.0);
            EXPECT_NEAR(state.wall_distance(c), original_state.wall_distance(c), 0.0);
            EXPECT_NEAR(state.strain_rate(c), original_state.strain_rate(c), 0.0);
        }
        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            EXPECT_NEAR(phi(f), original_phi(f), 0.0);
    });


    run_case("sst_relative_eddy_viscosity_change_is_dimensionless", [] {
        SSTProductionState previous(2);
        SSTProductionState current(2);
        for (std::size_t cell = 0; cell < 2; ++cell) {
            previous.k(cell) = current.k(cell) = 0.1;
            previous.omega(cell) = current.omega(cell) = 1.0;
            previous.wall_distance(cell) = current.wall_distance(cell) = 0.5;
            previous.strain_rate(cell) = current.strain_rate(cell) = 0.0;
        }
        previous.nu_t(0) = 1e-5;
        current.nu_t(0) = 2e-5;
        previous.nu_t(1) = current.nu_t(1) = 1e-3;

        EXPECT_NEAR(
            sst_relative_eddy_viscosity_change(previous, current, 1e-5),
            0.5, 1e-14);
        EXPECT_NEAR(
            sst_relative_eddy_viscosity_change(current, current, 1e-5),
            0.0, 0.0);
        EXPECT_THROW(
            sst_relative_eddy_viscosity_change(previous, current, 0.0),
            std::invalid_argument);
    });

    run_case("sst_production_runtime_requires_explicit_valid_state_and_excludes_prescribed_nu_t", [] {
        IncompressibleSolverControls controls;
        SSTProductionRuntime runtime;
        controls.sst_runtime = &runtime;

        EXPECT_THROW(validate_incompressible_controls(controls, 1), std::invalid_argument);

        SSTProductionState state(1);
        state.k(0) = 0.1;
        state.omega(0) = 1.0;
        state.nu_t(0) = 0.01;
        state.wall_distance(0) = 0.5;
        state.strain_rate(0) = 0.2;
        runtime.state = &state;
        runtime.transport.model = TurbulenceModel::SST;
        runtime.transport.density = controls.density;
        runtime.transport.molecular_viscosity = controls.kinematic_viscosity;
        validate_incompressible_controls(controls, 1);

        runtime.transport.density *= 2.0;
        EXPECT_THROW(validate_incompressible_controls(controls, 1), std::invalid_argument);
        runtime.transport.density = controls.density;

        controls.turbulent_viscosity = 1.0e-5;
        EXPECT_THROW(validate_incompressible_controls(controls, 1), std::invalid_argument);
    });

    run_case("sst_outer_iteration_transaction_rolls_back_after_inner_sst_commit", [] {
        using namespace cfdx::core;
        const auto mesh = make_sst_transaction_unit_cube();
        Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
        Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
        Field<double, Location::FACE> phi(mesh.n_faces(), "phi", "kg/s", 1);
        U.fill(0.0);
        U.component_data(0)[0] = 0.25;
        p.fill(3.0);
        phi.fill(0.0);

        SSTProductionState state(mesh.n_cells());
        state.k.fill(0.4);
        state.omega.fill(8.0);
        state.nu_t.fill(0.02);
        state.wall_distance.fill(0.1);
        state.strain_rate.fill(1.0);

        const auto original_U = U;
        const auto original_p = p;
        const auto original_phi = phi;
        const auto original_sst = snapshot_sst_production_state(state);

        // The outer transaction remains active across a successfully committed
        // inner SST transaction. A later carrier rejection must still restore
        // the entire coupled iterate to the pre-outer-iteration snapshot.
        SSTCoupledStateTransaction outer(U, p, phi, state);
        outer.begin();
        {
            SSTCoupledStateTransaction inner(U, p, phi, state);
            inner.begin();
            U.component_data(0)[0] = 1.25;
            p(0) = 9.0;
            phi(0) = 2.0;
            state.k(0) = 0.8;
            state.omega(0) = 12.0;
            state.nu_t(0) = 0.05;
            state.wall_distance(0) = 0.2;
            state.strain_rate(0) = 4.0;
            inner.commit();
        }

        EXPECT_NEAR(U.component_data(0)[0], 1.25, 0.0);
        EXPECT_NEAR(state.k(0), 0.8, 0.0);
        outer.rollback();

        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            for (std::size_t d = 0; d < 3; ++d)
                EXPECT_NEAR(U.component_data(d)[c], original_U.component_data(d)[c], 0.0);
            EXPECT_NEAR(p(c), original_p(c), 0.0);
            EXPECT_NEAR(state.k(c), original_sst.k(c), 0.0);
            EXPECT_NEAR(state.omega(c), original_sst.omega(c), 0.0);
            EXPECT_NEAR(state.nu_t(c), original_sst.nu_t(c), 0.0);
            EXPECT_NEAR(state.wall_distance(c), original_sst.wall_distance(c), 0.0);
            EXPECT_NEAR(state.strain_rate(c), original_sst.strain_rate(c), 0.0);
        }
        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            EXPECT_NEAR(phi(f), original_phi(f), 0.0);
    });

    return run_all();
}
