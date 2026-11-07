#include "cfdx/physics/steady_incompressible_solver.h"
#include "audit/sst_channel_fixture.h"
#include "cfdx/physics/pressure_velocity_system.h"
#include "cfdx/io/probe/probe_csv.h"
#include "cfdx/io/restart/dat_restart.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include "common/test_harness.h"
#include <limits>
#include <tuple>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_unit_cube()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i=0; i<8; ++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
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
    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,1,2,3,4,5};
    m.boundary().add_patch(wall);
    return m;
}


static Mesh make_two_cell_channel()
{
    Mesh m;
    m.points().resize(12);
    const double p[12][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i=0; i<12; ++i)
        m.points().set(i,p[i][0],p[i][1],p[i][2]);

    m.faces().push_face({0,3,7,4});       // x = 0
    m.faces().push_face({1,2,6,5});       // x = 1, internal
    m.faces().push_face({8,9,11,10});     // x = 2
    m.faces().push_face({0,1,5,4});       // y = 0, left
    m.faces().push_face({3,7,6,2});       // y = 1, left
    m.faces().push_face({0,3,2,1});       // z = 0, left
    m.faces().push_face({4,5,6,7});       // z = 1, left
    m.faces().push_face({1,8,10,5});      // y = 0, right
    m.faces().push_face({2,6,11,9});      // y = 1, right
    m.faces().push_face({1,2,9,8});       // z = 0, right
    m.faces().push_face({5,10,11,6});     // z = 1, right

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const std::size_t owner = (f == 2 || f >= 7) ? 1u : 0u;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);

    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
}

int main()
{
    run_case("n9_2_rhie_chow_pressure_sensitivity_and_flux_reference", [] {
        const Mesh m = make_two_cell_channel();
        const auto geometry = build_fv_geometry(m);

        Field<double,Location::CELL> U(2,"U","m/s",3);
        Field<double,Location::CELL> p0(2,"p","Pa",1);
        Field<double,Location::CELL> p1(2,"p","Pa",1);
        U.fill(0.0);
        U.set(0,2.0,0.0,0.0);
        U.set(1,2.0,0.0,0.0);
        p0.fill(0.0);
        p1.fill(0.0);
        p1(1) = 2.0;

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        const std::array<std::vector<double>,3> rAU{
            std::vector<double>{1.0,1.0},
            std::vector<double>{1.0,1.0},
            std::vector<double>{1.0,1.0}};

        const auto phi_cell = make_mass_flux(
            m, geometry, U, 1.0, ubc);
        const auto phi0 = make_rhie_chow_mass_flux(
            m, geometry, U, p0, rAU, 1.0, ubc, pbc);
        const auto phi1 = make_rhie_chow_mass_flux(
            m, geometry, U, p1, rAU, 1.0, ubc, pbc);

        EXPECT_NEAR(phi_cell(1), 2.0, 1e-12);
        EXPECT_NEAR(phi0(1), 2.0, 1e-12);
        EXPECT_NEAR(phi1(1), 0.0, 1e-12);
        EXPECT_TRUE(std::abs(phi0(1)) > std::abs(phi1(1)));
        EXPECT_NEAR(phi1(1), 0.0, 1e-12);
        EXPECT_NEAR(phi1(0), 0.0, 1e-14);
        EXPECT_NEAR(phi1(2), 0.0, 1e-14);

        // The authoritative flux is antisymmetric across the internal face:
        // owner cell sees phi, neighbour cell sees -phi.
        double corrected_div0 = phi1(1);
        double corrected_div1 = -phi1(1);
        EXPECT_NEAR(corrected_div0 + corrected_div1, 0.0, 1e-14);
        EXPECT_NEAR(std::abs(corrected_div0), 0.0, 1e-12);

        const auto contract = PressureVelocitySystemContract{
            6, 2,
            PressureGaugePolicy::REFERENCE_CELL,
            ContinuityFluxPolicy::RHIE_CHOW,
            PressureStabilizationPolicy::RHIE_CHOW,
            true, true, true};
        contract.validate();
    });


    run_case("actual_pressure_replay_is_post_gauge_and_opt_in", [] {
        const auto mesh = make_two_cell_channel();
        Field<double,Location::CELL> u(2,"U","m/s",3);
        Field<double,Location::CELL> p(2,"p","Pa",1);
        u.fill(0.0); p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        IncompressibleSolverControls controls;
        controls.convergence.max_iterations = 2;
        controls.pressure_reference_cell = 0;
        controls.pressure_replay_iterations = {1};
        const auto prefix = std::filesystem::temp_directory_path() / "cfdx_real_pressure_capture";
        controls.pressure_replay_prefix = prefix.string();
        const auto file = prefix.string() + "_i1_c0.cfdxlin";
        std::filesystem::remove(file);
        (void)solve_steady_incompressible(mesh,u,p,ubc,pbc,controls);
        EXPECT_TRUE(std::filesystem::exists(file));
        const auto captured = cfdx::core::read_linear_system_capture(file);
        // The two-cell pure Neumann operator is gauge-reduced to 1x1.
        EXPECT_TRUE(captured.matrix.n_rows() == 1);
        EXPECT_TRUE(captured.matrix.n_cols() == 1);
        EXPECT_TRUE(captured.rhs.size() == 1);
        EXPECT_NEAR(captured.initial_solution(0), 0.0, 0.0);
        EXPECT_TRUE(captured.problem == LinearProblemKind::PressurePoisson);
        std::filesystem::remove(file);
    });

    run_case("native_solver_consumes_dat_restart", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> seed_u(1,"U","m/s",3);
        Field<double,Location::CELL> seed_p(1,"p","Pa",1);
        seed_u.set(0,0.25,-0.15,0.05);
        seed_p(0) = 37.5;

        const auto dat = std::filesystem::temp_directory_path() / "cfdx_native_restart_test.dat";
        cfdx::io::write_dat_restart(dat.string(), m, seed_u, seed_p, 17, 2.5);

        Field<double,Location::CELL> loaded_u(1,"U","m/s",3);
        Field<double,Location::CELL> loaded_p(1,"p","Pa",1);
        const auto state = cfdx::io::read_dat_restart(dat.string(), m, loaded_u, loaded_p);
        EXPECT_TRUE(state.iteration == std::size_t{17});
        EXPECT_NEAR(state.time, 2.5, 1e-14);
        EXPECT_NEAR(loaded_u(0,0), 0.25, 1e-14);
        EXPECT_NEAR(loaded_u(0,1), -0.15, 1e-14);
        EXPECT_NEAR(loaded_u(0,2), 0.05, 1e-14);
        EXPECT_NEAR(loaded_p(0), 37.5, 1e-14);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls controls;
        controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
        controls.convergence.max_iterations = 1;
        controls.convergence.continuity_tolerance = 1e-12;
        controls.linear_tolerance = 1e-12;
        controls.pressure_reference_cell = 0;
        controls.pressure_reference_value = 0.0;

        Field<double,Location::CELL> expected_u = loaded_u;
        Field<double,Location::CELL> expected_p = loaded_p;
        (void)solve_steady_incompressible(m, expected_u, expected_p, ubc, pbc, controls);

        // Declare the checkpoint as the initialization strategy rather than
        // relying on the positional restart argument alone.
        IncompressibleSolverControls restart_controls = controls;
        restart_controls.initialization.mode = InitializationMode::Restart;
        restart_controls.initialization.restart_path = dat.string();

        Field<double,Location::CELL> restart_u(1,"U","m/s",3);
        Field<double,Location::CELL> restart_p(1,"p","Pa",1);
        restart_u.fill(0.0);
        restart_p.fill(0.0);
        (void)solve_steady_incompressible(
            m, restart_u, restart_p, ubc, pbc, restart_controls);

        EXPECT_NEAR(restart_u(0,0), expected_u(0,0), 1e-14);
        EXPECT_NEAR(restart_u(0,1), expected_u(0,1), 1e-14);
        EXPECT_NEAR(restart_u(0,2), expected_u(0,2), 1e-14);
        EXPECT_NEAR(restart_p(0), expected_p(0), 1e-14);

        std::filesystem::remove(dat);
    });

    run_case("native_solver_consumes_dat_v2_thermal_turbulence_restart", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> seed_u(1,"U","m/s",3), seed_p(1,"p","Pa",1);
        Field<double,Location::CELL> seed_T(1,"T","K",1), seed_k(1,"k","m2/s2",1), seed_omega(1,"omega","1/s",1);
        seed_u.set(0,0.2,0.0,0.0); seed_p(0)=12.0;
        seed_T(0)=350.0; seed_k(0)=0.25; seed_omega(0)=3.0;
        const auto dat = std::filesystem::temp_directory_path() / "cfdx_native_v2_restart_test.dat";
        cfdx::io::write_dat_restart_fields(
            dat.string(),m,seed_u,seed_p,{&seed_T,&seed_k,&seed_omega},9,0.75);

        Field<double,Location::CELL> U(1,"U","m/s",3), p(1,"p","Pa",1);
        Field<double,Location::CELL> T(1,"T","K",1), k(1,"k","m2/s2",1), omega(1,"omega","1/s",1);
        U.fill(0.0); p.fill(0.0); T.fill(0.0); k.fill(0.0); omega.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        IncompressibleSolverControls controls;
        controls.convergence.max_iterations = 1;
        controls.linear_tolerance = 1e-12;

        Field<double,Location::CELL> expected_u(1,"U","m/s",3), expected_p(1,"p","Pa",1);
        expected_u.set(0,0.2,0.0,0.0); expected_p(0)=12.0;
        (void)solve_steady_incompressible(
            m,expected_u,expected_p,ubc,pbc,controls);

        IncompressibleSolverControls restart_controls = controls;
        restart_controls.initialization.mode = InitializationMode::Restart;
        restart_controls.initialization.restart_path = dat.string();
        cfdx::io::DatRestartFields restart_fields{&T,&k,&omega};
        (void)solve_steady_incompressible(
            m,U,p,ubc,pbc,restart_controls,"",restart_fields);

        EXPECT_NEAR(U(0,0),expected_u(0,0),1e-14);
        EXPECT_NEAR(U(0,1),expected_u(0,1),1e-14);
        EXPECT_NEAR(U(0,2),expected_u(0,2),1e-14);
        EXPECT_NEAR(p(0),expected_p(0),1e-14);
        EXPECT_NEAR(T(0),350.0,1e-14);
        EXPECT_NEAR(k(0),0.25,1e-14);
        EXPECT_NEAR(omega(0),3.0,1e-14);
        std::filesystem::remove(dat);
    });

    run_case("native_point_probe_emission_is_iteration_synchronized", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3), p(1,"p","Pa",1);
        U.set(0,2.0,-1.0,2.0);
        p(0)=7.0;
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        IncompressibleSolverControls controls;
        controls.convergence.max_iterations = 2;
        controls.linear_tolerance = 1e-12;
        controls.probes = {
            {"pressure_probe", {0.5,0.5,0.5}, IncompressibleProbeField::PRESSURE},
            {"velocity_probe", {0.5,0.5,0.5}, IncompressibleProbeField::U_MAGNITUDE}
        };
        std::vector<IncompressibleProbeSample> samples;
        controls.probe_callback = [&](const IncompressibleProbeSample& sample) { samples.push_back(sample); };
        (void)solve_steady_incompressible(m,U,p,ubc,pbc,controls);
        EXPECT_TRUE(samples.size() == 4);
        EXPECT_TRUE(samples[0].name == "pressure_probe");
        EXPECT_TRUE(samples[0].iteration == 1);
        EXPECT_NEAR(samples[0].time, 0.0, 1e-14);
        EXPECT_TRUE(samples[1].name == "velocity_probe");
        EXPECT_TRUE(samples[1].iteration == 1);
        EXPECT_TRUE(samples[2].name == "pressure_probe");
        EXPECT_TRUE(samples[2].iteration == 2);
        EXPECT_TRUE(samples[3].name == "velocity_probe");
        EXPECT_TRUE(samples[3].iteration == 2);
        EXPECT_TRUE(std::isfinite(samples[0].value));
        EXPECT_TRUE(std::isfinite(samples[1].value));
    });

    run_case("native_point_probe_matches_zero_reference_solution", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3), p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        IncompressibleSolverControls controls;
        controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
        controls.convergence.max_iterations = 2;
        controls.convergence.continuity_tolerance = 1e-12;
        controls.linear_tolerance = 1e-12;
        controls.pressure_reference_cell = 0;
        controls.pressure_reference_value = 0.0;
        controls.probes = {
            {"pressure_reference", {0.5,0.5,0.5}, IncompressibleProbeField::PRESSURE},
            {"velocity_reference", {0.5,0.5,0.5}, IncompressibleProbeField::U_MAGNITUDE}
        };
        std::vector<IncompressibleProbeSample> samples;
        controls.probe_callback = [&](const IncompressibleProbeSample& sample) { samples.push_back(sample); };
        (void)solve_steady_incompressible(m,U,p,ubc,pbc,controls);
        EXPECT_TRUE(samples.size() == 4);
        for (const auto& sample : samples) {
            EXPECT_TRUE(sample.iteration >= 1 && sample.iteration <= 2);
            EXPECT_NEAR(sample.time, 0.0, 1e-14);
            EXPECT_NEAR(sample.value, 0.0, 1e-14);
        }
    });

    // Closes the "compare CSV against analytical/reference values" item for the
    // probe CSV export stream: a zero-velocity, zero-gradient-pressure field on a
    // unit cube with fixed-zero walls admits the uniform zero field as its
    // analytical steady solution. The probes are exported to CSV and the values
    // read back from disk are asserted against that analytical reference.
    run_case("native_point_probe_csv_matches_analytical_reference", []() {
        const Mesh m = make_unit_cube();
        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p(1, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        IncompressibleSolverControls controls;
        controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
        controls.convergence.max_iterations = 2;
        controls.convergence.continuity_tolerance = 1e-12;
        controls.linear_tolerance = 1e-12;
        controls.pressure_reference_cell = 0;
        controls.pressure_reference_value = 0.0;
        controls.probes = {
            {"pressure_reference", {0.5, 0.5, 0.5}, IncompressibleProbeField::PRESSURE},
            {"velocity_reference", {0.5, 0.5, 0.5}, IncompressibleProbeField::U_MAGNITUDE},
        };

        const auto result = solve_steady_incompressible(m, U, p, ubc, pbc, controls);
        EXPECT_TRUE(!result.probe_samples.empty());

        const std::string csv_path = unique_temp_path(".csv").string();
        cfdx::io::write_probe_csv(csv_path, result.probe_samples, controls.probes);

        std::ifstream in(csv_path);
        EXPECT_TRUE(in.is_open());
        std::string line;
        EXPECT_TRUE(std::getline(in, line) && line == "# cfdx-probe-csv v1");
        EXPECT_TRUE(std::getline(in, line) && line == "# columns: probe,iteration,time,value");
        // The export now carries self-describing probe metadata (# probe lines),
        // ordered by probe name, matching the application-layer writer.
        EXPECT_TRUE(std::getline(in, line) && line.rfind("# probe ", 0) == 0);
        EXPECT_TRUE(line.find("pressure_reference") != std::string::npos);
        EXPECT_TRUE(line.find("field=p") != std::string::npos);
        EXPECT_TRUE(line.find("unit=Pa") != std::string::npos);
        std::size_t data_rows = 0;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            if (line[0] == '#')
                continue; // skip remaining # probe metadata lines
            const std::size_t last_comma = line.find_last_of(',');
            EXPECT_TRUE(last_comma != std::string::npos);
            std::istringstream value_stream(line.substr(last_comma + 1));
            double value = 0.0;
            value_stream >> value;
            EXPECT_TRUE(std::isfinite(value));
            // Analytical reference: the analytical solution is the uniform zero
            // field, so every exported sample must reproduce it.
            EXPECT_NEAR(value, 0.0, 1e-12);
            ++data_rows;
        }
        EXPECT_TRUE(data_rows == result.probe_samples.size());
        std::filesystem::remove(csv_path);
    });

    run_case("steady_incompressible_sst_dat_restart_hydrates_runtime_before_first_predictor", [] {
        SSTChannelCase c;
        auto seed_u = c.U;
        auto seed_p = c.p;
        auto seed_k = c.state.k;
        auto seed_omega = c.state.omega;
        seed_u.fill(0.0);
        for (std::size_t cell = 0; cell < c.mesh.n_cells(); ++cell)
            seed_u(cell, 0) = 1.0;
        seed_p.fill(0.0);
        seed_k.fill(0.02);
        seed_omega.fill(4.0);
        const auto dat = unique_temp_path("_sst_runtime_restart.dat");
        cfdx::io::write_dat_restart_fields(dat.string(), c.mesh, seed_u, seed_p,
            cfdx::io::DatRestartFields{nullptr, &seed_k, nullptr, &seed_omega},
            4, 0.2);
        c.controls.initialization.mode = InitializationMode::Restart;
        c.controls.initialization.restart_path = dat.string();
        double first_predictor_mu = std::numeric_limits<double>::quiet_NaN();
        c.controls.momentum_predictor_callback =
            [&](std::size_t, const std::vector<double>& mu,
                const Field<double, Location::CELL>&,
                const Field<double, Location::FACE>&) {
                if (!std::isfinite(first_predictor_mu))
                    first_predictor_mu = mu.at(0);
            };
        const auto accepted = solve_steady_incompressible(
            c.mesh, c.U, c.p, c.ubc, c.pbc, c.controls);
        EXPECT_TRUE(accepted.iterations == 1);
        EXPECT_NEAR(first_predictor_mu,
            c.controls.density * (c.controls.kinematic_viscosity + 0.02/4.0), 1e-12);
        const auto accepted_U = c.U;
        const auto accepted_p = c.p;
        const auto accepted_state = snapshot_sst_production_state(c.state);

        Field<double, Location::CELL> conflicting_k(3, "other_k", "m2/s2", 1);
        conflicting_k.fill(7.0);
        cfdx::io::DatRestartFields conflicting_fields;
        conflicting_fields.k = &conflicting_k;
        EXPECT_THROW_WITH(solve_steady_incompressible(
            c.mesh, c.U, c.p, c.ubc, c.pbc, c.controls, "", conflicting_fields),
            std::invalid_argument, "active runtime state");
        for (std::size_t cell = 0; cell < 3; ++cell) {
            EXPECT_NEAR(conflicting_k(cell), 7.0, 0.0);
            EXPECT_NEAR(c.p(cell), accepted_p(cell), 0.0);
            EXPECT_NEAR(c.state.k(cell), accepted_state.k(cell), 0.0);
            EXPECT_NEAR(c.state.omega(cell), accepted_state.omega(cell), 0.0);
            for (std::size_t d = 0; d < 3; ++d)
                EXPECT_NEAR(c.U(cell,d), accepted_U(cell,d), 0.0);
        }
        Field<double, Location::CELL> conflicting_omega(3, "other_omega", "1/s", 1);
        conflicting_omega.fill(9.0);
        cfdx::io::DatRestartFields wrong_omega;
        wrong_omega.omega = &conflicting_omega;
        EXPECT_THROW_WITH(solve_steady_incompressible(
            c.mesh, c.U, c.p, c.ubc, c.pbc, c.controls, "", wrong_omega),
            std::invalid_argument, "active runtime state");

        const auto incomplete = unique_temp_path("_sst_missing_omega.dat");
        cfdx::io::write_dat_restart_fields(incomplete.string(),
            c.mesh, seed_u, seed_p,
            cfdx::io::DatRestartFields{nullptr, &seed_k, nullptr, nullptr}, 5, 0.25);
        c.controls.initialization.restart_path = incomplete.string();
        EXPECT_THROW_WITH(solve_steady_incompressible(
            c.mesh, c.U, c.p, c.ubc, c.pbc, c.controls),
            std::runtime_error, "requested omega but DAT has none");
        for (std::size_t cell = 0; cell < 3; ++cell) {
            EXPECT_NEAR(c.p(cell), accepted_p(cell), 0.0);
            EXPECT_NEAR(c.state.k(cell), accepted_state.k(cell), 0.0);
            EXPECT_NEAR(c.state.omega(cell), accepted_state.omega(cell), 0.0);
            EXPECT_NEAR(c.state.nu_t(cell), accepted_state.nu_t(cell), 0.0);
            for (std::size_t d = 0; d < 3; ++d)
                EXPECT_NEAR(c.U(cell,d), accepted_U(cell,d), 0.0);
        }
        std::filesystem::remove(dat);
        std::filesystem::remove(incomplete);
    });

    run_case("steady_incompressible_sst_continuous_and_dat_restarted_trajectories_match", [] {
        SSTChannelCase continuous;
        continuous.controls.convergence.max_iterations = 2;
        const auto full = solve_steady_incompressible(
            continuous.mesh, continuous.U, continuous.p,
            continuous.ubc, continuous.pbc, continuous.controls);
        EXPECT_TRUE(full.iterations == 2);
        SSTChannelCase split;
        const auto first = solve_steady_incompressible(
            split.mesh, split.U, split.p, split.ubc, split.pbc, split.controls);
        EXPECT_TRUE(first.iterations == 1);
        auto persisted_flux = first.authoritative_mass_flux;
        const auto dat = unique_temp_path("_sst_multicell_trajectory.dat");
        cfdx::io::write_dat_restart_fields(dat.string(), split.mesh,
            split.U, split.p,
            cfdx::io::DatRestartFields{nullptr, &split.state.k, nullptr,
                                        &split.state.omega, &persisted_flux}, 4, 0.1);
        SSTChannelCase restarted;
        restarted.controls.initialization.mode = InitializationMode::Restart;
        restarted.controls.initialization.restart_path = dat.string();
        const auto resumed = solve_steady_incompressible(
            restarted.mesh, restarted.U, restarted.p,
            restarted.ubc, restarted.pbc, restarted.controls);
        EXPECT_TRUE(resumed.iterations == 1);
        for (std::size_t cell = 0; cell < 3; ++cell) {
            for (std::size_t d = 0; d < 3; ++d)
                EXPECT_NEAR(restarted.U(cell,d), continuous.U(cell,d), 1e-12);
            EXPECT_NEAR(restarted.p(cell), continuous.p(cell), 1e-12);
            EXPECT_NEAR(restarted.state.k(cell), continuous.state.k(cell), 1e-12);
            EXPECT_NEAR(restarted.state.omega(cell), continuous.state.omega(cell), 1e-12);
            EXPECT_NEAR(restarted.state.nu_t(cell), continuous.state.nu_t(cell), 1e-12);
        }
        for (std::size_t face = 0; face < continuous.mesh.n_faces(); ++face)
            EXPECT_NEAR(resumed.authoritative_mass_flux(face),
                        full.authoritative_mass_flux(face), 1e-12);
        std::filesystem::remove(dat);
    });

    run_case("steady_incompressible_opt_in_sst_runtime_executes_multicell_path", [] {
        SSTChannelCase c;
        c.controls.convergence.max_iterations = 1;
        const auto result = solve_steady_incompressible(
            c.mesh, c.U, c.p, c.ubc, c.pbc, c.controls);
        EXPECT_TRUE(result.history.size() == 1);
        EXPECT_TRUE(!result.converged);
        EXPECT_TRUE(result.convergence_status ==
                    cfdx::core::ConvergenceStatus::MAX_ITERATIONS);
        EXPECT_TRUE(result.convergence_reason.find("SST transport/closure") !=
                    std::string::npos);
        EXPECT_TRUE(result.authoritative_mass_flux.size() == c.mesh.n_faces());
        EXPECT_TRUE(c.state.k(0) > 0.0 && c.state.omega(0) > 0.0);
        EXPECT_TRUE(std::isfinite(c.state.nu_t(0)) && c.state.nu_t(0) >= 0.0);
        EXPECT_TRUE(c.state.nu_t(0) != 9.0);
    });

    run_case("steady_incompressible_sst_state_rolls_back_after_post_transport_retry", [] {
        SSTChannelCase c;
        c.controls.nonlinear_retry.max_retries = 1;
        const auto original = snapshot_sst_production_state(c.state);
        std::size_t output_calls = 0;
        std::size_t predictor_calls = 0;
        Field<double, Location::CELL> first_U;
        Field<double, Location::FACE> first_flux;
        c.controls.momentum_predictor_callback =
            [&](std::size_t, const auto&, const auto& velocity, const auto& flux) {
                ++predictor_calls;
                if (predictor_calls == 1) {
                    first_U = velocity;
                    first_flux = flux;
                } else {
                    for (std::size_t cell = 0; cell < c.mesh.n_cells(); ++cell)
                        for (std::size_t d = 0; d < 3; ++d)
                            EXPECT_NEAR(velocity(cell,d), first_U(cell,d), 0.0);
                    for (std::size_t face = 0; face < c.mesh.n_faces(); ++face)
                        EXPECT_NEAR(flux(face), first_flux(face), 0.0);
                }
            };
        c.controls.iteration_output_callback =
            [&](std::size_t, double, const auto&, const auto&, const auto&) {
                if (++output_calls == 1) {
                    EXPECT_TRUE(c.state.k(0) != original.k(0));
                    throw NonlinearRetryableFailure("deliberate post-SST rejection");
                }
                return false;
            };
        const auto accepted = solve_steady_incompressible(
            c.mesh, c.U, c.p, c.ubc, c.pbc, c.controls);
        EXPECT_TRUE(output_calls == 2 && predictor_calls == 2);
        EXPECT_TRUE(accepted.history.size() == 1);
        SSTChannelCase reference;
        reference.controls.coupling.alpha_u = accepted.history[0].effective_alpha_u;
        reference.controls.coupling.alpha_p = accepted.history[0].effective_alpha_p;
        const auto clean = solve_steady_incompressible(
            reference.mesh, reference.U, reference.p,
            reference.ubc, reference.pbc, reference.controls);
        for (std::size_t cell = 0; cell < 3; ++cell) {
            for (std::size_t d = 0; d < 3; ++d)
                EXPECT_NEAR(c.U(cell,d), reference.U(cell,d), 1e-12);
            EXPECT_NEAR(c.p(cell), reference.p(cell), 1e-12);
            EXPECT_NEAR(c.state.k(cell), reference.state.k(cell), 1e-12);
            EXPECT_NEAR(c.state.omega(cell), reference.state.omega(cell), 1e-12);
            EXPECT_NEAR(c.state.nu_t(cell), reference.state.nu_t(cell), 1e-12);
        }
        for (std::size_t face = 0; face < c.mesh.n_faces(); ++face)
            EXPECT_NEAR(accepted.authoritative_mass_flux(face),
                        clean.authoritative_mass_flux(face), 1e-12);
    });

    run_case("steady_incompressible_sst_runtime_multicell_refreshes_closure", [] {
        SSTChannelCase c;
        std::size_t predictors = 0;
        c.controls.momentum_predictor_callback =
            [&](std::size_t, const auto& mu, const auto&, const auto& flux) {
                ++predictors;
                for (std::size_t cell = 0; cell < 3; ++cell)
                    EXPECT_NEAR(mu[cell],
                        c.controls.density * (c.controls.kinematic_viscosity + 0.005), 1e-13);
                EXPECT_NEAR(flux(0), -1.2, 1e-13);
            };
        const auto result = solve_steady_incompressible(
            c.mesh, c.U, c.p, c.ubc, c.pbc, c.controls);
        EXPECT_TRUE(predictors == 1);
        EXPECT_TRUE(result.history.size() == 1);
        for (std::size_t cell = 0; cell < 3; ++cell) {
            EXPECT_TRUE(std::isfinite(c.state.k(cell)) && c.state.k(cell) > 0.0);
            EXPECT_TRUE(std::isfinite(c.state.omega(cell)) && c.state.omega(cell) > 0.0);
            EXPECT_TRUE(std::isfinite(c.state.nu_t(cell)) && c.state.nu_t(cell) >= 0.0);
            EXPECT_TRUE(c.state.nu_t(cell) != 9.0);
        }
    });

    run_case("steady_incompressible_sst_closed_unforced_state_remains_rejected", [] {
        // A wall-only, zero-strain one-cell problem has no production or inflow.
        // Its k/omega destruction leads to floors: never mask this failure by
        // weakening the transactional finite-state SST gate.
        const Mesh mesh = make_unit_cube();
        Field<double, Location::CELL> U(1, "U", "m/s", 3);
        Field<double, Location::CELL> p(1, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0,0};
        IncompressibleSolverControls controls;
        controls.convergence.max_iterations = 1;
        controls.pressure_reference_cell = 0;
        controls.nonlinear_retry.max_retries = 1;
        SSTProductionState state(1);
        state.k(0) = 0.1;
        state.omega(0) = 1.0;
        state.nu_t(0) = 0.1;
        state.wall_distance(0) = 0.5;
        SSTProductionRuntime runtime;
        runtime.state = &state;
        runtime.transport.model = TurbulenceModel::SST;
        runtime.transport.density = controls.density;
        runtime.transport.molecular_viscosity = controls.kinematic_viscosity;
        runtime.k_boundary_conditions["wall"] =
            {ScalarBoundaryType::ZERO_GRADIENT,0,0};
        runtime.omega_boundary_conditions["wall"] =
            {ScalarBoundaryType::ZERO_GRADIENT,0,0};
        runtime.wall_patch_names = {"wall"};
        controls.sst_runtime = &runtime;
        EXPECT_THROW_WITH(solve_steady_incompressible(
            mesh, U, p, ubc, pbc, controls), NonlinearRetryableFailure,
            "SST update rejected");
        EXPECT_NEAR(state.k(0), 0.1, 0.0);
        EXPECT_NEAR(state.omega(0), 1.0, 0.0);
        EXPECT_NEAR(U(0,0), 0.0, 0.0);
        EXPECT_NEAR(p(0), 0.0, 0.0);
    });

    run_case("steady_incompressible_zero_state_is_fixed_point", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(0.0);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::SIMPLE;
        c.convergence.max_iterations = 5;
        c.convergence.continuity_tolerance = 1e-12;
        c.linear_tolerance = 1e-12;
        c.pressure_reference_cell = 0;
        c.pressure_reference_value = 0.0;

        const auto r = solve_steady_incompressible(m,U,p,ubc,pbc,c);
        EXPECT_TRUE(r.converged);
        EXPECT_TRUE(!r.history.empty());
        EXPECT_NEAR(U(0,0),0.0,1e-14);
        EXPECT_NEAR(U(0,1),0.0,1e-14);
        EXPECT_NEAR(U(0,2),0.0,1e-14);
        EXPECT_NEAR(p(0),0.0,1e-14);
        EXPECT_NEAR(r.history.back().continuity_linf,0.0,1e-14);
    });


    run_case("simplec_zero_state_uses_consistent_momentum_diagonal", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(0.0);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::SIMPLEC;
        c.convergence.max_iterations = 5;
        c.convergence.continuity_tolerance = 1e-12;
        c.linear_tolerance = 1e-12;
        c.pressure_reference_cell = 0;
        c.pressure_reference_value = 0.0;

        const auto r = solve_steady_incompressible(m,U,p,ubc,pbc,c);
        EXPECT_TRUE(r.converged);
        EXPECT_TRUE(!r.history.empty());
    });

    run_case("piso_zero_state_runs_multiple_pressure_corrections", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(0.0);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::PISO;
        c.coupling.n_pressure_correctors = 3;
        c.convergence.max_iterations = 5;
        c.convergence.continuity_tolerance = 1e-12;
        c.linear_tolerance = 1e-12;
        c.pressure_reference_cell = 0;
        c.pressure_reference_value = 0.0;

        const auto r = solve_steady_incompressible(m,U,p,ubc,pbc,c);
        EXPECT_TRUE(r.converged);
        EXPECT_TRUE(!r.history.empty());
        // The zero-state pressure RHS is already satisfied by the initial
        // guess. ReusableCgContext must return before AMG setup; all three
        // pressure-correction calls are still accounted as solves.
        EXPECT_TRUE(r.pressure_linear_context.full_setups == 0);
        EXPECT_TRUE(r.pressure_linear_context.numeric_updates == 0);
        EXPECT_TRUE(r.pressure_linear_context.unchanged_reuses == 0);
        EXPECT_TRUE(
            r.pressure_linear_context.solves == 3 * r.iterations);
    });

    run_case("piso_neighbour_correction_rejects_wrong_algorithm", [] {
        const Mesh m = make_two_cell_channel();
        Field<double, Location::CELL> U(2,"U","m/s",3), p(2,"p","Pa",1);
        U.fill(0.0); p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0,0};
        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::SIMPLE;
        c.coupling.piso_neighbor_correction = true;
        c.convergence.max_iterations = 1;
        bool rejected = false;
        try { (void)solve_steady_incompressible(m,U,p,ubc,pbc,c); }
        catch (const std::invalid_argument&) { rejected = true; }
        EXPECT_TRUE(rejected);
    });

    run_case("pcd_commutator_option_requires_coupled_pcd", [] {
        IncompressibleSolverControls c;
        EXPECT_TRUE(!c.coupling.pcd_lumped_mass_commutator_experimental);
        c.coupling.pcd_lumped_mass_commutator_experimental = true;
        c.algorithm = PressureVelocityAlgorithm::SIMPLE;
        bool rejected = false;
        try { validate_incompressible_controls(c, 2); }
        catch (const std::invalid_argument&) { rejected = true; }
        EXPECT_TRUE(rejected);
        c.algorithm = PressureVelocityAlgorithm::COUPLED;
        c.coupling.schur_model = CoupledSchurModel::SIMPLE;
        rejected = false;
        try { validate_incompressible_controls(c, 2); }
        catch (const std::invalid_argument&) { rejected = true; }
        EXPECT_TRUE(rejected);
        c.coupling.schur_model = CoupledSchurModel::PCD;
        // Other controls have valid defaults; the option is now admissible.
        validate_incompressible_controls(c, 2);
    });

    run_case("schwarz_rejects_wrong_algorithm_schur_and_config", [] {
        IncompressibleSolverControls c;
        EXPECT_TRUE(!c.coupling.coupled_schwarz_experimental);
        c.coupling.coupled_schwarz_experimental = true;
        c.coupling.schur_model = CoupledSchurModel::BlockLocal;
        c.algorithm = PressureVelocityAlgorithm::SIMPLE;
        EXPECT_THROW(validate_incompressible_controls(c, 4), std::invalid_argument);
        c.algorithm = PressureVelocityAlgorithm::COUPLED;
        c.coupling.schur_model = CoupledSchurModel::PCD;
        EXPECT_THROW(validate_incompressible_controls(c, 4), std::invalid_argument);
        c.coupling.schur_model = CoupledSchurModel::BlockLocal;
        c.coupling.coupled_schwarz_max_local_unknowns = 2;
        EXPECT_THROW(validate_incompressible_controls(c, 4), std::invalid_argument);
        c.coupling.coupled_schwarz_max_local_unknowns = 16;
        c.coupling.coupled_schwarz_cells_per_subdomain = 5;
        EXPECT_THROW(validate_incompressible_controls(c, 4), std::invalid_argument);
        c.coupling.coupled_schwarz_cells_per_subdomain = 4;
        validate_incompressible_controls(c, 4);
    });

    run_case("piso_neighbour_correction_is_opt_in_and_uses_inner_velocity", [] {
        const Mesh m = make_two_cell_channel();
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0,0};

        auto run = [&](bool correction, int ncorrectors) {
            Field<double, Location::CELL> U(2,"U","m/s",3), p(2,"p","Pa",1);
            U.fill(0.0); p.fill(0.0);
            U.set(0, 1.0, 0.3, 0.0);
            U.set(1,-0.5,-0.2, 0.0);
            p(1)=3.0;
            IncompressibleSolverControls c;
            c.algorithm = PressureVelocityAlgorithm::PISO;
            c.kinematic_viscosity = 0.1;
            c.coupling.alpha_p = 1.0;
            c.coupling.n_pressure_correctors = ncorrectors;
            c.coupling.piso_neighbor_correction = correction;
            c.convergence.max_iterations = 1;
            c.linear_tolerance = 1e-11;
            c.pressure_reference_cell = 0;
            const auto result = solve_steady_incompressible(m,U,p,ubc,pbc,c);
            EXPECT_TRUE(result.history.size() == 1);
            EXPECT_TRUE(std::isfinite(result.history.back().continuity_linf));
            EXPECT_TRUE(std::isfinite(result.history.back().momentum_equation_residual));
            return std::make_tuple(result.history.back(), U, p);
        };

        const auto [one_h,one_U,one_p] = run(false,1);
        const auto [one_enabled_h,one_enabled_U,one_enabled_p] = run(true,1);
        // With one corrector, the option cannot change the current solve.
        EXPECT_TRUE(one_h.piso_neighbor_predictor_rebuilds == 0);
        EXPECT_TRUE(one_enabled_h.piso_neighbor_predictor_rebuilds == 0);
        for (std::size_t c=0; c<2; ++c) {
            EXPECT_NEAR(one_U(c,0),one_enabled_U(c,0),1e-13);
            EXPECT_NEAR(one_p(c),one_enabled_p(c),1e-13);
        }
        const auto [legacy_h,legacy_U,legacy_p] = run(false,3);
        const auto [improved_h,improved_U,improved_p] = run(true,3);
        (void)legacy_U; (void)legacy_p; (void)improved_U; (void)improved_p;
        EXPECT_TRUE(legacy_h.pressure_correctors_used == 3);
        EXPECT_TRUE(legacy_h.piso_neighbor_predictor_rebuilds == 0);
        EXPECT_TRUE(improved_h.pressure_correctors_used == 3);
        EXPECT_TRUE(improved_h.piso_neighbor_predictor_rebuilds == 2);
        EXPECT_TRUE(improved_h.piso_neighbor_predictor_flux_delta_linf > 1e-12);
        // Reconstructing the pressure flux from HbyA/p must agree with the
        // authoritative mass flux; an inconsistency would mask mass defects.
        EXPECT_TRUE(std::isfinite(improved_h.rhie_chow_identity_linf));
        EXPECT_TRUE(improved_h.rhie_chow_identity_linf < 1e-10);
    });

    run_case("piso_adaptive_acceleration_reduces_inner_work", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0); p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::PISO;
        c.convergence.max_iterations = 5;
        c.convergence.continuity_tolerance = 1e-12;
        c.acceleration.adaptive_linear_tolerance = true;
        c.acceleration.adaptive_pressure_correctors = true;
        c.acceleration.pressure_correctors.min_correctors = 1;
        c.acceleration.pressure_correctors.max_correctors = 3;

        const auto r = solve_steady_incompressible(m,U,p,ubc,pbc,c);
        EXPECT_TRUE(r.converged);
        EXPECT_TRUE(r.history.size() == 2);
        EXPECT_TRUE(r.history[0].pressure_correctors_used == 3);
        EXPECT_TRUE(r.history[1].pressure_correctors_used == 1);
        EXPECT_NEAR(r.history[0].linear_tolerance_used,
                    c.acceleration.linear_forcing.eta_max,1e-14);
        EXPECT_NEAR(r.history[1].linear_tolerance_used,
                    c.acceleration.linear_forcing.eta_min,1e-14);
        EXPECT_TRUE(r.pressure_linear_context.solves == 4);
    });


    run_case("momentum_assembly_uses_matching_pressure_gradient_component", [] {
        const Mesh m = make_unit_cube();
        const auto geometry = build_fv_geometry(m);
        Field<double,Location::FACE> mass_flux(m.n_faces(),"phi","kg/s",1);
        mass_flux.fill(0.0);
        Field<double,Location::CELL> grad_p(1,"grad_p","Pa/m",3);
        grad_p.set(0,1.0,2.0,3.0);
        Field<double,Location::CELL> body(1,"body","N/m3",1);
        body.fill(0.0);
        ScalarBoundaryConditions ubc;
        ubc["wall"] = {ScalarBoundaryType::FIXED_VALUE,0.0,0.0};

        for (std::size_t component = 0; component < 3; ++component) {
            const auto eq = assemble_momentum_component(
                m, geometry, mass_flux, grad_p, body, 1.0, ubc, component, false);
            EXPECT_NEAR(eq.rhs(0), -static_cast<double>(component + 1), 1e-14);
        }
    });

    run_case("momentum_assembly_consumes_mathematical_constraints", [] {
        const Mesh m = make_unit_cube();
        const auto geometry = build_fv_geometry(m);
        Field<double,Location::FACE> mass_flux(m.n_faces(),"phi","kg/s",1);
        mass_flux.fill(0.0);
        Field<double,Location::CELL> grad_p(1,"grad_p","Pa/m",3);
        grad_p.set(0,1.0,2.0,3.0);
        Field<double,Location::CELL> body(1,"body","N/m3",1);
        body.fill(0.0);
        BoundaryConstraintMap constraints;
        auto zero = std::make_shared<ConstantValueProvider>(0.0);
        constraints["wall"] = {
            {"U.x", Dirichlet{zero}},
            {"U.y", Dirichlet{zero}},
            {"U.z", Dirichlet{zero}}};
        for (std::size_t component = 0; component < 3; ++component) {
            const auto eq = assemble_momentum_component(
                m, geometry, mass_flux, grad_p, body, 1.0, constraints, component, false);
            EXPECT_NEAR(eq.rhs(0), -static_cast<double>(component + 1), 1e-14);
        }
    });

    run_case("rhie_chow_rejects_nonfinite_inverse_diagonal", [] {
        const Mesh m = make_unit_cube();
        const auto geometry = build_fv_geometry(m);
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        EXPECT_THROW(
            make_rhie_chow_mass_flux(m, geometry, U, p,
                                      {std::numeric_limits<double>::quiet_NaN()},
                                      1.0, ubc),
            std::invalid_argument);
    });

    run_case("fixed_pressure_boundary_is_accepted_by_coupled_solver", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(100.0);

        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::FIXED_VALUE,100.0,0.0};

        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::SIMPLE;
        c.convergence.max_iterations = 5;
        c.convergence.continuity_tolerance = 1e-12;
        c.linear_tolerance = 1e-12;
        c.pressure_reference_cell = 0;
        c.pressure_reference_value = 100.0;

        const auto r = solve_steady_incompressible(m,U,p,ubc,pbc,c);
        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(p(0),100.0,1e-12);
    });


    run_case("dat_restart_roundtrip_thermal_and_turbulence_fields", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3), p(1,"p","Pa",1);
        Field<double,Location::CELL> T(1,"T","K",1), k(1,"k","m2/s2",1), second(1,"omega","1/s",1);
        U.set(0,1.0,2.0,3.0); p(0)=42.0; T(0)=315.0; k(0)=0.12; second(0)=4.5;
        const auto dat = std::filesystem::temp_directory_path() / "cfdx_thermal_turbulence_restart.dat";
        cfdx::io::DatRestartFields out_fields{&T,&k,&second};
        cfdx::io::write_dat_restart_fields(dat.string(),m,U,p,out_fields,23,1.75);

        Field<double,Location::CELL> loaded_u(1,"U","m/s",3), loaded_p(1,"p","Pa",1);
        Field<double,Location::CELL> loaded_T(1,"T","K",1), loaded_k(1,"k","m2/s2",1), loaded_second(1,"omega","1/s",1);
        cfdx::io::DatRestartFields in_fields{&loaded_T,&loaded_k,&loaded_second};
        const auto state=cfdx::io::read_dat_restart_fields(dat.string(),m,loaded_u,loaded_p,in_fields);
        EXPECT_TRUE(state.iteration==23);
        EXPECT_NEAR(state.time,1.75,1e-14);
        EXPECT_NEAR(loaded_T(0),315.0,1e-14);
        EXPECT_NEAR(loaded_k(0),0.12,1e-14);
        EXPECT_NEAR(loaded_second(0),4.5,1e-14);
        std::filesystem::remove(dat);
    });

    run_case("nonlinear_retry_controls_are_validated_by_the_solver", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0); p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::SIMPLE;
        c.convergence.max_iterations = 2;
        c.convergence.relative_tolerance = 1e-8;
        c.convergence.continuity_tolerance = 1e-8;
        c.linear_tolerance = 1e-10;
        c.pressure_reference_cell = 0;
        c.pressure_reference_value = 0.0;

        // The solver owns the nonlinear retry policy, so an invalid policy must
        // be rejected at the solver boundary and not only inside the controller.
        bool rejected = false;
        try {
            c.nonlinear_retry.max_retries = 0;
            solve_steady_incompressible(m, U, p, ubc, pbc, c);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        EXPECT_TRUE(rejected);

        rejected = false;
        try {
            c.nonlinear_retry.max_retries = 3;
            c.nonlinear_retry.relaxation_shrink = 1.5;
            solve_steady_incompressible(m, U, p, ubc, pbc, c);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        EXPECT_TRUE(rejected);
    });

    run_case("nonlinear_state_rollback_restores_fields_exactly", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0); p.fill(0.0);

        NonlinearStateRollback transaction(U, p);
        transaction.begin();
        U.component_data(0)[0] = 3.5;
        U.component_data(1)[0] = -1.25;
        p(0) = 42.0;
        transaction.reject();

        EXPECT_TRUE(U.component_data(0)[0] == 0.0);
        EXPECT_TRUE(U.component_data(1)[0] == 0.0);
        EXPECT_TRUE(p(0) == 0.0);
        EXPECT_TRUE(!transaction.active());
    });

    run_case("convergence_qoi_gate_is_enforced_by_the_solver", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0); p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::SIMPLE;
        c.convergence.max_iterations = 40;
        c.convergence.relative_tolerance = 1e-8;
        c.convergence.continuity_tolerance = 1e-8;
        c.linear_tolerance = 1e-10;
        c.pressure_reference_cell = 0;
        c.pressure_reference_value = 0.0;
        c.body_force = cfdx::core::Vec3{0.0,0.0,-1.0};
        // Gate a probe against a target the solve cannot reach. The residual
        // and continuity criteria are satisfiable on this mesh, so if the QoI
        // gate did not participate in the acceptance decision the solve would
        // report convergence and this case would fail.
        c.probes.push_back(IncompressiblePointProbe{
            "centre_p", cfdx::core::Vec3{0.5,0.5,0.5},
            IncompressibleProbeField::PRESSURE});
        c.qoi_gates.push_back(cfdx::core::QoIGate{
            "centre_p", 1e-6, 1e-9, 1.0e6, true});

        const auto result = solve_steady_incompressible(m, U, p, ubc, pbc, c);
        EXPECT_TRUE(!result.converged);
    });

    run_case("convergence_qoi_gate_must_name_a_configured_probe", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        U.fill(0.0); p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls c;
        c.algorithm = PressureVelocityAlgorithm::SIMPLE;
        c.convergence.max_iterations = 2;
        c.convergence.relative_tolerance = 1e-8;
        c.convergence.continuity_tolerance = 1e-8;
        c.linear_tolerance = 1e-10;
        c.pressure_reference_cell = 0;
        c.pressure_reference_value = 0.0;

        // A gate naming nothing would be silently unsatisfiable forever, so it
        // must be rejected instead of degrading the solve to never-converged.
        c.qoi_gates.push_back(cfdx::core::QoIGate{
            "absent_probe", 1e-12, 1e-14, 0.0, false});
        bool rejected = false;
        try {
            solve_steady_incompressible(m, U, p, ubc, pbc, c);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        EXPECT_TRUE(rejected);
    });


    run_case("momentum_assembly_uses_harmonic_cellwise_viscosity", [] {
        const Mesh m = make_two_cell_channel();
        const FvGeometry geometry = build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> grad_p(m.n_cells(),"grad_p","Pa/m",3);
        Field<double,Location::CELL> body(m.n_cells(),"body","N/m3",1);
        flux.fill(0.0);
        grad_p.fill(0.0);
        body.fill(0.0);
        ScalarBoundaryConditions bc;
        bc["wall"] = {ScalarBoundaryType::FIXED_VALUE,0.0,0.0};

        const std::vector<double> mu_cell{1.0,3.0};
        const auto variable = assemble_momentum_component(
            m,geometry,flux,grad_p,body,mu_cell,bc,0,true);
        // The internal face uses the harmonic coefficient 1.5 and has
        // unit area and unit centre distance on this two-cell mesh.
        EXPECT_NEAR(variable.matrix(0,1),-1.5,1e-12);
        EXPECT_NEAR(variable.matrix(1,0),-1.5,1e-12);

        const std::vector<double> uniform_mu{2.0,2.0};
        const auto variable_uniform = assemble_momentum_component(
            m,geometry,flux,grad_p,body,uniform_mu,bc,0,true);
        const auto scalar_uniform = assemble_momentum_component(
            m,geometry,flux,grad_p,body,2.0,bc,0,true);
        for (std::size_t row=0; row<m.n_cells(); ++row) {
            for (std::size_t col=0; col<m.n_cells(); ++col)
                EXPECT_NEAR(variable_uniform.matrix(row,col),
                            scalar_uniform.matrix(row,col),1e-12);
        }

        const std::vector<double> wrong_size{1.0};
        EXPECT_THROW(assemble_momentum_component(
            m,geometry,flux,grad_p,body,wrong_size,bc,0,true),
            std::invalid_argument);
        const std::vector<double> invalid_mu{
            1.0,std::numeric_limits<double>::quiet_NaN()};
        EXPECT_THROW(assemble_momentum_component(
            m,geometry,flux,grad_p,body,invalid_mu,bc,0,true),
            std::invalid_argument);
    });

    run_case("cellwise_turbulent_viscosity_controls_are_validated", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3), p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls controls;
        controls.cell_turbulent_viscosity = {1.0e-3, 2.0e-3};
        EXPECT_THROW(
            solve_steady_incompressible(m, U, p, ubc, pbc, controls),
            std::invalid_argument);

        controls.cell_turbulent_viscosity = {-1.0e-3};
        EXPECT_THROW(
            solve_steady_incompressible(m, U, p, ubc, pbc, controls),
            std::invalid_argument);

        controls.cell_turbulent_viscosity = {
            std::numeric_limits<double>::quiet_NaN()};
        EXPECT_THROW(
            solve_steady_incompressible(m, U, p, ubc, pbc, controls),
            std::invalid_argument);

        controls.cell_turbulent_viscosity = {1.0e-3};
        controls.turbulent_viscosity = 1.0e-3;
        EXPECT_THROW(
            solve_steady_incompressible(m, U, p, ubc, pbc, controls),
            std::invalid_argument);

        controls.cell_turbulent_viscosity.clear();
        controls.turbulent_viscosity = 0.0;
        controls.density = 1.0e308;
        controls.kinematic_viscosity = 1.0e308;
        controls.initialization.mode = InitializationMode::Uniform;
        controls.initialization.uniform_velocity = {2.0, 0.0, 0.0};
        U.fill(0.0);
        EXPECT_THROW(
            solve_steady_incompressible(m, U, p, ubc, pbc, controls),
            std::invalid_argument);
        EXPECT_NEAR(U.component_data(0)[0], 0.0, 0.0);

        // Finite scalar inputs can still overflow only after adding the
        // cellwise turbulent contribution. Reject that derived value before
        // the uniform initialization mutates the caller-owned velocity.
        controls.density = 1.0e308;
        controls.kinematic_viscosity = 1.0;
        controls.turbulent_viscosity = 0.0;
        controls.cell_turbulent_viscosity = {1.0e308};
        U.fill(0.0);
        EXPECT_THROW(
            solve_steady_incompressible(m, U, p, ubc, pbc, controls),
            std::invalid_argument);
        EXPECT_NEAR(U.component_data(0)[0], 0.0, 0.0);
    });

    run_case("cellwise_effective_viscosity_overflow_is_rejected_before_initialization", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::CELL> U(1,"U","m/s",3), p(1,"p","Pa",1);
        U.fill(0.0);
        p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{0.0,0.0,0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        IncompressibleSolverControls controls;
        controls.density = 1.0e308;
        controls.kinematic_viscosity = 0.0;
        controls.cell_turbulent_viscosity = {1.0e308};
        controls.initialization.mode = InitializationMode::Uniform;
        controls.initialization.uniform_velocity = {2.0, 0.0, 0.0};

        EXPECT_THROW(
            solve_steady_incompressible(m, U, p, ubc, pbc, controls),
            std::invalid_argument);
        EXPECT_NEAR(U.component_data(0)[0], 0.0, 0.0);
    });

    // Closes the probe-observer contract: a sample produced by an iteration that
    // the nonlinear retry controller rejects must never be published, and the
    // committed iteration's samples must reach probe_callback exactly once.
    run_case("probe_samples_publish_only_committed_iteration_across_nonlinear_retry", [] {
        const Mesh m = make_unit_cube();
        Field<double, Location::CELL> U(1, "U", "m/s", 3), p(1, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);
        VelocityBoundaryConditions ubc;
        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0.0, 0.0, 0.0}};
        ScalarBoundaryConditions pbc;
        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        IncompressibleSolverControls controls;
        controls.convergence.max_iterations = 1;
        controls.pressure_reference_cell = 0;
        controls.nonlinear_retry.max_retries = 1;
        controls.probes = {
            {"retry_probe", {0.5, 0.5, 0.5}, IncompressibleProbeField::PRESSURE},
        };
        std::vector<IncompressibleProbeSample> published_probe_samples;
        controls.probe_callback = [&](const IncompressibleProbeSample& sample) {
            published_probe_samples.push_back(sample);
        };
        std::size_t output_calls = 0;
        controls.iteration_output_callback =
            [&](std::size_t, double, const auto&, const auto&, const auto&) {
                if (++output_calls == 1)
                    throw NonlinearRetryableFailure("deliberate probe-publication rejection");
                return false;
            };
        const auto accepted = solve_steady_incompressible(m, U, p, ubc, pbc, controls);
        EXPECT_TRUE(output_calls == 2);
        EXPECT_TRUE(accepted.history.size() == 1);
        EXPECT_TRUE(accepted.probe_samples.size() == 1);
        EXPECT_TRUE(published_probe_samples.size() == accepted.probe_samples.size());
        EXPECT_TRUE(published_probe_samples[0].iteration == accepted.probe_samples[0].iteration);
        EXPECT_NEAR(published_probe_samples[0].value, accepted.probe_samples[0].value, 1e-14);
    });

    return run_all();
}
