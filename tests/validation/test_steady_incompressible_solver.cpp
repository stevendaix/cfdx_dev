#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/io/restart/dat_restart.h"
#include <filesystem>
#include "common/test_harness.h"
#include <limits>

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

int main()
{
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

        Field<double,Location::CELL> restart_u(1,"U","m/s",3);
        Field<double,Location::CELL> restart_p(1,"p","Pa",1);
        restart_u.fill(0.0);
        restart_p.fill(0.0);
        (void)solve_steady_incompressible(
            m, restart_u, restart_p, ubc, pbc, controls, dat.string());

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

        (void)solve_steady_incompressible(
            m,U,p,ubc,pbc,controls,dat.string(),{&T,&k,&omega});

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
        controls.probe_callback = [&](const IncompressibleProbeSample& sample) {
            samples.push_back(sample);
        };

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

    return run_all();
}
