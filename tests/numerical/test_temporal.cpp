// M0.13-T02/T03 — Tests for Temporal Discretization (Local Time Stepping, Adaptive Time Step)

#include "cfdx/core/numerics/temporal.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/boundary.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

// Helper: build a simple 1D mesh with 4 cells (thin 3D slab)
Mesh make_test_mesh() {
    Mesh m;

    m.points().resize(10);
    m.points().set(0, 0.0, 0.0, 0.0);
    m.points().set(1, 0.25, 0.0, 0.0);
    m.points().set(2, 0.5, 0.0, 0.0);
    m.points().set(3, 0.75, 0.0, 0.0);
    m.points().set(4, 1.0, 0.0, 0.0);
    m.points().set(5, 0.0, 0.0, 0.1);
    m.points().set(6, 0.25, 0.0, 0.1);
    m.points().set(7, 0.5, 0.0, 0.1);
    m.points().set(8, 0.75, 0.0, 0.1);
    m.points().set(9, 1.0, 0.0, 0.1);

    m.faces().push_face({0, 5, 6, 1});
    m.faces().push_face({1, 6, 7, 2});
    m.faces().push_face({2, 7, 8, 3});
    m.faces().push_face({3, 8, 9, 4});
    m.faces().push_face({0, 1, 6, 5});
    m.faces().push_face({4, 9, 8, 3});

    m.ownership().resize(6);
    m.ownership().set_owner(0, 0); m.ownership().set_neighbour(0, 1);
    m.ownership().set_owner(1, 1); m.ownership().set_neighbour(1, 2);
    m.ownership().set_owner(2, 2); m.ownership().set_neighbour(2, 3);
    m.ownership().set_owner(3, 3); m.ownership().set_neighbour(3, FaceOwnership::BOUNDARY);
    m.ownership().set_owner(4, 0); m.ownership().set_neighbour(4, FaceOwnership::BOUNDARY);
    m.ownership().set_owner(5, 3); m.ownership().set_neighbour(5, FaceOwnership::BOUNDARY);

    m.cells().push_cell({0, 4});
    m.cells().push_cell({0, 1});
    m.cells().push_cell({1, 2});
    m.cells().push_cell({2, 3, 5});

    BoundaryPatches bp;
    bp.add_patch("inlet", PatchType::INLET);
    bp.patch(0).face_ids = {4};
    bp.add_patch("outlet", PatchType::OUTLET);
    bp.patch(1).face_ids = {5};
    m.set_boundary(bp);

    return m;
}

// Simple RHS for testing: dφ/dt = -φ (decay)
void decay_rhs(const Field<double, Location::CELL>& phi, Field<double, Location::CELL>& rhs) {
    const std::size_t n = phi.size();
    const std::size_t dim = phi.dimension();
    for (std::size_t c = 0; c < n; ++c) {
        for (std::size_t d = 0; d < dim; ++d) {
            rhs.component_data(d)[c] = -phi.component_data(d)[c];
        }
    }
}

int main() {
    run_case("euler_explicit_decay", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(1.0);

        double dt = 0.1;
        auto phi_new = advance_time(TimeScheme::EULER_EXPLICIT, phi, dt, decay_rhs);

        // Exact: φ(t) = φ0 * exp(-t)
        // Euler: φ^{n+1} = φ^n * (1 - dt)
        EXPECT_NEAR(phi_new(0), 1.0 * (1.0 - dt), 1e-12);
    });

    run_case("crank_nicolson_decay", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(1.0);

        double dt = 0.1;
        auto phi_new = advance_time(TimeScheme::CRANK_NICOLSON, phi, dt, decay_rhs);

        // CN predictor-corrector:
        // Predictor: φ* = φ^n + Δt * RHS(φ^n) = 1.0 + 0.1 * (-1.0) = 0.9
        // Corrector: φ^{n+1} = φ^n + 0.5*Δt * (RHS(φ^n) + RHS(φ*))
        // = 1.0 + 0.5*0.1 * (-1.0 + -0.9) = 1.0 - 0.05 * 1.9 = 1.0 - 0.095 = 0.905
        double phi_star = 1.0 + dt * (-1.0);
        double expected = 1.0 + 0.5 * dt * (-1.0 + -phi_star);
        EXPECT_NEAR(phi_new(0), expected, 1e-12);
    });

    run_case("bdf2_requires_history", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(1.0);

        TimeIntegrationContext ctx(4, 1, "phi");
        ctx.initialize(phi);

        double dt = 0.1;
        // First step (falls back to Euler)
        auto phi1 = advance_time(TimeScheme::BDF2, phi, dt, decay_rhs, &ctx);
        EXPECT_NEAR(phi1(0), 1.0 * (1.0 - dt), 1e-12);

        // Second step (uses BDF2 explicit):
        // φ2 = (4φ1 - φ0 + 2Δt * RHS(φ1)) / 3
        // RHS(φ1) = -φ1
        auto phi2 = advance_time(TimeScheme::BDF2, phi1, dt, decay_rhs, &ctx);
        double expected = (4.0 * phi1(0) - 1.0 + 2.0 * dt * (-phi1(0))) / 3.0;
        EXPECT_NEAR(phi2(0), expected, 1e-12);
    });

    run_case("compute_cfl_time_step", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> velocity(4, "U", "m/s", 3);
        // Set uniform velocity of 10 m/s in x direction
        for (std::size_t c = 0; c < 4; ++c) {
            velocity(c, 0) = 10.0;  // ux
            velocity(c, 1) = 0.0;   // uy
            velocity(c, 2) = 0.0;   // uz
        }

        double dt = compute_cfl_time_step(mesh, velocity, 0.5, 0.0);

        // dt should be positive and finite (mesh dependent)
        EXPECT_TRUE(dt > 0.0);
        EXPECT_TRUE(std::isfinite(dt));
    });

    run_case("adaptive_time_stepper", [&]() {
        AdaptiveTimeStepper stepper;
        stepper.dt = 0.1;
        stepper.dt_min = 1e-6;
        stepper.dt_max = 1.0;
        stepper.cfl_target = 0.5;

        // If CFL achieved is higher than target, reduce dt
        stepper.update(0.8);  // CFL = 0.8 > 0.5
        EXPECT_TRUE(stepper.dt < 0.1);

        // If CFL achieved is lower than target, increase dt
        stepper.dt = 0.1;
        stepper.update(0.2);  // CFL = 0.2 < 0.5
        EXPECT_TRUE(stepper.dt > 0.1);

        // Test reject
        stepper.dt = 0.1;
        stepper.reject();
        EXPECT_NEAR(stepper.dt, 0.05, 1e-12);  // shrink_factor = 0.5
    });

    run_case("local_time_stepping", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> velocity(4, "U", "m/s", 3);
        // Varying velocities
        velocity(0, 0) = 10.0; velocity(0, 1) = 0.0; velocity(0, 2) = 0.0;
        velocity(1, 0) = 20.0; velocity(1, 1) = 0.0; velocity(1, 2) = 0.0;
        velocity(2, 0) = 5.0;  velocity(2, 1) = 0.0; velocity(2, 2) = 0.0;
        velocity(3, 0) = 1.0;  velocity(3, 1) = 0.0; velocity(3, 2) = 0.0;

        auto dt_local = compute_local_time_steps(mesh, velocity, 0.5);

        // All dt should be positive and finite
        for (std::size_t c = 0; c < 4; ++c) {
            EXPECT_TRUE(dt_local[c] > 0.0);
            EXPECT_TRUE(std::isfinite(dt_local[c]));
        }

        // Apply local time stepping
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(1.0);
        Field<double, Location::CELL> phi_new(4, "phi_new", "1", 1);

        apply_local_time_stepping(phi, dt_local, decay_rhs, phi_new);

        // Each cell updated with its own dt: φ_new = φ + dt_local * (-φ)
        for (std::size_t c = 0; c < 4; ++c) {
            double expected = 1.0 * (1.0 - dt_local[c]);
            EXPECT_NEAR(phi_new(c), expected, 1e-12);
        }
    });

    run_case("time_scheme_enum", [&]() {
        EXPECT_TRUE(std::string(to_string(TimeScheme::EULER_EXPLICIT)) == "euler_explicit");
        EXPECT_TRUE(std::string(to_string(TimeScheme::EULER_IMPLICIT)) == "euler_implicit");
        EXPECT_TRUE(std::string(to_string(TimeScheme::CRANK_NICOLSON)) == "crank_nicolson");
        EXPECT_TRUE(std::string(to_string(TimeScheme::BDF2)) == "bdf2");

        EXPECT_TRUE(time_scheme_from_string("euler_explicit") == TimeScheme::EULER_EXPLICIT);
        EXPECT_TRUE(time_scheme_from_string("euler_implicit") == TimeScheme::EULER_IMPLICIT);
        EXPECT_TRUE(time_scheme_from_string("crank_nicolson") == TimeScheme::CRANK_NICOLSON);
        EXPECT_TRUE(time_scheme_from_string("bdf2") == TimeScheme::BDF2);
    });

    return run_all();
}