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

// Four closed hexahedra: dx=1/4, dy=dz=1/10, V=1/400 each.
// The former coplanar fixture had zero volume and only exercised the 1e30
// fallback in the CFL/local-step routines.
Mesh make_test_mesh() {
    constexpr std::size_t cells=4;
    Mesh mesh;
    mesh.points().resize(4*(cells+1));
    for(std::size_t i=0;i<=cells;++i) {
        const auto b=4*i;
        const double x=.25*static_cast<double>(i);
        mesh.points().set(b,  x,0.0,0.0);
        mesh.points().set(b+1,x,.1,0.0);
        mesh.points().set(b+2,x,.1,.1);
        mesh.points().set(b+3,x,0.0,.1);
        if(i==0) mesh.faces().push_face({b,b+3,b+2,b+1});
        else mesh.faces().push_face({b,b+1,b+2,b+3});
    }
    for(std::size_t i=0;i<cells;++i) {
        const auto b=4*i;
        const auto side=mesh.n_faces();
        mesh.faces().push_face({b,b+4,b+7,b+3});
        mesh.faces().push_face({b+1,b+2,b+6,b+5});
        mesh.faces().push_face({b,b+1,b+5,b+4});
        mesh.faces().push_face({b+3,b+7,b+6,b+2});
        mesh.cells().push_cell({i,i+1,side,side+1,side+2,side+3});
    }
    mesh.ownership().resize(mesh.n_faces());
    for(std::size_t i=0;i<=cells;++i) {
        mesh.ownership().set_owner(i,i==0 ? 0 : i-1);
        mesh.ownership().set_neighbour(i,i>0 && i<cells
            ? static_cast<CellIndex>(i) : FaceOwnership::BOUNDARY);
    }
    for(std::size_t i=0;i<cells;++i)
        for(std::size_t side=0;side<4;++side) {
            const auto face=cells+1+4*i+side;
            mesh.ownership().set_owner(face,i);
            mesh.ownership().set_neighbour(face,FaceOwnership::BOUNDARY);
        }
    BoundaryPatches patches;
    patches.add_patch("inlet",PatchType::INLET);
    patches.patch(0).face_ids={0};
    patches.add_patch("outlet",PatchType::OUTLET);
    patches.patch(1).face_ids={cells};
    patches.add_patch("walls",PatchType::WALL);
    for(std::size_t face=cells+1;face<mesh.n_faces();++face)
        patches.patch(2).face_ids.push_back(face);
    mesh.set_boundary(patches);
    return mesh;
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
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(1.0);

        double dt = 0.1;
        auto phi_new = advance_time(phi, dt, decay_rhs, TimeScheme::EULER_EXPLICIT);

        // Exact: φ(t) = φ0 * exp(-t)
        // Euler: φ^{n+1} = φ^n * (1 - dt)
        EXPECT_NEAR(phi_new(0), 1.0 * (1.0 - dt), 1e-12);
    });

    run_case("euler_implicit_decay", []() {
        Field<double,Location::CELL> phi(4,"phi","1",2);
        for(std::size_t cell=0;cell<phi.size();++cell) {
            phi(cell,0)=1.0+.1*static_cast<double>(cell);
            phi(cell,1)=.5+.2*static_cast<double>(cell);
        }
        constexpr double dt=.1;
        const auto result=advance_time(phi,dt,decay_rhs,TimeScheme::EULER_IMPLICIT);
        for(std::size_t cell=0;cell<phi.size();++cell)
            for(std::size_t component=0;component<phi.dimension();++component)
                EXPECT_NEAR(result(cell,component),phi(cell,component)/(1.0+dt),1e-12);
    });

    run_case("crank_nicolson_decay", [&]() {
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(1.0);

        double dt = 0.1;
        auto phi_new = advance_time(phi, dt, decay_rhs, TimeScheme::CRANK_NICOLSON);

        // Solve (1+dt/2) phi^{n+1} = (1-dt/2) phi^n.
        // The former oracle was Heun's explicit predictor-corrector, not CN.
        const double expected=(1.0-.5*dt)/(1.0+.5*dt);
        for(std::size_t cell=0;cell<phi.size();++cell) {
            EXPECT_NEAR(phi_new(cell),expected,1e-12);
            EXPECT_NEAR(phi_new(cell)-phi(cell)+.5*dt*(phi(cell)+phi_new(cell)),0.0,1e-12);
        }
    });

    run_case("bdf2_requires_history", [&]() {
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(1.0);

        TimeIntegrationContext ctx(4, 1, "phi");
        ctx.initialize(phi);

        double dt = 0.1;
        // Bootstrap solves (1+dt) phi1 = phi0 by implicit Euler.
        EXPECT_FALSE(ctx.has_prev);
        const auto phi1=advance_time(phi,dt,decay_rhs,TimeScheme::BDF2,&ctx);
        EXPECT_TRUE(ctx.has_prev);
        EXPECT_NEAR(ctx.dt_prev,dt,1e-15);
        for(std::size_t cell=0;cell<phi.size();++cell) {
            EXPECT_NEAR(phi1(cell),1.0/(1.0+dt),1e-12);
            EXPECT_NEAR(ctx.phi_prev(cell),phi(cell),1e-15);
            EXPECT_NEAR(ctx.phi_curr(cell),phi1(cell),1e-15);
        }

        // BDF2: (3+2dt) phi2 = 4 phi1 - phi0.
        const auto phi2=advance_time(phi1,dt,decay_rhs,TimeScheme::BDF2,&ctx);
        for(std::size_t cell=0;cell<phi.size();++cell) {
            const double expected=(4.0*phi1(cell)-phi(cell))/(3.0+2.0*dt);
            EXPECT_NEAR(phi2(cell),expected,1e-12);
            EXPECT_NEAR((3.0+2.0*dt)*phi2(cell)-4.0*phi1(cell)+phi(cell),0.0,1e-12);
            EXPECT_NEAR(ctx.phi_prev(cell),phi1(cell),1e-15);
            EXPECT_NEAR(ctx.phi_curr(cell),phi2(cell),1e-15);
        }
    });

    run_case("bdf2_variable_step_implicit_decay", []() {
        Field<double,Location::CELL> phi(4,"phi","1",1);
        phi.fill(1.0);
        TimeIntegrationContext ctx(4,1,"phi");
        ctx.initialize(phi);
        const auto first=advance_time(phi,.1,decay_rhs,TimeScheme::BDF2,&ctx);
        // dt_new / dt_old = 2: 5*phi2 - 9*phi1 + 4*phi0 = .6*RHS(phi2).
        // This closed-form oracle does not call bdf2_coefficients().
        const auto second=advance_time(first,.2,decay_rhs,TimeScheme::BDF2,&ctx);
        for(std::size_t cell=0;cell<phi.size();++cell)
            EXPECT_NEAR(second(cell),(9.0*first(cell)-4.0*phi(cell))/5.6,1e-12);
        EXPECT_NEAR(ctx.dt_prev,.2,1e-15);
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

        // The existing API uses h=V^(1/3), with V=.25*.1*.1.
        const double expected=.5*std::cbrt(.0025)/10.0;
        EXPECT_NEAR(dt,expected,1e-12);
        EXPECT_NEAR(compute_cfl_time_step(mesh,velocity,.5,2.0),
                    .5*std::cbrt(.0025)/12.0,1e-12);
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

        // Equal positive volumes make the expected local dt independent of
        // the geometry routine under test. A 1e30 fallback must fail.
        for(std::size_t c=0;c<4;++c)
            EXPECT_NEAR(dt_local[c],.5*std::cbrt(.0025)/velocity(c,0),1e-12);

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