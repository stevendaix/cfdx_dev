#include "cfdx/physics/local_time_stepping.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh cube(double length)
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    for (std::size_t i = 0; i < 8; ++i)
        m.points().set(i, length*p[i][0], length*p[i][1], length*p[i][2]);
    m.faces().push_face({0, 3, 2, 1});
    m.faces().push_face({4, 5, 6, 7});
    m.faces().push_face({0, 1, 5, 4});
    m.faces().push_face({3, 7, 6, 2});
    m.faces().push_face({0, 4, 7, 3});
    m.faces().push_face({1, 2, 6, 5});
    m.ownership().resize(6);
    for (std::size_t i = 0; i < 6; ++i) {
        m.ownership().set_owner(i, 0);
        m.ownership().set_neighbour(i, FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0, 1, 2, 3, 4, 5});
    return m;
}

static double run_cube(double length)
{
    const Mesh m = cube(length);
    Field<double, Location::FACE> flux(6, "mass_flux", "kg/s");
    flux.fill(0.0);
    flux(4) = length * length;  // rho=1, u=1 m/s, area=L^2
    flux(5) = length * length;
    Field<double, Location::CELL> rho(1, "rho", "kg/m3");
    rho(0) = 1.0;
    Field<double, Location::CELL> dt;
    LocalTimeStepControls c;
    c.cfl = 0.5;
    compute_local_time_step(m, flux, rho, dt, c);
    return dt(0);
}

int main()
{
    run_case("volume_scaled_convective_CFL", [] {
        // With two throughflow faces: CFL = dt*2*U/L.
        EXPECT_NEAR(run_cube(1.0), 0.25, 1e-13);
        EXPECT_NEAR(run_cube(2.0), 0.50, 1e-13);
        EXPECT_NEAR(run_cube(0.5), 0.125, 1e-13);
    });
    run_case("explicit_volumes_use_mass_flux_not_volumetric_flux", [] {
        const Mesh m = cube(2.0);
        Field<double, Location::FACE> flux(6); flux.fill(0.0);
        flux(4) = 8.0; flux(5) = 8.0;
        Field<double, Location::CELL> rho(1); rho(0) = 3.0;
        Field<double, Location::CELL> dt;
        LocalTimeStepControls c; c.cfl = 0.5;
        compute_local_time_step(m, flux, rho, std::vector<double>{8.0}, dt, c);
        EXPECT_NEAR(dt(0), 0.75, 1e-13);
    });
    run_case("invalid_volume_does_not_mutate_previous_dt", [] {
        const Mesh m = cube(2.0);
        Field<double, Location::FACE> flux(6); flux.fill(1.0);
        Field<double, Location::CELL> rho(1); rho(0) = 1.0;
        Field<double, Location::CELL> dt(1); dt(0) = 0.4;
        EXPECT_THROW(compute_local_time_step(m, flux, rho,
            std::vector<double>{0.0}, dt), std::invalid_argument);
        EXPECT_NEAR(dt(0), 0.4, 0.0);
        EXPECT_THROW(compute_local_time_step(m, flux, rho,
            std::vector<double>{std::numeric_limits<double>::quiet_NaN()}, dt),
            std::invalid_argument);
        EXPECT_NEAR(dt(0), 0.4, 0.0);
    });
    run_case("invalid_mass_flux_does_not_mutate_previous_dt", [] {
        const Mesh m = cube(1.0);
        Field<double, Location::FACE> flux(6); flux.fill(1.0);
        flux(5) = std::numeric_limits<double>::quiet_NaN();
        Field<double, Location::CELL> rho(1); rho(0) = 1.0;
        Field<double, Location::CELL> dt(1); dt(0) = 0.3;
        EXPECT_THROW(compute_local_time_step(m, flux, rho, dt), std::invalid_argument);
        EXPECT_NEAR(dt(0), 0.3, 0.0);
    });
    run_case("global_pseudo_dt_rejects_nonfinite_and_nonpositive", [] {
        Field<double, Location::CELL> dt(2); dt(0) = 1.0; dt(1) = 0.2;
        EXPECT_NEAR(global_pseudo_time_step(dt), 0.2, 1e-14);
        dt(1) = std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(global_pseudo_time_step(dt), std::invalid_argument);
        dt(1) = -1.0;
        EXPECT_THROW(global_pseudo_time_step(dt), std::invalid_argument);
    });
    return run_all();
}
