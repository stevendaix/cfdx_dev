#include "cfdx/physics/transport_models.h"
#include "common/test_harness.h"

#include <limits>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Field<double, Location::CELL> scalar(std::size_t n, double value)
{
    Field<double, Location::CELL> out(n);
    out.fill(value);
    return out;
}

int main()
{
    run_case("constant_density_invalid_late_cell_is_transactional", [] {
        auto T = scalar(2, 300.0);
        T(1) = -1.0;
        auto mu = scalar(2, 9.0);
        auto k  = scalar(2, 8.0);
        auto D  = scalar(2, 7.0);
        EXPECT_THROW(compute_transport_fields(T, mu, k, D, 1.0), std::invalid_argument);
        EXPECT_NEAR(mu(0), 9.0, 0.0);
        EXPECT_NEAR(k(0), 8.0, 0.0);
        EXPECT_NEAR(D(0), 7.0, 0.0);
    });
    run_case("EOS_invalid_late_cell_is_transactional", [] {
        IdealGasEOS eos;
        auto T = scalar(2, 300.0);
        auto p = scalar(2, 101325.0);
        p(1) = -1.0;
        auto mu = scalar(2, 9.0);
        auto k  = scalar(2, 8.0);
        auto D  = scalar(2, 7.0);
        EXPECT_THROW(compute_transport_fields(T, mu, k, D, eos, p), std::invalid_argument);
        EXPECT_NEAR(mu(0), 9.0, 0.0);
        EXPECT_NEAR(k(0), 8.0, 0.0);
        EXPECT_NEAR(D(0), 7.0, 0.0);
    });
    run_case("EOS_pressure_and_output_shapes_rejected_before_writes", [] {
        IdealGasEOS eos;
        auto T = scalar(2, 300.0);
        auto p = scalar(1, 101325.0);
        auto mu = scalar(2, 9.0);
        auto k  = scalar(2, 8.0);
        auto D  = scalar(2, 7.0);
        EXPECT_THROW(compute_transport_fields(T, mu, k, D, eos, p), std::invalid_argument);
        EXPECT_NEAR(mu(0), 9.0, 0.0);
        auto small = scalar(1, 2.0);
        EXPECT_THROW(compute_transport_fields(T, mu, small, D, 1.0), std::invalid_argument);
        EXPECT_NEAR(mu(0), 9.0, 0.0);
    });
    run_case("wrong_scalar_dimensions_and_aliases_rejected", [] {
        Field<double, Location::CELL> T(2, "T", "K", 2);
        T(0) = T(1) = 300.0;
        auto mu = scalar(2, 9.0);
        auto k = scalar(2, 8.0);
        auto D = scalar(2, 7.0);
        EXPECT_THROW(compute_transport_fields(T, mu, k, D, 1.0), std::invalid_argument);
        auto Tscalar = scalar(2, 300.0);
        EXPECT_THROW(compute_transport_fields(Tscalar, mu, mu, D, 1.0), std::invalid_argument);
        EXPECT_THROW(compute_transport_fields(Tscalar, Tscalar, k, D, 1.0), std::invalid_argument);
    });
    run_case("overflowed_transport_coefficients_do_not_escape", [] {
        EXPECT_THROW(prandtl_conductivity(1.0e308, 1000.0, 0.71), std::runtime_error);
        EXPECT_THROW(schmidt_diffusivity(1.0e308, 1.0e-308, 1.0), std::runtime_error);
    });
    run_case("healthy_transport_properties_unchanged", [] {
        auto T = scalar(2, 300.0);
        T(1) = 350.0;
        auto mu = scalar(2, -1.0);
        auto k  = scalar(2, -1.0);
        auto D  = scalar(2, -1.0);
        compute_transport_fields(T, mu, k, D, 1.2);
        EXPECT_NEAR(mu(0), sutherland_viscosity(300.0), 1e-18);
        EXPECT_NEAR(k(1), prandtl_conductivity(mu(1), 1004.6667184201457, 0.71), 1e-5);
        EXPECT_TRUE(D(1) > 0.0);
    });
    return run_all();
}
