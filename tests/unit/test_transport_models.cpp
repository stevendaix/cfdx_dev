// M0.14-T02 — Transport Models (Viscosity, thermal conductivity, diffusivity)
// Unit tests for transport property models
#include "cfdx/physics/transport_models.h"
#include "cfdx/physics/equation_of_state.h"
#include "cfdx/core/field/field.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("eos_rejects_invalid_parameters", []() {
        IdealGasParams ideal;
        ideal.gamma = 0.5;
        EXPECT_THROW(IdealGasEOS(ideal), std::invalid_argument);
        IncompressibleParams inc;
        inc.rho = -1.0;
        EXPECT_THROW(IncompressibleEOS(inc), std::invalid_argument);
    });

    run_case("ideal_gas_density", []() {
        IdealGasParams params;
        params.M = 0.02896546;
        params.gamma = 1.4;
        params.T_ref = 300.0;
        params.p_ref = 101325.0;
        
        IdealGasEOS eos(params);
        double rho = eos.density(101325.0, 300.0);
        double R = 8.314462618 / 0.02896546;
        double expected = 101325.0 / (R * 300.0);
        EXPECT_NEAR(rho, expected, 1e-6);
        EXPECT_NEAR(rho, 1.167, 0.1);
    });

    run_case("ideal_gas_enthalpy", []() {
        IdealGasParams params;
        params.M = 0.02896546;
        params.gamma = 1.4;
        params.T_ref = 300.0;
        
        IdealGasEOS eos(params);
        double h = eos.enthalpy(101325.0, 300.0);
        EXPECT_NEAR(h, 0.0, 1e-6);
    });

    run_case("ideal_gas_speed_of_sound", []() {
        IdealGasParams params;
        params.M = 0.02896546;
        params.gamma = 1.4;
        params.T_ref = 300.0;
        
        IdealGasEOS eos(params);
        double cs = eos.speed_of_sound(101325.0, 300.0);
        double R = 8.314462618 / 0.02896546;
        double expected = std::sqrt(1.4 * R * 300.0);
        EXPECT_NEAR(cs, expected, 0.1);
        EXPECT_TRUE(cs > 0);
    });

    run_case("ideal_gas_density_comparison", []() {
        IdealGasParams params;
        params.M = 0.02896546;
        params.gamma = 1.4;
        params.T_ref = 300.0;
        params.p_ref = 101325.0;
        
        IdealGasEOS eos(params);
        double rho = eos.density(101325.0, 300.0);
        EXPECT_NEAR(rho, 1.167, 0.1);
    });

    run_case("constant_viscosity", []() {
        double mu = 1.8e-5;
        double mu_val = constant_viscosity(mu);
        EXPECT_NEAR(mu_val, mu, 1e-15);
    });

    run_case("sutherland_viscosity", []() {
        double T = 300.0;
        double mu = sutherland_viscosity(T);
        double expected = 1.716e-5 * std::pow(T / 273.15, 1.5) * (273.15 + 110.4) / (T + 110.4);
        EXPECT_NEAR(mu, expected, 1e-12);
        EXPECT_TRUE(mu > 0);
    });

    run_case("power_law_viscosity", []() {
        double mu0 = 1.8e-5;
        double T0 = 273.15;
        double n = 0.7;
        double mu_val = power_law_viscosity(T0, mu0, T0, n);
        EXPECT_NEAR(mu_val, mu0, 1e-15);
    });

    run_case("prandtl_conductivity", []() {
        double mu = 1.8e-5;
        double Cp = 1004.5;
        double Pr = 0.71;
        double k = prandtl_conductivity(mu, Cp, Pr);
        double expected = mu * Cp / Pr;
        EXPECT_NEAR(k, expected, 1e-10);
    });

    run_case("compute_transport_with_rho", []() {
        double T = 300.0;
        double rho = 1.2;
        TransportProperties tp = compute_transport(T, rho);
        EXPECT_TRUE(tp.mu > 0);
        EXPECT_TRUE(tp.k > 0);
        EXPECT_TRUE(tp.D > 0);
        EXPECT_TRUE(tp.Cp > 0);
    });

    run_case("compute_transport_with_eos", []() {
        IdealGasEOS eos;
        eos.set_params(IdealGasParams{0.02896546, 1.4, 300.0, 101325.0});
        double T = 300.0;
        double p = 101325.0;
        TransportProperties tp = compute_transport(T, eos, p);
        EXPECT_TRUE(tp.mu > 0);
        EXPECT_TRUE(tp.k > 0);
        EXPECT_TRUE(tp.D > 0);
        EXPECT_TRUE(tp.rho > 0);
    });

    run_case("incompressible_eos", []() {
        IncompressibleEOS eos;
        IncompressibleParams params;
        params.rho = 1000.0;
        params.Cp = 4182.0;
        eos.set_params(params);
        
        double rho = eos.density(101325.0, 293.15);
        EXPECT_NEAR(rho, 1000.0, 1e-6);
        
        double h = eos.enthalpy(101325.0, 300.0);
        EXPECT_NEAR(h, 4182.0 * (300.0 - 300.0), 1e-3);
    });

    run_case("transport_fields_on_field", []() {
        Field<double, Location::CELL> T(3, "T", "K", 1);
        T(0) = 280.0; T(1) = 290.0; T(2) = 310.0;
        
        Field<double, Location::CELL> mu(3, "mu", "Pa.s", 1);
        Field<double, Location::CELL> k(3, "k", "W/m/K", 1);
        Field<double, Location::CELL> D(3, "D", "m2/s", 1);
        
        double rho = 1.2;
        compute_transport_fields(T, mu, k, D, rho);
        
        EXPECT_TRUE(mu(0) > 0);
        EXPECT_TRUE(k(0) > 0);
        EXPECT_TRUE(D(0) > 0);
    });

    return run_all();
}