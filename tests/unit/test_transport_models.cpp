// M0.14-T02 — Transport Models (Viscosity, thermal conductivity, diffusivity)
// Unit tests for transport property models
#include "transport_models.h"
#include <gtest/gtest.h>
#include <memory>

TEST(TEST_TRANSPORT_MODELS, IdealGasEOS), 
    "Ideal gas EOS should compute correct properties"
{
    // Arrange
    IdealGasEOS eos;
    eos.set_params(SutherlandParams{mu0 = 1.716e-5, T0 = 273.15, S = 110.4});
    
    // Act & Assert
    TEST_F(TEST_TRANSPORT_MODELS, ideal_gas_density) {
        double rho = eos.density(101325.0, 300.0);
        EXPECT_FLOAT_EQ(101325.0, rho, 1.0); // Should be close to reference density
    }
    
    TEST_F(TEST_TRANSPORT_MODELS, ideal_gas_enthalpy) {
        double h = eos.enthalpy(101325.0, 300.0);
        EXPECT_FLOAT_EQ(1004.5 * (300.0 - 273.15), h);
    }
    
    TEST_F(TEST_TRANSPORT_MODELS, ideal_gas_speed_of_sound) {
        double cs = eos.speed_of_sound(101325.0, 300.0);
        EXPECT_FLOAT_EQ(340.0, cs, 0.001); // Speed of sound in air ~340 m/s
    }
}

TEST(TEST_TRANSPORT_MODELS, ConstantViscosity), 
    "Constant viscosity model should work"
{
    // Arrange
    double mu = 1.8e-5; // Dynamic viscosity of air
    double T0 = 273.15;
    
    // Act
    double mu_val = constant_viscosity(mu);
    EXPECT_DOUBLE_EQ(mu, mu_val);
}

TEST(TEST_TRANSPORT_MODELS, PowerLawViscosity), 
    "Power law viscosity model should work"
{
    // Arrange
    double mu0 = 1.8e-5;
    double T0 = 273.15;
    double n = 0.7; // Exponent
    
    // Act
    double mu_val = power_law_viscosity(T0, mu0, T0, n);
    
    // Assert - approximate check
    EXPECT_GREATER_THAN(1e-10, mu_val);
}

TEST(TEST_TRANSPORT_PROPERTIES, ComputeTransportFields), 
    "Transport properties should be computed correctly"
{
    // Arrange
    double T = 300.0; // Temperature in K
    SutherlandParams suth;
    suth.mu0 = 1.716e-5;
    suth.T0 = 273.15;
    suth.S = 110.4;
    
    // Act
    double mu = compute_transport(T, suth);
    double k = compute_transport(T, suth, 0.71, 0.7, 1004.5);
    double D = compute_transport(T, suth, 0.71, 0.7, 2.0e-5);
    
    // Assert
    EXPECT_GREATER_THAN(0.0, mu);
    EXPECT_GREATER_THAN(0.0, k);
    EXPECT_GREATER_THAN(0.0, D);
}

TEST(TEST_TRANSPORT_PROPERTIES, ComputeTransportOnField), 
    "Transport properties should be computed on a field"
{
    // Arrange
    double T_field[3] = {280.0, 290.0, 310.0}; // Array of temperatures
    double mu_field[3];
    double k_field[3];
    double D_field[3];
    
    // Act
    compute_transport_fields(T_field, mu_field, k_field, D_field);
    
    // Assert - spot check
    EXPECT_GREATER_THAN(0.0, mu_field[0]);
    EXPECT_GREATER_THAN(0.0, k_field[0]);
    EXPECT_GREATER_THAN(0.0, D_field[0]);
}
