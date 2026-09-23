#include "cfdx/physics/thermophysical_models.h"
#include "cfdx/physics/turbulence_models.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("scalar_property_models", [] {
        ScalarPropertyControls c;
        c.reference_value=10.0;
        EXPECT_NEAR(evaluate_scalar_property(c,500.0),10.0,1e-12);
        c.model=ScalarPropertyModel::LINEAR; c.coefficients={0.1};
        EXPECT_NEAR(evaluate_scalar_property(c,400.0),20.0,1e-12);
        c.model=ScalarPropertyModel::POLYNOMIAL; c.coefficients={10.0,0.1,0.001};
        EXPECT_NEAR(evaluate_scalar_property(c,400.0),30.0,1e-12);
        c.model=ScalarPropertyModel::POWER_LAW; c.reference_value=4.0; c.exponent=1.0;
        EXPECT_NEAR(evaluate_scalar_property(c,600.0),8.0,1e-12);
        c.model=ScalarPropertyModel::EXPONENTIAL; c.reference_value=2.0; c.exponent=0.01;
        EXPECT_NEAR(evaluate_scalar_property(c,400.0),2.0*std::exp(0.01*(400.0-300.0)/300.0),1e-12);
        c.model=ScalarPropertyModel::ARRHENIUS; c.reference_value=2.0; c.activation_temperature=1000.0;
        EXPECT_NEAR(evaluate_scalar_property(c,600.0),2.0*std::exp(1000.0*(1.0/300.0-1.0/600.0)),1e-10);
    });
    run_case("tabulated_property_policies", [] {
        ScalarPropertyControls c; c.model=ScalarPropertyModel::TABLE;
        c.table.temperature={300,400,500}; c.table.value={10,20,30};
        EXPECT_NEAR(evaluate_scalar_property(c,450.0),25.0,1e-12);
        EXPECT_NEAR(evaluate_scalar_property(c,200.0),10.0,1e-12);
        c.table.extrapolation=ExtrapolationPolicy::LINEAR;
        EXPECT_NEAR(evaluate_scalar_property(c,250.0),5.0,1e-12);
        c.table.extrapolation=ExtrapolationPolicy::REJECT;
        bool threw=false; try {(void)evaluate_scalar_property(c,200.0);} catch(...) {threw=true;}
        EXPECT_TRUE(threw);
    });
    run_case("density_and_enthalpy", [] {
        DensityControls d; d.rho0=1000.0; d.model=DensityModel::BOUSSINESQ; d.beta=1e-3;
        EXPECT_NEAR(evaluate_density(d,310.0),990.0,1e-12);
        d.model=DensityModel::BOUSSINESQ;
        bool boussinesq_threw=false; try {(void)evaluate_density(d,450.0);} catch(...) {boussinesq_threw=true;}
        EXPECT_TRUE(boussinesq_threw);
        d.model=DensityModel::IDEAL_GAS; d.gas_constant=287.05;
        EXPECT_NEAR(evaluate_density(d,300.0,101325.0),101325.0/(287.05*300.0),1e-12);
        ThermophysicalProperties p;
        p.heat_capacity.reference_value=1000.0;
        EXPECT_NEAR(p.enthalpy(400.0,300.0,64),100000.0,1e-8);
        p.heat_capacity.model=ScalarPropertyModel::LINEAR;
        p.heat_capacity.coefficients={2.0};
        EXPECT_NEAR(p.enthalpy(400.0,300.0,8),110000.0,1e-8);
        p.heat_capacity.model=ScalarPropertyModel::TABLE;
        p.heat_capacity.table.temperature={300,400,500};
        p.heat_capacity.table.value={1000,1200,1400};
        EXPECT_NEAR(p.enthalpy(500.0,300.0),120000.0,1e-8);
        p.heat_capacity.model=ScalarPropertyModel::LINEAR;
        p.heat_capacity.reference_value=1000.0;
        p.heat_capacity.coefficients={2.0};
        EXPECT_NEAR((p.cp(400.0)-p.cp(300.0))/100.0,2.0,1e-12);
    });
    run_case("thermal_turbulence_coupling", [] {
        EXPECT_NEAR(turbulent_thermal_diffusivity(0.01,1.0,0.9),0.0111111111111111,1e-12);
        EXPECT_NEAR(k_epsilon_eddy_viscosity(4.0,2.0),0.72,1e-12);
        EXPECT_NEAR(rng_kepsilon_eddy_viscosity(4.0,2.0),0.676,1e-12);
        EXPECT_NEAR(komega_eddy_viscosity(2.0,4.0),0.5,1e-12);
        EXPECT_NEAR(sst_eddy_viscosity(1.0,2.0,1.0),0.31,1e-12);
        EXPECT_NEAR(spalart_allmaras_nu_t(1e-4,0.01,1e-5),1e-4*std::pow(10.0,3)/(std::pow(10.0,3)+std::pow(7.1,3)),1e-12);
        EXPECT_NEAR(spalart_allmaras_nu_t(-1e-4,0.01,1e-5),0.0,1e-15);
        bool sa_threw=false; try {(void)spalart_allmaras_nu_t(-1e-4,0.01,1e-5,NegativeNuTildePolicy::REJECT);} catch(...) {sa_threw=true;}
        EXPECT_TRUE(sa_threw);
        EXPECT_NEAR(smagorinsky_nu_t(0.17,0.1,10.0),0.00289,1e-12);
        EXPECT_NEAR(wale_nu_t(0.5,0.1,10.0),0.025,1e-12);
        EXPECT_EQ(static_cast<int>(implementation_kind(AdvancedTurbulenceModel::SST)),
                  static_cast<int>(TurbulenceImplementationKind::CLOSURE));
        EXPECT_NEAR(des_length_scale(0.5,0.2,0.65),0.13,1e-12);
        EXPECT_NEAR(ddes_length_scale_from_rd(0.5,0.2,0.65,0.0),0.13,1e-12);
        EXPECT_NEAR(ddes_length_scale_from_rd(0.5,0.2,0.65,1.0),0.5,1e-12);
        EXPECT_NEAR(iddes_length_scale_from_rd(0.5,0.2,0.65,0.0,1.0),0.13,1e-12);
        EXPECT_NEAR(iddes_length_scale_from_rd(0.5,0.2,0.65,1.0,1.0),0.5,1e-12);
        EXPECT_NEAR(ddes_shielding(0.0),1.0,1e-12);
        EXPECT_LT(ddes_shielding(1.0),1e-10);
        EXPECT_NEAR(ddes_length_scale(0.5,0.2,0.65,1.0),0.13,1e-12);
        EXPECT_NEAR(ddes_length_scale(0.5,0.2,0.65,0.0),0.5,1e-12);
    });
    run_case("invalid_models_are_rejected", [] {
        bool threw=false; try { ScalarPropertyControls c; c.reference_value=-1; (void)evaluate_scalar_property(c,300); } catch(...) {threw=true;} EXPECT_TRUE(threw);
        threw=false; try { DensityControls d; d.model=DensityModel::BOUSSINESQ; d.beta=1; (void)evaluate_density(d,1000); } catch(...) {threw=true;} EXPECT_TRUE(threw);
        threw=false; try { TabulatedProperty t; t.temperature={300,300}; t.value={1,2}; (void)t.evaluate(300); } catch(...) {threw=true;} EXPECT_TRUE(threw);
    });
    return run_all();
}
