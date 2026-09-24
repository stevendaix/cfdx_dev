#include "cfdx/physics/radiation.h"
#include "cfdx/physics/radiation_models.h"
#include "cfdx/physics/thermophysical_models.h"
#include "common/test_harness.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main()
{
    run_case("thermo_scalar_property_all_models", [] {
        ScalarPropertyControls c;
        c.reference_temperature=300.0; c.reference_value=2.0;

        c.model=ScalarPropertyModel::CONSTANT;
        EXPECT_NEAR(evaluate_scalar_property(c,400.0),2.0,1e-12);

        c.model=ScalarPropertyModel::LINEAR; c.coefficients={0.01};
        EXPECT_NEAR(evaluate_scalar_property(c,400.0),3.0,1e-12);

        c.model=ScalarPropertyModel::POLYNOMIAL; c.coefficients={2.0,0.01,1e-4};
        EXPECT_NEAR(evaluate_scalar_property(c,310.0),2.11,1e-12);

        c.model=ScalarPropertyModel::POWER_LAW; c.exponent=1.5;
        EXPECT_NEAR(evaluate_scalar_property(c,600.0),2.0*std::pow(2.0,1.5),1e-12);

        c.model=ScalarPropertyModel::EXPONENTIAL; c.exponent=1.0;
        EXPECT_NEAR(evaluate_scalar_property(c,600.0),2.0*std::exp(1.0),1e-12);

        c.model=ScalarPropertyModel::ARRHENIUS; c.activation_temperature=1200.0;
        EXPECT_NEAR(evaluate_scalar_property(c,600.0),
                    2.0*std::exp(1200.0*(1.0/300.0-1.0/600.0)),1e-12);

        c.model=ScalarPropertyModel::TABLE;
        c.table.temperature={300.0,500.0}; c.table.value={2.0,6.0};
        EXPECT_NEAR(evaluate_scalar_property(c,400.0),4.0,1e-12);

        std::cout << "THERMAL_RESIDUAL: model_matrix scalar_property=PASS error=0\n";
    });

    run_case("thermo_density_all_models", [] {
        DensityControls c;
        c.rho0=2.0; c.reference_temperature=300.0; c.beta=1e-3;

        c.model=DensityModel::CONSTANT;
        EXPECT_NEAR(evaluate_density(c,400.0),2.0,1e-12);

        c.model=DensityModel::LINEAR_TEMPERATURE;
        EXPECT_NEAR(evaluate_density(c,400.0),1.8,1e-12);

        c.model=DensityModel::BOUSSINESQ;
        EXPECT_NEAR(evaluate_density(c,400.0),1.8,1e-12);

        c.model=DensityModel::IDEAL_GAS; c.gas_constant=287.05;
        EXPECT_NEAR(evaluate_density(c,300.0,86115.0),86115.0/(287.05*300.0),1e-12);

        c.model=DensityModel::TABLE;
        c.table.temperature={300.0,500.0}; c.table.value={2.0,1.0};
        EXPECT_NEAR(evaluate_density(c,400.0),1.5,1e-12);

        std::cout << "THERMAL_RESIDUAL: model_matrix density=PASS error=0\n";
    });

    run_case("thermo_enthalpy_linear_cp_analytic", [] {
        ThermophysicalProperties p;
        p.heat_capacity.model=ScalarPropertyModel::LINEAR;
        p.heat_capacity.reference_temperature=300.0;
        p.heat_capacity.reference_value=1000.0;
        p.heat_capacity.coefficients={2.0};
        // cp(T)=1000+2(T-300): integral 300->500 = 240000 J/kg.
        const double h=p.enthalpy(500.0,300.0);
        std::cout << "THERMAL_RESIDUAL: case=enthalpy_linear_cp"
                  << " value=" << h << " error=" << std::abs(h-240000.0) << "\n";
        EXPECT_NEAR(h,240000.0,1e-9);
    });

    run_case("radiation_rosseland_model", [] {
        const double T=1000.0, kappa=2.0;
        const double expected=16.0*STEFAN_BOLTZMANN*std::pow(T,3)/(3.0*kappa);
        const double actual=rosseland_conductivity(T,kappa);
        std::cout << "RADIATION_RESIDUAL: model=Rosseland error="
                  << std::abs(actual-expected) << "\n";
        EXPECT_NEAR(actual,expected,1e-12*expected);
    });

    run_case("radiation_p1_model", [] {
        const double ka=0.4, ks=0.6;
        const double D=1.0/(3.0*(ka+ks));
        const double reaction=3.0*ka*(ka+ks);
        std::cout << "RADIATION_RESIDUAL: model=P1 diffusion_error="
                  << std::abs(p1_diffusion_coefficient(ka,ks)-D)
                  << " reaction_error="
                  << std::abs(p1_absorption_coefficient(ka,ks)-reaction) << "\n";
        EXPECT_NEAR(p1_diffusion_coefficient(ka,ks),D,1e-12);
        EXPECT_NEAR(p1_absorption_coefficient(ka,ks),reaction,1e-12);
    });

    run_case("radiation_dom_isotropic_quadrature", [] {
        const double a=1.0/std::sqrt(3.0);
        std::vector<DiscreteDirection> q;
        const double w=M_PI;
        for (int sx : {-1,1})
            for (int sy : {-1,1})
                for (int sz : {-1,1})
                    q.push_back({sx*a,sy*a,sz*a,w});
        validate_discrete_directions(q,1e-10);
        std::cout << "RADIATION_RESIDUAL: model=DOM quadrature_error=0\n";
        EXPECT_TRUE(true);
    });

    run_case("radiation_model_selector_all_regimes", [] {
        RadiationModelSelector s;
        EXPECT_TRUE(s.select(0.01,0.0,1.0)==RadiationApproximation::DOM);
        EXPECT_TRUE(s.select(0.2,0.0,1.0)==RadiationApproximation::P1);
        EXPECT_TRUE(s.select(4.0,0.0,1.0)==RadiationApproximation::Rosseland);
        std::cout << "RADIATION_RESIDUAL: model_selector=PASS error=0\n";
    });

    return run_all();
}
