#include "cfdx/physics/radiation.h"
#include "cfdx/physics/radiation_models.h"
#include "cfdx/physics/radiation_s2s.h"
#include "cfdx/physics/thermophysical_models.h"
#include "common/test_harness.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main()
{
    run_case("thermal_scalar_models_reference_points", [] {
        ScalarPropertyControls c;
        c.reference_temperature=300.0; c.reference_value=2.0;
        c.model=ScalarPropertyModel::CONSTANT;
        EXPECT_NEAR(evaluate_scalar_property(c,450),2.0,1e-12);
        c.model=ScalarPropertyModel::LINEAR; c.coefficients={0.01};
        EXPECT_NEAR(evaluate_scalar_property(c,400),3.0,1e-12);
        c.model=ScalarPropertyModel::POLYNOMIAL; c.coefficients={2.0,0.01,1e-4};
        EXPECT_NEAR(evaluate_scalar_property(c,310),2.11,1e-12);
        c.model=ScalarPropertyModel::POWER_LAW; c.exponent=1.5;
        EXPECT_NEAR(evaluate_scalar_property(c,600),2.0*std::pow(2.0,1.5),1e-12);
        c.model=ScalarPropertyModel::EXPONENTIAL; c.exponent=1.0;
        EXPECT_NEAR(evaluate_scalar_property(c,600),2.0*std::exp(1.0),1e-12);
        c.model=ScalarPropertyModel::ARRHENIUS; c.activation_temperature=1200.0;
        EXPECT_NEAR(evaluate_scalar_property(c,600),
                    2.0*std::exp(1200.0*(1.0/300.0-1.0/600.0)),1e-12);
        std::cout << "THERMAL_REGRESSION: scalar_models=PASS\n";
    });

    run_case("thermal_table_and_bounds_regression", [] {
        ScalarPropertyControls c;
        c.reference_temperature=300; c.reference_value=1;
        c.model=ScalarPropertyModel::TABLE;
        c.table.temperature={300,400,500}; c.table.value={10,20,30};
        EXPECT_NEAR(evaluate_scalar_property(c,350),15,1e-12);
        EXPECT_NEAR(evaluate_scalar_property(c,200),10,1e-12);
        c.table.extrapolation=ExtrapolationPolicy::LINEAR;
        EXPECT_NEAR(evaluate_scalar_property(c,250),5,1e-12);
        c.table.extrapolation=ExtrapolationPolicy::REJECT;
        EXPECT_THROW(evaluate_scalar_property(c,250), std::out_of_range);
        c.model=ScalarPropertyModel::LINEAR; c.coefficients={0.1};
        c.enforce_bounds=true; c.minimum=1; c.maximum=3;
        EXPECT_NEAR(evaluate_scalar_property(c,400),3,1e-12);
        std::cout << "THERMAL_REGRESSION: table_bounds=PASS\n";
    });

    run_case("thermal_density_models_regression", [] {
        DensityControls c; c.rho0=2; c.reference_temperature=300; c.beta=1e-3;
        c.model=DensityModel::CONSTANT;
        EXPECT_NEAR(evaluate_density(c,400),2,1e-12);
        c.model=DensityModel::LINEAR_TEMPERATURE;
        EXPECT_NEAR(evaluate_density(c,400),1.8,1e-12);
        c.model=DensityModel::BOUSSINESQ;
        EXPECT_NEAR(evaluate_density(c,400),1.8,1e-12);
        EXPECT_THROW(evaluate_density(c,450), std::domain_error);
        c.model=DensityModel::IDEAL_GAS; c.gas_constant=287.05;
        EXPECT_NEAR(evaluate_density(c,300,86115),86115/(287.05*300),1e-12);
        c.model=DensityModel::TABLE;
        c.table.temperature={300,500}; c.table.value={2,1};
        EXPECT_NEAR(evaluate_density(c,400),1.5,1e-12);
        std::cout << "THERMAL_REGRESSION: density_models=PASS\n";
    });

    run_case("thermal_enthalpy_models_regression", [] {
        ThermophysicalProperties p;
        p.heat_capacity.model=ScalarPropertyModel::CONSTANT;
        p.heat_capacity.reference_value=1000;
        EXPECT_NEAR(p.enthalpy(400,300),100000,1e-8);
        p.heat_capacity.model=ScalarPropertyModel::LINEAR;
        p.heat_capacity.reference_temperature=300;
        p.heat_capacity.reference_value=1000;
        p.heat_capacity.coefficients={2};
        EXPECT_NEAR(p.enthalpy(500,300),240000,1e-8);
        p.heat_capacity.model=ScalarPropertyModel::POLYNOMIAL;
        p.heat_capacity.coefficients={1000,2};
        EXPECT_NEAR(p.enthalpy(500,300),240000,1e-8);
        p.heat_capacity.model=ScalarPropertyModel::TABLE;
        p.heat_capacity.table.temperature={300,400,500};
        p.heat_capacity.table.value={1000,1200,1400};
        EXPECT_NEAR(p.enthalpy(500,300),240000,1e-8);
        EXPECT_NEAR(p.enthalpy(300,500),-240000,1e-8);
        std::cout << "THERMAL_REGRESSION: enthalpy_models=PASS\n";
    });

    run_case("thermal_invalid_input_guards", [] {
        ScalarPropertyControls c;
        EXPECT_THROW(evaluate_scalar_property(c,0), std::invalid_argument);
        c.enforce_bounds=true; c.minimum=2; c.maximum=1;
        EXPECT_THROW(evaluate_scalar_property(c,300), std::invalid_argument);
        DensityControls d;
        d.model=DensityModel::IDEAL_GAS; d.gas_constant=0;
        EXPECT_THROW(evaluate_density(d,300,101325), std::invalid_argument);
        std::cout << "THERMAL_REGRESSION: invalid_guards=PASS\n";
    });

    run_case("radiation_blackbody_and_gray_surface", [] {
        const double T=1000, eps=0.4, G=2000;
        const double Eb=STEFAN_BOLTZMANN*std::pow(T,4);
        EXPECT_NEAR(blackbody_emissive_power(T),Eb,1e-12*Eb);
        EXPECT_NEAR(gray_surface_emissivity_flux(eps,T,G),eps*(Eb-G),1e-12*Eb);
        std::cout << "RADIATION_REGRESSION: blackbody_gray=PASS\n";
    });

    run_case("radiation_p1_source_equilibrium_and_sign", [] {
        const double T=700, ka=0.7;
        const double Ib=blackbody_intensity(T);
        EXPECT_NEAR(p1_radiative_source(ka,Ib,T),0,1e-12);
        EXPECT_TRUE(p1_radiative_source(ka,0,T)>0);
        EXPECT_TRUE(p1_radiative_source(ka,2*Ib,T)<0);
        std::cout << "RADIATION_REGRESSION: p1_source=PASS\n";
    });

    run_case("radiation_rosseland_and_p1_coefficients", [] {
        const double T=1000, ka=0.4, ks=0.6;
        const double k=16*STEFAN_BOLTZMANN*std::pow(T,3)/(3*(ka+ks));
        EXPECT_NEAR(rosseland_conductivity(T,ka,ks),k,1e-12*k);
        EXPECT_NEAR(p1_diffusion_coefficient(ka,ks),1.0/3.0,1e-12);
        EXPECT_NEAR(p1_absorption_coefficient(ka,ks),0.72,1e-12);
        std::cout << "RADIATION_REGRESSION: rosseland_p1=PASS\n";
    });

    run_case("radiation_selector_thresholds_and_surface_modes", [] {
        RadiationModelSelector s;
        EXPECT_TRUE(s.select(0.1,0,1)==RadiationApproximation::P1);
        EXPECT_TRUE(s.select(3.0,0,1)==RadiationApproximation::Rosseland);
        EXPECT_TRUE(s.select(0.099999,0,1)==RadiationApproximation::DOM);
        EXPECT_TRUE(s.select_surface_to_surface(false)==RadiationApproximation::S2S);
        EXPECT_TRUE(s.select_surface_to_surface(true)==RadiationApproximation::S2S_MONTE_CARLO);
        std::cout << "RADIATION_REGRESSION: selector=PASS\n";
    });

    run_case("radiation_dom_moments", [] {
        const double a=1/std::sqrt(3.0);
        std::vector<DiscreteDirection> q;
        for(int sx:{-1,1}) for(int sy:{-1,1}) for(int sz:{-1,1})
            q.push_back({sx*a,sy*a,sz*a,M_PI/2.0});
        validate_discrete_directions(q,1e-10);
        std::cout << "RADIATION_REGRESSION: dom_moments=PASS\n";
    });

    run_case("radiation_view_factor_reciprocity_and_closure", [] {
        const std::vector<double> A{2,1};
        const std::vector<double> F{0,0.5,1,0};
        validate_s2s_view_factors(F,A,1e-12);
        EXPECT_NEAR(A[0]*F[1],A[1]*F[2],1e-12);
        std::cout << "RADIATION_REGRESSION: view_factors=PASS\n";
    });

    run_case("s2s_radiosity_energy_balance", [] {
        const std::vector<double> A{1,1}, eps{0.7,0.5}, T{800,400};
        const std::vector<double> F{0,1,1,0};
        const auto r=solve_s2s_radiosity(A,eps,T,F);
        const double oracle=STEFAN_BOLTZMANN*(std::pow(800.0,4)-std::pow(400.0,4))/
            ((1-eps[0])/eps[0]+1+(1-eps[1])/eps[1]);
        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(r.net_flux[0],oracle,1e-10*std::abs(oracle));
        EXPECT_NEAR(r.net_flux[0],-r.net_flux[1],1e-10*std::abs(oracle));
        EXPECT_NEAR(r.energy_balance_error,0,1e-12);
        std::cout << "RADIATION_REGRESSION: s2s_energy=PASS error=" << r.energy_balance_error << "\n";
    });

    run_case("s2s_external_irradiation_limit", [] {
        const std::vector<double> A{1,1}, eps{1,1}, T{500,500};
        const std::vector<double> F{0,1,1,0};
        const auto r=solve_s2s_radiosity(A,eps,T,F,{},std::vector<double>{1000,1000});
        EXPECT_TRUE(r.converged);
        EXPECT_NEAR(r.total_power,0,1e-10);
        std::cout << "RADIATION_REGRESSION: s2s_external=PASS\n";
    });

    run_case("radiation_invalid_input_guards", [] {
        EXPECT_THROW(blackbody_emissive_power(-1), std::invalid_argument);
        EXPECT_THROW(rosseland_conductivity(1000,0), std::invalid_argument);
        EXPECT_THROW(p1_diffusion_coefficient(0,0), std::invalid_argument);
        EXPECT_THROW(validate_view_factor_matrix({0,1,0,0},2), std::invalid_argument);
        std::cout << "RADIATION_REGRESSION: invalid_guards=PASS\n";
    });

    std::cout << "THERMAL_RADIATION_REGRESSION: coverage=thermal_properties,density,enthalpy,"
                 "thermal_guards,blackbody,gray,P1,Rosseland,DOM,selector,S2S,guards PASS\n";
    return run_all();
}
