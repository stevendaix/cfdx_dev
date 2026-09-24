#include "cfdx/physics/radiation.h"
#include "cfdx/physics/radiation_models.h"
#include "common/test_harness.h"
#include <cmath>
using namespace cfdx::physics; using namespace cfdx::testing;
int main() {
    run_case("M4_diffuse_gray_wall_reflection_uses_G_over_pi", [] {
        const double eps=0.4, T=500.0, G=1000.0;
        const double expected=eps*blackbody_intensity(T)+(1.0-eps)*G/M_PI;
        EXPECT_NEAR(gray_diffuse_wall_intensity(eps,T,G),expected,1e-12*std::max(1.0,std::abs(expected)));
    });
    run_case("M4_p1_radiative_equilibrium_source_is_zero", [] {
        const double T=600.0;
        const double G=4.0*M_PI*blackbody_intensity(T);
        EXPECT_NEAR(p1_radiative_source(0.7,G/(4.0*M_PI),T),0.0,1e-12);
    });
    run_case("M4_radiation_model_selector_respects_optical_regimes", [] {
        RadiationModelSelector s;
        EXPECT_TRUE(s.select(0.01,0.0,1.0)==RadiationApproximation::DOM);
        EXPECT_TRUE(s.select(0.2,0.0,1.0)==RadiationApproximation::P1);
        EXPECT_TRUE(s.select(4.0,0.0,1.0)==RadiationApproximation::Rosseland);
    });
    return run_all();
}
