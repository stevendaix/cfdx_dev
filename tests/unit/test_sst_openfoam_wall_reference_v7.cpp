#include "cfdx/physics/sst_openfoam_wall_reference_v7.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>
using namespace cfdx::physics;
using namespace cfdx::testing;
int main() {
    run_case("openfoam_omega_viscous_limit_is_six_nu_beta_y2", [] {
        const auto w=openfoam_omega_wall_reference(0.0,1e-5,1e-3,0.0);
        EXPECT_NEAR(w.reference_omega,6e-5/(0.075e-6),1e-10);
        EXPECT_NEAR(w.y_plus_from_k,0.0,0.0);
    });
    run_case("openfoam_omega_log_limit_uses_k_based_friction_velocity", [] {
        const double k=100.0,nu=1e-5,y=1e-3;
        const auto w=openfoam_omega_wall_reference(k,nu,y,20.0);
        const double expected=std::sqrt(std::sqrt(0.09)*k)/
            (std::sqrt(0.09)*0.41*y);
        EXPECT_NEAR(w.reference_omega,expected,1e-12*expected);
        EXPECT_TRUE(w.y_plus_from_k>11.06);
    });
    run_case("openfoam_omega_blended_wall_has_both_branches", [] {
        const auto w=openfoam_omega_wall_reference(0.0004,1e-5,0.002,30.0,
            ReferenceOmegaWallBlend::LaminarWeighted);
        EXPECT_TRUE(w.laminar_fraction>0.0 && w.laminar_fraction<1.0);
        const double mixed=w.laminar_fraction*w.viscous_omega+
            (1.0-w.laminar_fraction)*w.logarithmic_omega;
        EXPECT_NEAR(w.reference_omega,mixed,1e-12);
    });
    run_case("openfoam_omega_ref_rejects_nonphysical_inputs", [] {
        EXPECT_THROW(openfoam_omega_wall_reference(-1.0,1e-5,0.001,1.0),std::invalid_argument);
        EXPECT_THROW(openfoam_omega_wall_reference(1.0,0.0,0.001,1.0),std::invalid_argument);
        EXPECT_THROW(openfoam_omega_wall_reference(1.0,1e-5,0.0,1.0),std::invalid_argument);
        EXPECT_THROW(openfoam_omega_wall_reference(1.0,1e-5,0.001,
            std::numeric_limits<double>::quiet_NaN()),std::invalid_argument);
    });
    return run_all();
}
