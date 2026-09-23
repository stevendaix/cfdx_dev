#include "cfdx/physics/turbulence_solver.h"
#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/spalart_allmaras.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>

using namespace cfdx::physics;
using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("M2_SA_canonical_closures", [] {
        SpalartAllmarasModel sa;
        EXPECT_NEAR(sa.fv1(1.0), 1.0/(1.0+7.1*7.1*7.1), 1e-14);
        EXPECT_TRUE(sa.fv2(0.0) >= 0.0);
        EXPECT_TRUE(sa.fv2(1.0) >= 0.0 && sa.fv2(1.0) <= 1.0);
        EXPECT_TRUE(sa.ft2(1.0) > 0.0 && sa.ft2(1.0) < 1.2);
        EXPECT_NEAR(sa.cw1, sa.cb1/(sa.kappa*sa.kappa)+(1.0+sa.cb2)/sa.sigma, 1e-14);
        EXPECT_TRUE(sa.turbulent_viscosity(1.0,1e-4,1e-5) > 0.0);
    });

    run_case("M2_wall_treatment_identities", [] {
        const double k=0.25, y=0.01;
        EXPECT_NEAR(wall_epsilon_from_k(k,y),
                    std::pow(0.09,0.75)*std::pow(k,1.5)/y,1e-14);
        EXPECT_NEAR(wall_omega_from_k(k,y),
                    std::sqrt(k)/(std::sqrt(0.09)*y),1e-14);
        EXPECT_TRUE(wall_epsilon_from_k(k,y) > 0.0);
        EXPECT_TRUE(wall_omega_from_k(k,y) > 0.0);
    });

    run_case("M2_SST_blending_bounds_and_limits", [] {
        const auto near_wall=compute_sst_blending(1.0,10.0,1e-4,1e-5,0.09);
        const auto far_wall=compute_sst_blending(1.0,1.0,10.0,1e-5,0.09);
        EXPECT_TRUE(near_wall.first>=0.0 && near_wall.first<=1.0);
        EXPECT_TRUE(near_wall.second>=0.0 && near_wall.second<=1.0);
        EXPECT_TRUE(far_wall.first>=0.0 && far_wall.first<=1.0);
        EXPECT_TRUE(far_wall.second>=0.0 && far_wall.second<=1.0);
        EXPECT_TRUE(near_wall.first >= far_wall.first);
    });

    run_case("M2_controls_reject_nonfinite_values", [] {
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::KEPSILON;
        c.k_min=std::numeric_limits<double>::quiet_NaN();
        bool rejected=false;
        try { validate_turbulence_controls(c); } catch(const std::invalid_argument&) { rejected=true; }
        EXPECT_TRUE(rejected);
    });

    run_case("M2_wall_invalid_inputs_rejected", [] {
        bool rejected=false;
        try { (void)wall_epsilon_from_k(1.0,0.0); }
        catch(const std::invalid_argument&) { rejected=true; }
        EXPECT_TRUE(rejected);
    });

    return run_all();
}
