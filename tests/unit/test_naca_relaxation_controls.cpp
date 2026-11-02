#include "../validation/naca_relaxation_controls.h"
#include "common/test_harness.h"
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

using namespace cfdx::testing;

static void setenv_checked(const char* key, const char* value) {
    if (::setenv(key, value, 1) != 0) throw std::runtime_error("setenv failed");
}
static void clear_env() {
    ::unsetenv("CFDX_NACA_AB_ALPHA_U");
    ::unsetenv("CFDX_NACA_AB_ALPHA_P");
}

int main() {
    clear_env();
    run_case("naca_ab_preserves_default", [] {
        cfdx::physics::IncompressibleSolverControls controls;
        controls.coupling.alpha_u = 0.6;
        controls.coupling.alpha_p = 0.25;
        cfdx::validation::apply_naca_relaxation_ab(controls);
        EXPECT_TRUE(std::abs(controls.coupling.alpha_u - 0.6) < 1e-15);
        EXPECT_TRUE(std::abs(controls.coupling.alpha_p - 0.25) < 1e-15);
    });
    run_case("naca_ab_experiment_is_opt_in", [] {
        setenv_checked("CFDX_NACA_AB_ALPHA_U", "0.35");
        setenv_checked("CFDX_NACA_AB_ALPHA_P", "0.05");
        cfdx::physics::IncompressibleSolverControls controls;
        cfdx::validation::apply_naca_relaxation_ab(controls);
        EXPECT_TRUE(std::abs(controls.coupling.alpha_u - 0.35) < 1e-15);
        EXPECT_TRUE(std::abs(controls.coupling.alpha_p - 0.05) < 1e-15);
        clear_env();
    });
    run_case("naca_ab_rejects_invalid_settings", [] {
        const char* invalid[] = {"0", "-0.5", "1.01", "nan", "inf", "0.35garbage", ""};
        for (const char* value : invalid) {
            setenv_checked("CFDX_NACA_AB_ALPHA_U", value);
            bool rejected = false;
            try { (void)cfdx::validation::read_naca_relaxation("CFDX_NACA_AB_ALPHA_U", 0.6); }
            catch (const std::invalid_argument&) { rejected = true; }
            EXPECT_TRUE(rejected);
        }
        clear_env();
    });
    return run_all();
}
