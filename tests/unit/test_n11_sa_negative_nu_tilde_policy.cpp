// N11 F1 (#480): route the existing NegativeNuTildePolicy knob through the
// Spalart-Allmaras transport path. CLAMP_ZERO is the historical default and
// must stay unchanged; REJECT is the fail-closed qualification mode.
#include "cfdx/physics/turbulence_transport.h"

#include "common/test_harness.h"

#include <cmath>
#include <limits>
#include <stdexcept>

using namespace cfdx::testing;

namespace {

cfdx::physics::TurbulenceTransportControls sa_controls(
    cfdx::physics::NegativeNuTildePolicy policy)
{
    cfdx::physics::TurbulenceTransportControls c;
    c.model = cfdx::physics::TurbulenceModel::SPALART_ALLMARAS;
    c.molecular_viscosity = 1.0e-5;
    c.negative_nu_tilde_policy = policy;
    return c;
}

// Independent transcription of fiche 09, section 2 (SA):
// nu_t = nu_tilde * fv1, fv1 = chi^3/(chi^3 + cv1^3), chi = nu_tilde/nu.
double sa_reference_formula(double nu_tilde, double molecular_nu)
{
    const double chi = nu_tilde / molecular_nu;
    const double chi3 = chi * chi * chi;
    const double cv1 = 7.1;
    return nu_tilde * chi3 / (chi3 + cv1 * cv1 * cv1);
}

void default_policy_is_clamp_zero()
{
    cfdx::physics::TurbulenceTransportControls c;
    EXPECT_TRUE(c.negative_nu_tilde_policy ==
                cfdx::physics::NegativeNuTildePolicy::CLAMP_ZERO);
}

void clamp_zero_clamps_negative_nu_tilde()
{
    const auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::CLAMP_ZERO);
    const double nu = cfdx::physics::turbulence_nu_t(0.0, -1.0e-4, 0.0, 1.0, c);
    EXPECT_TRUE(std::isfinite(nu));
    EXPECT_NEAR(nu, 0.0, 1e-15);
}

void clamp_zero_keeps_positive_formula()
{
    const auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::CLAMP_ZERO);
    const double nu = cfdx::physics::turbulence_nu_t(0.0, 2.0e-5, 0.0, 1.0, c);
    EXPECT_NEAR(nu, sa_reference_formula(2.0e-5, 1.0e-5), 1e-18);
}

void reject_throws_domain_error_on_negative()
{
    const auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::REJECT);
    EXPECT_THROW_WITH(
        cfdx::physics::turbulence_nu_t(0.0, -1.0e-4, 0.0, 1.0, c),
        std::domain_error, "Spalart-Allmaras negative nu_tilde");
}

void reject_throws_invalid_argument_on_nonfinite()
{
    const auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::REJECT);
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW_WITH(
        cfdx::physics::turbulence_nu_t(0.0, nan_value, 0.0, 1.0, c),
        std::invalid_argument, "Spalart-Allmaras invalid inputs");
}

void reject_keeps_positive_path_unchanged()
{
    const auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::REJECT);
    const double nu = cfdx::physics::turbulence_nu_t(0.0, 2.0e-5, 0.0, 1.0, c);
    EXPECT_NEAR(nu, sa_reference_formula(2.0e-5, 1.0e-5), 1e-18);
}

void molecular_viscosity_guard_intact()
{
    auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::CLAMP_ZERO);
    c.molecular_viscosity = 0.0;
    EXPECT_THROW_WITH(
        cfdx::physics::turbulence_nu_t(0.0, 2.0e-5, 0.0, 1.0, c),
        std::invalid_argument, "SA eddy viscosity requires molecular viscosity");
}

void consistent_with_audited_closure_wrapper()
{
    const auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::CLAMP_ZERO);
    const double transport = cfdx::physics::turbulence_nu_t(0.0, 2.0e-5, 0.0, 1.0, c);
    const double wrapper = cfdx::physics::spalart_allmaras_nu_t(
        2.0e-5, 1.0, 1.0e-5, cfdx::physics::NegativeNuTildePolicy::CLAMP_ZERO);
    EXPECT_NEAR(transport, wrapper, 1e-15);
}

} // namespace

int main()
{
    run_case("SA transport: default policy is CLAMP_ZERO (historical)", default_policy_is_clamp_zero);
    run_case("SA transport: CLAMP_ZERO clamps negative nu_tilde to zero", clamp_zero_clamps_negative_nu_tilde);
    run_case("SA transport: CLAMP_ZERO keeps the positive formula unchanged", clamp_zero_keeps_positive_formula);
    run_case("SA transport: REJECT throws domain_error on negative nu_tilde", reject_throws_domain_error_on_negative);
    run_case("SA transport: REJECT throws invalid_argument on non-finite nu_tilde", reject_throws_invalid_argument_on_nonfinite);
    run_case("SA transport: REJECT keeps the positive path unchanged", reject_keeps_positive_path_unchanged);
    run_case("SA transport: molecular viscosity guard intact", molecular_viscosity_guard_intact);
    run_case("SA transport: consistent with the audited closure wrapper", consistent_with_audited_closure_wrapper);
    return run_all();
}
