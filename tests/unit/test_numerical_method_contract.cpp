#include "cfdx/core/numerics/numerical_method_contract.h"

#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::core;

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    NumericalMethodContract tvd{
        "convection.tvd.minmod", "TVD MinMod",
        NumericalMethodFamily::Convection, VerificationStatus::Implemented,
        ConservationContract::LocalFaceConservative, true, true, false, 1, 0,
        "MUSCL face reconstruction with MinMod limiter",
        "numerics.convection.tvd.minmod",
        {"boundedness", "conservation", "smooth-field refinement", "steep-gradient regression"}
    };
    validate_numerical_method_contract(tvd);
    require(std::string(to_string(tvd.family)) == "convection", "family metadata mismatch");
    require(std::string(to_string(tvd.status)) == "implemented", "status metadata mismatch");
    require(std::string(to_string(tvd.conservation)) == "local_face_conservative",
            "conservation metadata mismatch");

    NumericalMethodContract bad = tvd;
    bad.id.clear();
    bool rejected = false;
    try { validate_numerical_method_contract(bad); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "invalid contract was not rejected");

    bad = tvd;
    bad.status = VerificationStatus::Validated;
    bad.verification_requirements.clear();
    rejected = false;
    try { validate_numerical_method_contract(bad); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "validated contract without evidence requirements was accepted");

    std::cout << "numerical method contract checks: PASS\n";
    return 0;
}
