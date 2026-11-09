#include "../validation/naca_coupled_selection.h"
#include "common/test_harness.h"
#include <stdexcept>
#include <string>

using namespace cfdx::testing;
using namespace cfdx::physics;
using namespace cfdx::validation;

int main() {
    run_case("default_stays_simple_and_existing_pcd_default", [] {
        const auto s = parse_naca_coupling_selection(nullptr, nullptr);
        EXPECT_TRUE(s.algorithm == PressureVelocityAlgorithm::SIMPLE);
        EXPECT_TRUE(s.schur_model == CoupledSchurModel::PCD);
        EXPECT_TRUE(std::string(naca_algorithm_name(s.algorithm)) == "SIMPLE");
    });
    run_case("simplec_preserved", [] {
        const auto s = parse_naca_coupling_selection("SIMPLEC", nullptr);
        EXPECT_TRUE(s.algorithm == PressureVelocityAlgorithm::SIMPLEC);
        EXPECT_TRUE(std::string(naca_algorithm_name(s.algorithm)) == "SIMPLEC");
    });
    run_case("all_six_schur_models_explicit_and_reproducible", [] {
        for (const char* name : {"block_local", "pcd", "lsc", "bfbt", "simple", "simplec"}) {
            const auto s = parse_naca_coupling_selection("COUPLED", name);
            EXPECT_TRUE(s.algorithm == PressureVelocityAlgorithm::COUPLED);
            EXPECT_TRUE(std::string(to_string(s.schur_model)) == name);
            EXPECT_TRUE(std::string(naca_algorithm_name(s.algorithm)) == "COUPLED");
        }
    });
    run_case("unsupported_or_implicit_selection_is_rejected", [] {
        EXPECT_THROW(parse_naca_coupling_selection("COUPLED", nullptr), std::invalid_argument);
        EXPECT_THROW(parse_naca_coupling_selection("COUPLED", "invalid"), std::invalid_argument);
        EXPECT_THROW(parse_naca_coupling_selection("SIMPLE", "pcd"), std::invalid_argument);
        EXPECT_THROW(parse_naca_coupling_selection("SIMPLEC", "lsc"), std::invalid_argument);
        EXPECT_THROW(parse_naca_coupling_selection("PISO", nullptr), std::invalid_argument);
        EXPECT_THROW(parse_naca_coupling_selection("", nullptr), std::invalid_argument);
    });
    return run_all();
}