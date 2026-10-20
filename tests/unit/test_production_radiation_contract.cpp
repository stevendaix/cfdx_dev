#include "cfdx/physics/production_radiation_contract.h"
#include "common/test_harness.h"

#include <stdexcept>
#include <string>

using namespace cfdx::io;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {

// Declaring a model that the production driver does not orchestrate yet
// must throw the fail-closed "not coupled" rejection.
bool rejected_not_coupled(const char* model)
{
    CaseSetup setup;
    setup.radiation_model = model;
    try {
        validate_production_radiation_selection(setup);
    } catch (const std::invalid_argument& error) {
        return std::string(error.what()).find("not coupled") !=
               std::string::npos;
    }
    return false;
}

}  // namespace

int main()
{
    run_case("production_radiation_accepts_radiation_free_selection", [] {
        CaseSetup setup;
        // Unspecified declaration: legacy case files may omit the block.
        setup.radiation_model = "";
        validate_production_radiation_selection(setup);

        setup.radiation_model = "none";
        validate_production_radiation_selection(setup);

        setup.radiation_model = "off";
        validate_production_radiation_selection(setup);

        setup.radiation_model = "disabled";
        validate_production_radiation_selection(setup);
    });

    run_case("production_radiation_accepts_normalized_free_keys", [] {
        // Same normalization contract as the turbulence selection guard:
        // case and surrounding whitespace must not flip a radiation-free
        // declaration into a spurious reject.
        CaseSetup setup;
        setup.radiation_model = "NONE";
        validate_production_radiation_selection(setup);
        setup.radiation_model = " none ";
        validate_production_radiation_selection(setup);
        setup.radiation_model = "Off";
        validate_production_radiation_selection(setup);
        setup.radiation_model = "DISABLED";
        validate_production_radiation_selection(setup);
    });

    run_case("production_radiation_rejects_schema_declared_models", [] {
        // Every declaration the case schema documents is rejected while
        // the production driver does not consume its radiative
        // flux/source.
        EXPECT_TRUE(rejected_not_coupled("S2S"));
        EXPECT_TRUE(rejected_not_coupled("DO"));
        EXPECT_TRUE(rejected_not_coupled("solar"));
        // Lowercase spellings must be rejected too: normalization must
        // not turn a declared model into a radiation-free case.
        EXPECT_TRUE(rejected_not_coupled("s2s"));
        EXPECT_TRUE(rejected_not_coupled("do"));
    });

    run_case("production_radiation_rejects_available_physics_models", [] {
        // The radiation bricks exist (fiche 11: Rosseland, P1, S2S, DOM)
        // but their coupling to the production energy equation is not
        // integrated: declaring them must fail closed instead of
        // silently running radiation-free.
        EXPECT_TRUE(rejected_not_coupled("P1"));
        EXPECT_TRUE(rejected_not_coupled("p1"));
        EXPECT_TRUE(rejected_not_coupled("DOM"));
        EXPECT_TRUE(rejected_not_coupled("dom"));
        EXPECT_TRUE(rejected_not_coupled("Rosseland"));
        EXPECT_TRUE(rejected_not_coupled("rosseland"));
    });

    run_case("production_radiation_rejects_unknown_model", [] {
        EXPECT_TRUE(rejected_not_coupled("unrecognized-model"));
    });

    return run_all();
}
