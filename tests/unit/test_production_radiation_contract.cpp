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

// Declaring S2S without the conjugate energy model must be rejected with
// an explicit "conjugate" pointer at the energy stage that carries the
// radiative source/flux.
bool rejected_without_conjugate(const char* model, const char* energy_model)
{
    CaseSetup setup;
    setup.radiation_model = model;
    setup.energy_model = energy_model;
    try {
        validate_production_radiation_selection(setup);
    } catch (const std::invalid_argument& error) {
        return std::string(error.what()).find("conjugate") !=
               std::string::npos;
    }
    return false;
}

// S2S with the conjugate energy model is accepted (issue #885 stage 2).
bool accepted_s2s(const char* model)
{
    CaseSetup setup;
    setup.radiation_model = model;
    setup.energy_model = "conjugate";
    try {
        validate_production_radiation_selection(setup);
        return true;
    } catch (...) {
        return false;
    }
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

    run_case("production_radiation_s2s_requires_conjugate_energy_model", [] {
        // S2S is the first coupled radiation declaration (issue #885
        // stage 2): accepted only when the conjugate energy model carries
        // it, rejected with an explicit "conjugate" pointer otherwise.
        EXPECT_TRUE(accepted_s2s("S2S"));
        EXPECT_TRUE(accepted_s2s("s2s"));
        EXPECT_TRUE(accepted_s2s(" s2s "));
        EXPECT_TRUE(rejected_without_conjugate("S2S", ""));
        EXPECT_TRUE(rejected_without_conjugate("s2s", "none"));
        EXPECT_TRUE(rejected_without_conjugate("S2S", "isothermal"));
    });

    run_case("production_radiation_rejects_schema_declared_models", [] {
        // Remaining declared models are rejected while the production
        // driver does not consume their radiative flux/source.
        EXPECT_TRUE(rejected_not_coupled("DO"));
        EXPECT_TRUE(rejected_not_coupled("solar"));
        // Lowercase spellings must be rejected too: normalization must
        // not turn a declared model into a radiation-free case.
        EXPECT_TRUE(rejected_not_coupled("do"));
    });

    run_case("production_radiation_rejects_available_physics_models", [] {
        // The radiation bricks exist (fiche 11: Rosseland, P1, DOM) but
        // their coupling to the production energy equation is not
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
