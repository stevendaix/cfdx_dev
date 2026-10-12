// Issue #461 — N1 case -> selection -> report wiring.
//
// resolve_case_numerics consumes an explicit case numerics configuration and
// produces an auditable report: every entry is either resolved against the
// registry (with method id and verification status) or recorded as an error,
// with no silent fallback. The report is deterministic.

#include "cfdx/core/numerics/numerical_method_selection.h"
#include "cfdx/core/numerics/case_numerics_config.h"
#include "common/test_harness.h"

#include <algorithm>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("case_report_valid_config_is_deterministic", [] {
        // One valid entry per family present in the registry.
        std::vector<CaseNumericsEntry> config;
        for (const auto family : {
                 NumericalMethodFamily::Gradient,
                 NumericalMethodFamily::Convection,
                 NumericalMethodFamily::Diffusion,
                 NumericalMethodFamily::Interpolation,
                 NumericalMethodFamily::Temporal}) {
            const auto keys = configuration_keys_for_family(family);
            EXPECT_TRUE(!keys.empty());
            config.push_back({family, keys.front()});
        }

        const auto report = resolve_case_numerics("channel", config);
        EXPECT_TRUE(report.errors.empty());
        EXPECT_TRUE(report.valid());
        EXPECT_TRUE(report.resolved.size() == config.size());
        for (const auto& s : report.resolved) {
            EXPECT_TRUE(s.resolved);
            EXPECT_TRUE(!s.method_id.empty());
            EXPECT_TRUE(!s.method_name.empty());
        }
        const std::string a = format_numerics_report(report);
        const std::string b = format_numerics_report(report);
        EXPECT_TRUE(a == b);   // deterministic
        EXPECT_TRUE(a.find("case=channel") != std::string::npos);
        EXPECT_TRUE(a.find("scheme[") != std::string::npos);
    });

    run_case("case_report_resolves_explicit_schur_family", [] {
        const std::vector<CaseNumericsEntry> config = {
            {NumericalMethodFamily::PressureVelocity, "pressure_velocity.coupled"},
            {NumericalMethodFamily::Schur, "schur.pcd"},
        };
        const auto report = resolve_case_numerics("coupled", config);
        EXPECT_TRUE(report.errors.empty());
        EXPECT_TRUE(report.valid());
        EXPECT_TRUE(report.resolved.size() == 2);
        EXPECT_TRUE(report.resolved[1].family == NumericalMethodFamily::Schur);
        EXPECT_TRUE(report.resolved[1].method_id == "schur.pcd");
        EXPECT_TRUE(format_numerics_report(report).find("scheme[schur]=schur.pcd") != std::string::npos);
    });

    run_case("case_report_unknown_and_empty_keys_are_errors", [] {
        std::vector<CaseNumericsEntry> config = {
            {NumericalMethodFamily::Convection, "numerics.convection.does_not_exist"},
            {NumericalMethodFamily::Diffusion, ""},
        };
        const auto report = resolve_case_numerics("broken", config);
        EXPECT_TRUE(!report.valid());
        EXPECT_TRUE(report.errors.size() == 2);
        EXPECT_TRUE(report.resolved.empty());
        const std::string text = format_numerics_report(report);
        EXPECT_TRUE(text.find("error:") != std::string::npos);
    });

    run_case("case_report_wrong_family_is_an_error", [] {
        // A valid key of the wrong family must be rejected by the family check.
        std::vector<CaseNumericsEntry> config = {
            {NumericalMethodFamily::Diffusion, "numerics.convection.tvd.minmod"},
        };
        const auto report = resolve_case_numerics("mismatch", config);
        EXPECT_TRUE(!report.valid());
        EXPECT_TRUE(report.errors.size() == 1);
    });

    run_case("case_report_round_trips_every_registry_method", [] {
        // A config requesting every registered method (by its own key and
        // family) must resolve cleanly — a structural registry round trip at
        // the case level.
        std::vector<CaseNumericsEntry> config;
        for (const auto& method : numerical_method_registry())
            config.push_back({method.family, method.configuration_key});
        const auto report = resolve_case_numerics("all", config);
        EXPECT_TRUE(report.errors.empty());
        EXPECT_TRUE(report.resolved.size() == config.size());
    });

    return run_all();
}