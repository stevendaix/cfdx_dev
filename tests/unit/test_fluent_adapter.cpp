// M0.10-T03: Tests for Fluent .cas/.dat adapter
#include "cfdx/io/fluent/fluent_adapter.h"
#include "cfdx/io/cfdx_io/io_interface.h"
#include "cfdx/io/cfdx_io/case_schema.h"
#include "common/test_harness.h"

#include <string>

using namespace cfdx::io::fluent;
using namespace cfdx::testing;

int main() {
    run_case("adapter_metadata", []() {
        FluentAdapter adapter;
        EXPECT_TRUE(std::string(adapter.solver_name()) == "Fluent");
        EXPECT_TRUE(std::string(adapter.format_name()) == "cas_dat");
    });

    run_case("detect_source", []() {
        FluentAdapter adapter;
        cfdx::io::SourceInfo info;
        bool detected = adapter.detect_source("tests/data/fluent/cavity.cas", info);
        EXPECT_TRUE(detected);
        EXPECT_TRUE(info.solver == "Fluent");
        EXPECT_TRUE(info.format == "cas_dat");
    });

    run_case("detect_source_not_found", []() {
        FluentAdapter adapter;
        cfdx::io::SourceInfo info;
        bool detected = adapter.detect_source("tests/data/fluent/nonexistent.cas", info);
        EXPECT_FALSE(detected);
    });

    run_case("parse_cas", []() {
        FluentAdapter adapter;
        bool ok = adapter.parse_cas("tests/data/fluent/cavity.cas");
        EXPECT_TRUE(ok);
        EXPECT_TRUE(adapter.zones().size() >= 8);
        EXPECT_TRUE(adapter.setup().materials.size() == 3);

        const auto& mat = adapter.setup().materials;
        bool has_air = false;
        for (const auto& m : mat) {
            if (m.name == "air") has_air = true;
        }
        EXPECT_TRUE(has_air);
    });

    run_case("convert_full", []() {
        FluentAdapter adapter;
        cfdx::io::ConversionResult result;
        bool ok = adapter.convert("tests/data/fluent/cavity.cas", result);
        EXPECT_TRUE(ok);
        EXPECT_FALSE(result.gap_report.has_blocking());
        EXPECT_TRUE(result.setup.boundary_conditions.size() >= 8);
        EXPECT_TRUE(result.setup.materials.size() == 3);
    });

    run_case("import_results_unsupported", []() {
        FluentAdapter adapter;
        cfdx::io::ConversionResult result;
        bool ok = adapter.import_results("tests/data/fluent/cavity.dat", result);
        EXPECT_FALSE(ok);
        EXPECT_TRUE(result.gap_report.has_blocking());
    });

    return run_all();
}
