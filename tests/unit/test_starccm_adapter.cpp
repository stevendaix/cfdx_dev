// M0.10-T03: Tests for STAR-CCM+ .sim adapter
#include "cfdx/io/starccm/starccm_adapter.h"
#include "cfdx/io/cfdx_io/io_interface.h"
#include "common/test_harness.h"

#include <string>

using namespace cfdx::io::starccm;
using namespace cfdx::testing;

int main() {
    run_case("adapter_metadata", []() {
        StarCCMAdapter adapter;
        EXPECT_TRUE(std::string(adapter.solver_name()) == "STAR-CCM+");
        EXPECT_TRUE(std::string(adapter.format_name()) == "sim");
    });

    run_case("detect_source", []() {
        StarCCMAdapter adapter;
        cfdx::io::SourceInfo info;
        bool detected = adapter.detect_source("tests/data/starccm/test_case.sim", info);
        EXPECT_TRUE(detected);
        EXPECT_TRUE(info.solver == "STAR-CCM+");
        EXPECT_TRUE(info.format == "sim");
    });

    run_case("convert_blocks_on_binary_topology", []() {
        StarCCMAdapter adapter;
        cfdx::io::ConversionResult result;
        bool ok = adapter.convert("tests/data/starccm/test_case.sim", result);
        EXPECT_FALSE(ok);
        EXPECT_TRUE(result.gap_report.has_blocking());
    });

    run_case("zones_detected", []() {
        StarCCMAdapter adapter;
        cfdx::io::ConversionResult result;
        adapter.convert("tests/data/starccm/test_case.sim", result);
        EXPECT_TRUE(adapter.zones().size() > 0);
    });

    run_case("import_results_unsupported", []() {
        StarCCMAdapter adapter;
        cfdx::io::ConversionResult result;
        bool ok = adapter.import_results("tests/data/starccm/test_case.sim", result);
        EXPECT_FALSE(ok);
        EXPECT_TRUE(result.gap_report.has_blocking());
    });

    return run_all();
}
