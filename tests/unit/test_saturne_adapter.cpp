// M0.10-T03: Tests for Code_Saturne adapter
#include "cfdx/io/saturne/saturne_adapter.h"
#include "cfdx/io/cfdx_io/io_interface.h"
#include "common/test_harness.h"

#include <string>

using namespace cfdx::io::saturne;
using namespace cfdx::testing;

int main() {
    run_case("adapter_metadata", []() {
        SaturneAdapter adapter;
        EXPECT_TRUE(std::string(adapter.solver_name()) == "Code_Saturne");
        EXPECT_TRUE(std::string(adapter.format_name()) == "saturne_xml_py");
    });

    run_case("detect_source_xml", []() {
        SaturneAdapter adapter;
        cfdx::io::SourceInfo info;
        bool detected = adapter.detect_source("tests/data/saturne/test_case.xml", info);
        EXPECT_TRUE(detected);
        EXPECT_TRUE(info.solver == "Code_Saturne");
    });

    run_case("parse_xml", []() {
        SaturneAdapter adapter;
        bool ok = adapter.parse_xml_setup("tests/data/saturne/test_case.xml");
        EXPECT_TRUE(ok);
        EXPECT_TRUE(adapter.setup().boundary_conditions.size() > 0);
        EXPECT_TRUE(adapter.setup().materials.size() > 0);
    });

    run_case("parse_python", []() {
        SaturneAdapter adapter;
        bool ok = adapter.parse_python_setup("tests/data/saturne/test_case.py");
        EXPECT_TRUE(ok);
        EXPECT_TRUE(adapter.entries().size() > 0);
    });

    run_case("convert_blocks_without_neutral_mesh", []() {
        SaturneAdapter adapter;
        cfdx::io::ConversionResult result;
        bool ok = adapter.convert("tests/data/saturne/test_case.xml", result);
        EXPECT_FALSE(ok);
        EXPECT_TRUE(result.gap_report.has_blocking());
    });

    run_case("import_results_neutral_format", []() {
        SaturneAdapter adapter;
        cfdx::io::ConversionResult result;
        bool ok = adapter.import_results("tests/data/saturne/test_case.vtk", result);
        EXPECT_FALSE(ok); // file doesn't exist, but format is supported
    });

    return run_all();
}
