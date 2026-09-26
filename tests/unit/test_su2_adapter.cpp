// M0.10-T03: Tests for SU2 adapter
#include "cfdx/io/su2/su2_adapter.h"
#include "cfdx/io/cfdx_io/io_interface.h"
#include "cfdx/io/cfdx_io/case_schema.h"
#include "common/test_harness.h"

#include <string>

using namespace cfdx::io::su2;
using namespace cfdx::testing;

int main() {
    run_case("adapter_metadata", []() {
        Su2Adapter adapter;
        EXPECT_TRUE(std::string(adapter.solver_name()) == "SU2");
        EXPECT_TRUE(std::string(adapter.format_name()) == "su2");
    });

    run_case("detect_source", []() {
        Su2Adapter adapter;
        cfdx::io::SourceInfo info;
        bool detected = adapter.detect_source("tests/data/su2/mesh_NACA0012_inv.su2", info);
        EXPECT_TRUE(detected);
        EXPECT_TRUE(info.solver == "SU2");
        EXPECT_TRUE(info.format == "su2");
    });

    run_case("detect_source_not_found", []() {
        Su2Adapter adapter;
        cfdx::io::SourceInfo info;
        bool detected = adapter.detect_source("tests/data/su2/nonexistent.su2", info);
        EXPECT_FALSE(detected);
    });

    run_case("parse_mesh", []() {
        Su2Adapter adapter;
        bool ok = adapter.parse_mesh("tests/data/su2/mesh_NACA0012_inv.su2");
        EXPECT_TRUE(ok);
        EXPECT_TRUE(adapter.su2_points().size() == 5233);
        EXPECT_TRUE(adapter.su2_boundaries().size() == 2);

        const auto& bounds = adapter.su2_boundaries();
        EXPECT_TRUE(bounds[0].name == "farfield" || bounds[1].name == "farfield");
        EXPECT_TRUE(bounds[0].name == "airfoil" || bounds[1].name == "airfoil");

        const auto& setup = adapter.setup();
        EXPECT_TRUE(setup.mesh_info.n_vertices == 5233);
        EXPECT_TRUE(setup.mesh_info.dimension == 2);
        EXPECT_TRUE(setup.mesh_info.n_cells == 10216);
    });

    run_case("parse_config", []() {
        Su2Adapter adapter;
        bool ok = adapter.parse_config("tests/data/su2/inv_NACA0012_basic.cfg");
        EXPECT_TRUE(ok);
        const auto& setup = adapter.setup();
        EXPECT_TRUE(setup.physics_model != "");
    });

    run_case("convert_full", []() {
        Su2Adapter adapter;
        cfdx::io::ConversionResult result;
        bool ok = adapter.convert("tests/data/su2/mesh_NACA0012_inv.su2", result);
        EXPECT_TRUE(ok);
        EXPECT_FALSE(result.gap_report.has_blocking());
        EXPECT_TRUE(result.mesh.n_vertices() == 5233);
        EXPECT_TRUE(result.setup.mesh_info.n_vertices == 5233);
    });

    return run_all();
}
