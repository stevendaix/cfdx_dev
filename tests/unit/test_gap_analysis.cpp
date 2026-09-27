// M0.10-T03: Tests for GapAnalysis engine
#include "cfdx/io/cfdx_io/io_interface.h"
#include "common/test_harness.h"

using namespace cfdx::io;
using namespace cfdx::testing;

int main() {
    run_case("default_empty", []() {
        GapAnalysis gap;
        EXPECT_FALSE(gap.has_blocking());
        EXPECT_TRUE(gap.findings().empty());
        EXPECT_TRUE(gap.n_supported() == 0);
    });

    run_case("add_supported", []() {
        GapAnalysis gap;
        gap.supported("mesh", "points", "Parsed 100 vertices");
        EXPECT_FALSE(gap.has_blocking());
        EXPECT_TRUE(gap.n_supported() == 1);
    });

    run_case("add_blocking", []() {
        GapAnalysis gap;
        gap.unsupported_blocking("mesh", "topology", "Not supported", "Use export");
        EXPECT_TRUE(gap.has_blocking());
        EXPECT_TRUE(gap.n_unsupported_blocking() == 1);
    });

    run_case("add_nonblocking", []() {
        GapAnalysis gap;
        gap.unsupported_nonblocking("features", "x", "Not supported", "Workaround");
        EXPECT_FALSE(gap.has_blocking());
        EXPECT_TRUE(gap.n_unsupported_nonblocking() == 1);
    });

    run_case("add_approximated", []() {
        GapAnalysis gap;
        gap.approximated("physics", "model", "Approximated", "Verify manually");
        EXPECT_FALSE(gap.has_blocking());
        EXPECT_TRUE(gap.n_approximated() == 1);
        // Approximated counts as supported for n_supported
        EXPECT_TRUE(gap.n_supported() == 1);
    });

    run_case("add_unavailable", []() {
        GapAnalysis gap;
        gap.unavailable("config", "file", "File not found");
        EXPECT_FALSE(gap.has_blocking());
        EXPECT_TRUE(gap.n_unavailable() == 1);
    });

    run_case("mixed_findings", []() {
        GapAnalysis gap;
        gap.supported("mesh", "points", "Parsed");
        gap.approximated("physics", "model", "Approx", "Fix");
        gap.unsupported_nonblocking("features", "x", "Not supported", "Fix");
        gap.unsupported_blocking("mesh", "topo", "Blocked", "Export");
        EXPECT_TRUE(gap.has_blocking());
        EXPECT_TRUE(gap.n_supported() == 2);  // supported + approximated
        EXPECT_TRUE(gap.n_approximated() == 1);
        EXPECT_TRUE(gap.n_unsupported_nonblocking() == 1);
        EXPECT_TRUE(gap.n_unsupported_blocking() == 1);
        EXPECT_TRUE(gap.n_unavailable() == 0);
    });

    return run_all();
}
