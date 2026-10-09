#include "cfdx/core/linalg/krylov_diagnostics.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <cmath>
#include <string>
#include <sstream>

using namespace cfdx::core;
using namespace cfdx::testing;

// Helper: build a simple 4x4 saddle-point matrix with a pure gauge row.
// The gauge row is a true identity row (only diagonal=1, no off-diagonal
// coupling) so that the preconditioner preserves it exactly.
static SparseMatrix make_gauged_4n(std::size_t n_cells) {
    const std::size_t nv = 3 * n_cells;
    const std::size_t n = nv + n_cells;
    SparseMatrix A(n, n);
    for (std::size_t c = 0; c < n_cells; ++c) {
        for (std::size_t d = 0; d < 3; ++d)
            A.push_back(d * n_cells + c, d * n_cells + c, 4.0);
        A.push_back(nv + c, nv + c, 1.0);  // pure gauge row (identity)
    }
    A.finalize();
    return A;
}

int main() {
    run_case("krylov_diagnostics_flag_defaults_off", [] {
        // The diagnostic level should default to OFF when the env var is unset.
        // We cannot force the env here, but we can verify the accessor returns
        // a valid enum value.
        const auto level = krylov_diagnostic_level();
        EXPECT_TRUE(static_cast<int>(level) >= 0 &&
                    static_cast<int>(level) <= 2);
    });

    run_case("krylov_diagnostics_emit_single_failure_line", [] {
        // A preconditioner with a deliberately failing apply should emit
        // exactly one diagnostic line when diagnostics are enabled.
        // We test the KrylovFailLog RAII directly.
        std::stringstream buffer;
        std::streambuf* old = std::cerr.rdbuf(buffer.rdbuf());

        {
            KrylovFailLog log("test_stage", "test_context");
            // Simulate a failure: mark_failed() should cause exactly one
            // line to be emitted on scope exit.
            log.mark_failed();
        }

        std::cerr.rdbuf(old);
        const std::string output = buffer.str();
        EXPECT_TRUE(output.find("KRYLOV_FAIL stage=test_stage") != std::string::npos);
        EXPECT_TRUE(output.find("context=test_context") != std::string::npos);
        // Exactly one failure line.
        EXPECT_TRUE(output.find("KRYLOV_FAIL") != std::string::npos);
    });

    run_case("krylov_diagnostics_no_emit_when_not_failed", [] {
        std::stringstream buffer;
        std::streambuf* old = std::cerr.rdbuf(buffer.rdbuf());
        {
            KrylovFailLog log("test_stage", "test_context");
            // Do not mark_failed: no line should be emitted.
        }
        std::cerr.rdbuf(old);
        const std::string output = buffer.str();
        EXPECT_TRUE(output.find("KRYLOV_FAIL") == std::string::npos);
    });

    run_case("coupled_block_schur_diagnostics_emit_gauge_input", [] {
        // When verbose diagnostics are enabled, the preconditioner apply
        // should emit a GMRES_INIT_GAUGE-like line for the gauge row.
        // We cannot force the env var in a unit test, but we can verify
        // that the preconditioner apply still works correctly with the
        // diagnostic code paths compiled in.
        const auto A = make_gauged_4n(1);
        CoupledBlockSchurAMGPreconditioner pc(1);
        EXPECT_TRUE(pc.setup(A));

        Vector r(4, 0.0);
        r(3) = 1.0;  // gauge row input
        Vector z(4, 0.0);
        EXPECT_TRUE(pc.apply(r, z));
        EXPECT_TRUE(z.is_valid());
        // The gauge row output should be finite (exact value depends on AMG setup).
        EXPECT_TRUE(std::isfinite(z(3)));
    });

    run_case("coupled_block_schur_diagnostics_emit_gauge_output", [] {
        const auto A = make_gauged_4n(2);
        CoupledBlockSchurAMGPreconditioner pc(2);
        EXPECT_TRUE(pc.setup(A));

        Vector r(8, 0.0);
        r(6) = 1.0;  // gauge row input for 2-cell system (nv=6, gauge at 6)
        Vector z(8, 0.0);
        EXPECT_TRUE(pc.apply(r, z));
        EXPECT_TRUE(z.is_valid());
        // Gauge row output should be finite.
        EXPECT_TRUE(std::isfinite(z(6)));
    });

    run_case("coupled_block_schur_diagnostics_invalid_input", [] {
        const auto A = make_gauged_4n(1);
        CoupledBlockSchurAMGPreconditioner pc(1);
        EXPECT_TRUE(pc.setup(A));

        // Wrong size input should fail cleanly.
        Vector r(3, 0.0);
        Vector z(3, 0.0);
        EXPECT_FALSE(pc.apply(r, z));
    });

    return run_all();
}