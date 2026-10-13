#include "cfdx/core/linalg/amg_preconditioner.h"
#include "cfdx/core/linalg/amg_setup_footprint.h"
#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/linear_operator.h"
#include "common/test_harness.h"

#include <cstddef>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

// N11 B4/B5 AMG setup footprint and native AMG contract checks (#480):
// - B4: MatrixFreeVcyclePreconditioner reports, through setup_footprint(),
//   how many singleton smoothed-aggregation aggregates, fallback
//   interpolation rows and dropped Galerkin entries the latest setup()
//   produced. The counters are reporting-only: no tolerance or acceptance
//   criterion reads them, and a failed setup() keeps its counters visible.
// - B5: NativeAMGPreconditioner documents the real last_error() contract
//   (apply() is const but mutates the message through the pimpl) and the
//   AMGSetupFootprint type is forwarded into cfdx::core::linalg.

static SparseMatrix make_1d_laplacian(std::size_t n) {
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A.push_back(i, i, 2.0);
        if (i > 0) A.push_back(i, i - 1, -1.0);
        if (i + 1 < n) A.push_back(i, i + 1, -1.0);
    }
    A.finalize();
    return A;
}

static FunctionalLinearOperator make_operator(const SparseMatrix& A) {
    return FunctionalLinearOperator(
        A.n_rows(), [&A](const Vector& x, Vector& y) {
            const auto values = A.matvec(x);
            for (std::size_t i = 0; i < values.size(); ++i) y(i) = values[i];
        });
}

static void expect_footprint(const MatrixFreeVcyclePreconditioner& amg,
                             std::size_t singletons,
                             std::size_t fallbacks,
                             std::size_t drops) {
    const AMGSetupFootprint& f = amg.setup_footprint();
    EXPECT_TRUE(f.singleton_aggregates == singletons);
    EXPECT_TRUE(f.fallback_interpolation_rows == fallbacks);
    EXPECT_TRUE(f.galerkin_dropped_entries == drops);
}

int main() {
    run_case("b4_directcf_clean_setup_reports_zero_footprint", [] {
        // Direct C/F interpolation on a clean 1-D Laplacian: Direct C/F
        // never builds smoothed-aggregation aggregates, every Galerkin
        // entry is O(1), and the footprint must read exactly zero.
        const SparseMatrix A = make_1d_laplacian(32);
        const auto op = make_operator(A);
        MatrixFreeVcyclePreconditioner amg(
            op, 0.7, 4, 4, 0.25, 25, AMGInterpolationPolicy::DirectCF);
        EXPECT_TRUE(amg.setup(A));
        expect_footprint(amg, 0, 0, 0);
        EXPECT_TRUE(amg.coarse_size() > 0);
    });

    run_case("b4_sa_counts_isolated_point_as_singleton", [] {
        // 1-D Laplacian on 17 coupled points plus one fully isolated row:
        // uncoupled aggregation must isolate the point into its own
        // aggregate, visible as exactly one singleton. The coupled part
        // aggregates into triples plus one final pair, so no other
        // singleton can appear.
        SparseMatrix A(18, 18);
        for (std::size_t i = 0; i < 17; ++i) {
            A.push_back(i, i, 2.0);
            if (i > 0) A.push_back(i, i - 1, -1.0);
            if (i + 1 < 17) A.push_back(i, i + 1, -1.0);
        }
        A.push_back(17, 17, 1.0);
        A.finalize();
        const auto op = make_operator(A);
        MatrixFreeVcyclePreconditioner amg(
            op, 0.7, 4, 4, 0.25, 25, AMGInterpolationPolicy::SmoothedAggregation);
        EXPECT_TRUE(amg.setup(A));
        expect_footprint(amg, 1, 0, 0);
        // Cross-check through the public first-level aggregation: the
        // aggregate id of the isolated point must have exactly one member.
        const std::size_t id = amg.aggregate_of(17);
        std::size_t members = 0;
        for (std::size_t i = 0; i < 18; ++i)
            if (amg.aggregate_of(i) == id) ++members;
        EXPECT_TRUE(members == 1);
    });

    run_case("b4_sa_counts_normalized_zero_row_fallback", [] {
        // Rows 0..2 reproduce a 3x3 pattern whose middle row sums its
        // smoothed weights to ~6e-17: below the 1e-14 cut the weight is
        // excluded from the row sum, the sum reads zero and the row falls
        // back to {own aggregate, 1.0}. Rows 3..16 are a clean Laplacian
        // chain, so this fallback is the only one and the setup stays valid.
        SparseMatrix A(17, 17);
        A.push_back(0, 0, 2.0); A.push_back(0, 1, -1.0);
        A.push_back(1, 0, 1.0); A.push_back(1, 1, 2.0);
        A.push_back(2, 1, -1.0); A.push_back(2, 2, 2.0);
        for (std::size_t i = 3; i < 17; ++i) {
            A.push_back(i, i, 2.0);
            if (i > 3) A.push_back(i, i - 1, -1.0);
            if (i + 1 < 17) A.push_back(i, i + 1, -1.0);
        }
        A.finalize();
        const auto op = make_operator(A);
        MatrixFreeVcyclePreconditioner amg(
            op, 0.7, 4, 4, 0.25, 25, AMGInterpolationPolicy::SmoothedAggregation);
        EXPECT_TRUE(amg.setup(A));
        expect_footprint(amg, 0, 1, 0);
    });

    run_case("b4_sa_counts_dropped_zero_galerkin_entry_and_fails", [] {
        // Rows 0..1 form [[2,-2],[-2,2]]: both rows share one aggregate and
        // normalize to exactly 1.0, so the Galerkin product of the block is
        // exactly 0.0. The single coarse entry is dropped, the coarse row
        // stays empty and setup() must fail closed while keeping the
        // footprint visible. Rows 2..16 keep the matrix large enough for
        // the hierarchy loop to run (n > 16).
        SparseMatrix A(17, 17);
        A.push_back(0, 0, 2.0); A.push_back(0, 1, -2.0);
        A.push_back(1, 0, -2.0); A.push_back(1, 1, 2.0);
        for (std::size_t i = 2; i < 17; ++i) {
            A.push_back(i, i, 2.0);
            if (i > 2) A.push_back(i, i - 1, -1.0);
            if (i + 1 < 17) A.push_back(i, i + 1, -1.0);
        }
        A.finalize();
        const auto op = make_operator(A);
        MatrixFreeVcyclePreconditioner amg(
            op, 0.7, 4, 4, 0.25, 25, AMGInterpolationPolicy::SmoothedAggregation);
        EXPECT_TRUE(!amg.setup(A));
        expect_footprint(amg, 0, 0, 1);
        // Failed setup leaves the object safe: apply() must refuse.
        Vector r(17, 1.0), z(17, 0.0);
        EXPECT_TRUE(!amg.apply(r, z));
    });

    run_case("b4_new_setup_resets_previous_footprint", [] {
        // The same SA instance: a failed setup with one Galerkin drop,
        // then a clean setup on a valid Laplacian must reset the counters
        // rather than accumulate. setup() validates the matrix against
        // the operator dimensions only, so the operator built from the
        // clean matrix also accepts the bad 17x17 matrix for the first
        // (failing) attempt; apply() is never called on it.
        SparseMatrix bad(17, 17);
        bad.push_back(0, 0, 2.0); bad.push_back(0, 1, -2.0);
        bad.push_back(1, 0, -2.0); bad.push_back(1, 1, 2.0);
        for (std::size_t i = 2; i < 17; ++i) {
            bad.push_back(i, i, 2.0);
            if (i > 2) bad.push_back(i, i - 1, -1.0);
            if (i + 1 < 17) bad.push_back(i, i + 1, -1.0);
        }
        bad.finalize();
        const SparseMatrix clean = make_1d_laplacian(17);
        const auto op = make_operator(clean);
        MatrixFreeVcyclePreconditioner amg(
            op, 0.7, 4, 4, 0.25, 25, AMGInterpolationPolicy::SmoothedAggregation);
        EXPECT_TRUE(!amg.setup(bad));
        EXPECT_TRUE(amg.setup_footprint().galerkin_dropped_entries == 1);
        EXPECT_TRUE(amg.setup(clean));
        expect_footprint(amg, 0, 0, 0);
        EXPECT_TRUE(amg.coarse_size() > 0);
    });

    run_case("b4_numeric_update_commits_only_successful_diagnostics", [] {
        const SparseMatrix clean = make_1d_laplacian(32);
        const auto op = make_operator(clean);
        MatrixFreeVcyclePreconditioner amg(
            op, 0.7, 4, 4, 0.25, 25, AMGInterpolationPolicy::SmoothedAggregation);
        EXPECT_TRUE(amg.setup(clean));
        const AMGSetupFootprint before = amg.setup_footprint();

        // A valid numeric refresh commits its hierarchy and keeps the
        // diagnostic stable when no Galerkin entries are dropped.
        EXPECT_TRUE(amg.update_values(clean));
        expect_footprint(amg, before.singleton_aggregates,
                         before.fallback_interpolation_rows,
                         before.galerkin_dropped_entries);

        // An invalid diagonal fails before a refreshed hierarchy is
        // committed. Its temporary attempt must not change the active
        // hierarchy's diagnostic counters.
        SparseMatrix invalid = clean;
        auto* values = invalid.values_data();
        for (std::size_t i = 0; i < invalid.nnz(); ++i) values[i] = 0.0;
        EXPECT_TRUE(!amg.update_values(invalid));
        expect_footprint(amg, before.singleton_aggregates,
                         before.fallback_interpolation_rows,
                         before.galerkin_dropped_entries);
        EXPECT_TRUE(amg.coarse_size() > 0);
    });

    run_case("b5_native_amg_last_error_contract_and_forwarding", [] {
        // last_error() is the status of the most recent lifecycle
        // operation: empty after construction, set by a failed
        // setup()/apply(), and apply() is const but still mutates the
        // message through the pimpl.
        NativeBoomerAMGPreconditioner p;
        EXPECT_TRUE(!p.is_ready());
        EXPECT_TRUE(p.last_error().empty());

        SparseMatrix empty(0, 0);
        empty.finalize();
        EXPECT_TRUE(!p.setup(empty));
        EXPECT_TRUE(!p.last_error().empty());

        Vector r(32, 1.0), z(32, 0.0);
        EXPECT_TRUE(!p.apply(r, z));
        EXPECT_TRUE(!p.last_error().empty());

        // A successful lifecycle operation clears the previous diagnostic.
        const SparseMatrix valid = make_1d_laplacian(32);
        EXPECT_TRUE(p.setup(valid));
        EXPECT_TRUE(p.last_error().empty());
        EXPECT_TRUE(p.update_values(valid));
        EXPECT_TRUE(p.last_error().empty());
        EXPECT_TRUE(p.apply(r, z));
        EXPECT_TRUE(p.last_error().empty());

        // Type forwarding: the footprint is reachable from the linalg
        // namespace used by the solver-facing headers.
        cfdx::core::linalg::AMGSetupFootprint f{};
        EXPECT_TRUE(f.singleton_aggregates == 0);
        EXPECT_TRUE(f.fallback_interpolation_rows == 0);
        EXPECT_TRUE(f.galerkin_dropped_entries == 0);
    });

    return run_all();
}
