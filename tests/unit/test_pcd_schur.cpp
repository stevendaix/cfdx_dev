#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/pcd_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <cmath>
#include <initializer_list>
#include <tuple>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_sparse(
    std::size_t rows, std::size_t cols,
    const std::initializer_list<std::tuple<std::size_t, std::size_t, double>>& entries) {
    SparseMatrix A(rows, cols);
    for (const auto& [i, j, value] : entries) A.push_back(i, j, value);
    A.finalize();
    return A;
}

static std::vector<std::vector<double>> dense(const SparseMatrix& A) {
    std::vector<std::vector<double>> out(A.n_rows(), std::vector<double>(A.n_cols(), 0.0));
    for (std::size_t i = 0; i < A.n_rows(); ++i)
        for (std::size_t k = A.row_offsets_data()[i]; k < A.row_offsets_data()[i + 1]; ++k)
            out[i][A.columns_data()[k]] += A.values_data()[k];
    return out;
}

static std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> A) {
    const std::size_t n = A.size();
    std::vector<std::vector<double>> I(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) I[i][i] = 1.0;

    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i)
            if (std::abs(A[i][k]) > std::abs(A[pivot][k])) pivot = i;
        if (std::abs(A[pivot][k]) <= 1e-14) throw std::runtime_error("singular test matrix");
        std::swap(A[k], A[pivot]);
        std::swap(I[k], I[pivot]);

        const double d = A[k][k];
        for (std::size_t j = 0; j < n; ++j) {
            A[k][j] /= d;
            I[k][j] /= d;
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            const double f = A[i][k];
            for (std::size_t j = 0; j < n; ++j) {
                A[i][j] -= f * A[k][j];
                I[i][j] -= f * I[k][j];
            }
        }
    }
    return I;
}

static Vector dense_matvec(const std::vector<std::vector<double>>& A, const Vector& x) {
    Vector y(A.size(), 0.0);
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j)
            y(i) += A[i][j] * x(j);
    return y;
}

static void assert_close(const Vector& a, const Vector& b, double tol = 1e-12) {
    EXPECT_TRUE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i)
        EXPECT_TRUE(std::abs(a(i) - b(i)) <= tol);
}

int main() {
    const auto Auu = make_sparse(2, 2, {
        {0, 0, 4.0}, {0, 1, 0.5}, {1, 0, -0.25}, {1, 1, 3.0}
    });
    const auto G = make_sparse(2, 2, {{0, 0, 1.0}, {1, 1, 1.0}});
    const auto D = make_sparse(2, 2, {{0, 0, 1.0}, {1, 1, 1.0}});
    const auto C = make_sparse(2, 2, {});

    const auto Mp = make_sparse(2, 2, {
        {0, 0, 2.0}, {0, 1, 0.2}, {1, 0, 0.1}, {1, 1, 1.5}
    });
    const auto Kp = make_sparse(2, 2, {
        {0, 0, 3.0}, {0, 1, -0.4}, {1, 0, 0.2}, {1, 1, 2.0}
    });
    const auto Fp = make_sparse(2, 2, {
        {0, 0, 1.7}, {0, 1, 0.3}, {1, 0, -0.15}, {1, 1, 2.4}
    });

    const BlockOperator blocks(Auu, G, D, C);
    blocks.validate();

    const auto Kinv = inverse(dense(Kp));
    const auto Finv = inverse(dense(Fp));

    auto solve_K = [Kinv](const Vector& rhs, Vector& x) {
        if (rhs.size() != 2) return false;
        x = dense_matvec(Kinv, rhs);
        return true;
    };
    auto solve_F = [Finv](const Vector& rhs, Vector& x) {
        if (rhs.size() != 2) return false;
        x = dense_matvec(Finv, rhs);
        return true;
    };

    run_case("pcd_setup_and_action", [&] {
        PcdSchurApproximation pcd(Mp, Kp, Fp, solve_K, solve_F);
        EXPECT_TRUE(pcd.setup(blocks));
        Vector rhs(2, 0.0);
        rhs(0) = 1.0;
        rhs(1) = -0.7;
        Vector actual(2, 0.0);
        EXPECT_TRUE(pcd.apply(rhs, actual));

        const Vector z = dense_matvec(Kinv, rhs);
        const Vector y = dense(Mp).size() ? dense_matvec(dense(Mp), z) : Vector{};
        Vector expected = dense_matvec(Finv, y);
        expected *= -1.0;
        assert_close(actual, expected);
    });

    run_case("pcd_reports_name", [&] {
        PcdSchurApproximation pcd(Mp, Kp, Fp, solve_K, solve_F);
        EXPECT_TRUE(std::string(pcd.name()) == "pcd_schur");
    });

    run_case("pcd_rejects_dimension_mismatch", [&] {
        const auto bad_F = make_sparse(3, 3, {{0, 0, 1.0}, {1, 1, 1.0}, {2, 2, 1.0}});
        PcdSchurApproximation pcd(Mp, Kp, bad_F, solve_K, solve_F);
        EXPECT_TRUE(!pcd.setup(blocks));
    });

    run_case("pcd_rejects_graph_change", [&] {
        PcdSchurApproximation pcd(Mp, Kp, Fp, solve_K, solve_F);
        EXPECT_TRUE(pcd.setup(blocks));
        const auto changed_K = make_sparse(2, 2, {
            {0, 0, 3.0}, {0, 1, -0.4}, {1, 1, 2.0}
        });
        const BlockOperator same_blocks(Auu, G, D, C);
        const PcdSchurApproximation changed(
            Mp, changed_K, Fp, solve_K, solve_F);
        EXPECT_TRUE(!changed.setup(same_blocks));
        EXPECT_TRUE(!pcd.update_values(same_blocks));
    });

    run_case("pcd_accepts_value_update_same_graph", [&] {
        PcdSchurApproximation pcd(Mp, Kp, Fp, solve_K, solve_F);
        EXPECT_TRUE(pcd.setup(blocks));
        auto updated_F = make_sparse(2, 2, {
            {0, 0, 1.9}, {0, 1, 0.3}, {1, 0, -0.15}, {1, 1, 2.2}
        });
        PcdSchurApproximation updated(Mp, Kp, updated_F, solve_K, solve_F);
        EXPECT_TRUE(updated.setup(same_blocks));
    });

    run_case("pcd_requires_pressure_nullspace_compatibility", [&] {
        const auto null_space = NullSpaceProjector::constant(2);
        PcdSchurApproximation pcd(Mp, Kp, Fp, solve_K, solve_F, null_space);
        EXPECT_TRUE(pcd.setup(blocks));

        Vector incompatible(2, 1.0);
        Vector out(2, 0.0);
        EXPECT_TRUE(!pcd.apply(incompatible, out));

        Vector compatible(2, 0.0);
        compatible(0) = 1.0;
        compatible(1) = -1.0;
        EXPECT_TRUE(pcd.apply(compatible, out));
        EXPECT_TRUE(std::abs(out(0) + out(1)) <= 1e-12);
    });
}
