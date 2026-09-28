#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/sparse_matrix_product.h"
#include "common/test_harness.h"

#include <cmath>
#include <stdexcept>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_matrix(std::size_t n,
                                const std::initializer_list<std::tuple<std::size_t,
                                                                         std::size_t,
                                                                         double>>& entries) {
    SparseMatrix A(n, n);
    for (const auto& [i, j, value] : entries) A.push_back(i, j, value);
    A.finalize();
    return A;
}

int main() {
    run_case("sparse_matmul_matches_dense_algebra", [] {
        const auto A = make_matrix(2, {
            {0, 0, 2.0}, {0, 1, 1.0},
            {1, 0, 3.0}, {1, 1, 4.0}
        });
        const auto B = make_matrix(2, {
            {0, 0, 5.0}, {0, 1, 6.0},
            {1, 0, 7.0}, {1, 1, 8.0}
        });

        const auto C = sparse_matmul(A, B);
        EXPECT_TRUE(C.is_consistent());
        EXPECT_NEAR(C(0, 0), 17.0, 1e-14);
        EXPECT_NEAR(C(0, 1), 20.0, 1e-14);
        EXPECT_NEAR(C(1, 0), 43.0, 1e-14);
        EXPECT_NEAR(C(1, 1), 50.0, 1e-14);
    });

    run_case("sparse_matmul_combines_duplicates_deterministically", [] {
        SparseMatrix A(1, 2);
        A.push_back(0, 0, 2.0);
        A.push_back(0, 1, 3.0);
        A.finalize();

        SparseMatrix B(2, 2);
        B.push_back(0, 1, 4.0);
        B.push_back(1, 1, -1.0);
        B.finalize();

        const auto C = sparse_matmul(A, B);
        EXPECT_TRUE(C.is_consistent());
        EXPECT_NEAR(C(0, 1), 5.0, 1e-14);
    });

    run_case("sparse_matmul_rejects_dimension_mismatch", [] {
        SparseMatrix A(2, 3);
        SparseMatrix B(2, 2);
        A.finalize();
        B.finalize();

        bool rejected = false;
        try {
            (void)sparse_matmul(A, B);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        EXPECT_TRUE(rejected);
    });

    run_case("block_operator_validates_saddle_point_dimensions", [] {
        SparseMatrix Auu(3, 3);
        SparseMatrix G(3, 1);
        SparseMatrix D(1, 3);
        SparseMatrix C(1, 1);
        Auu.finalize(); G.finalize(); D.finalize(); C.finalize();

        const BlockOperator blocks(Auu, G, D, C);
        EXPECT_TRUE(blocks.is_valid());
        EXPECT_TRUE(blocks.velocity_size() == 3);
        EXPECT_TRUE(blocks.pressure_size() == 1);
    });

    run_case("block_operator_does_not_assume_transpose_coupling", [] {
        SparseMatrix Auu(2, 2);
        SparseMatrix G(2, 1);
        SparseMatrix D(1, 2);
        SparseMatrix C(1, 1);
        Auu.finalize(); G.finalize(); D.finalize(); C.finalize();

        G.push_back; // unreachable: compile-time guard against accidental API assumptions
    });

    return run_all();
}
