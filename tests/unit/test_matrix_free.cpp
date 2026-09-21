#include "cfdx/core/numerics/matrix_free.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

int main()
{
    run_case("matrix_free_matches_sparse_matvec", [] {
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 4.0);
        A.push_back(0, 1, -1.0);
        A.push_back(1, 0, -1.0);
        A.push_back(1, 1, 4.0);
        A.push_back(1, 2, -1.0);
        A.push_back(2, 1, -1.0);
        A.push_back(2, 2, 4.0);
        A.finalize();

        const auto op = matrix_free_from_sparse(A);
        Vector x(3);
        x(0)=1.0; x(1)=2.0; x(2)=3.0;
        const auto y_sparse = A.matvec(x);
        const auto y_free = op(x);
        for (std::size_t i=0; i<3; ++i)
            EXPECT_NEAR(y_free(i), y_sparse[i], 1e-14);
    });

    run_case("matrix_free_rejects_dimension_mismatch", [] {
        MatrixFreeOperator op(2, 3, [](const Vector&, Vector& y) { y(0)=0.0; y(1)=0.0; });
        Vector x(2);
        EXPECT_THROW(op(x), std::invalid_argument);
    });

    return run_all();
}
