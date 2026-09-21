#include "cfdx/core/linalg/amg_preconditioner.h"
#include "cfdx/core/linalg/chebyshev_smoother.h"
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/linalg/sparse_matrix.h"

#include <cmath>
#include <cstddef>
#include <iostream>

int main() {
    using namespace cfdx::core;

    SparseMatrix A(4, 4);
    A.push_back(0, 0, 2.0);
    A.push_back(0, 1, -1.0);
    A.push_back(1, 0, -1.0);
    A.push_back(1, 1, 2.0);
    A.push_back(1, 2, -1.0);
    A.push_back(2, 1, -1.0);
    A.push_back(2, 2, 2.0);
    A.push_back(2, 3, -1.0);
    A.push_back(3, 2, -1.0);
    A.push_back(3, 3, 2.0);
    A.finalize();

    FunctionalLinearOperator op(
        4,
        [&A](const Vector& x, Vector& y) {
            const auto result = A.matvec(x);
            for (std::size_t i = 0; i < result.size(); ++i) {
                y(i) = result[i];
            }
        });

    MatrixFreeVcyclePreconditioner amg(op);
    if (!amg.setup(A)) {
        return 1;
    }

    if (amg.coarse_size() != 2 ||
        amg.aggregate_of(0) != amg.aggregate_of(1) ||
        amg.aggregate_of(2) != amg.aggregate_of(3) ||
        amg.aggregate_of(0) == amg.aggregate_of(2)) {
        return 2;
    }

    Vector rhs(4, 1.0);
    Vector correction;
    if (!amg.apply(rhs, correction) || !correction.is_valid()) {
        return 3;
    }

    ChebyshevSmoother::Controls controls;
    controls.lambda_max = 0.0;
    controls.spectral_iterations = 4;

    ChebyshevSmoother smoother(controls);
    const double lambda_max = smoother.lambda_max(op);

    if (!(lambda_max > 0.0) || !std::isfinite(lambda_max)) {
        return 4;
    }

    Vector x(4, 0.0);
    smoother.apply(op, rhs, x);

    if (!x.is_valid()) {
        return 5;
    }

    std::cout << "AMG/Chebyshev runtime: PASS\n";
    return 0;
}
