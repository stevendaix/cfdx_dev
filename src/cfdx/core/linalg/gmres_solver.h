#pragma once
#include <vector>
#include "sparse_matrix.h"
#include "vector.h"

namespace cfdx::core::linalg {

bool gmres_solve(const SparseMatrix& A,
                 const std::vector<double>& b,
                 std::vector<double>& x,
                 int restart = 30,
                 double tolerance = 1e-6,
                 int max_iter = 100);
} // namespace
