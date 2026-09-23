#pragma once

#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"

#include <cstddef>

namespace cfdx::runtime::gpu {

bool cuda_poisson_available();

cfdx::core::SolverResult solve_poisson_cuda(
    const cfdx::core::SparseMatrix& matrix,
    const cfdx::core::Vector& rhs,
    cfdx::core::Vector& solution,
    std::size_t max_iterations = 2000,
    double tolerance = 1e-10);

} // namespace cfdx::runtime::gpu
