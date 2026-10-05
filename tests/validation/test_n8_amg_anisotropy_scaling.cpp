#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

class TestSparseOperator final : public LinearOperatorBase {
public:
    explicit TestSparseOperator(const SparseMatrix& A) : A_(A) {}
    std::size_t rows() const noexcept override { return A_.n_rows(); }