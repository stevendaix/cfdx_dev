#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/vector.h"

#include <string>

namespace cfdx::core {

// Common contract for algebraic Schur approximations:
//
//     S = C - D Auu^{-1} G
//
// Implementations are responsible for defining the approximation used for
// Auu^{-1}, the pressure-side operator, and its inner/preconditioning solve.
// The interface intentionally separates the Schur approximation from the
// outer block factorisation and from the pressure Krylov solver.
class SchurApproximation {
public:
    virtual ~SchurApproximation() = default;

    virtual const char* name() const noexcept = 0;

    // Build all value-dependent state from the current block operator.
    // Implementations must reject an invalid block layout.
    virtual bool setup(const BlockOperator& blocks) = 0;

    // Apply the configured approximation to a pressure right-hand side.
    virtual bool apply(const Vector& rhs_p, Vector& pressure) const = 0;

    // Refresh numerical values when the block sparsity patterns are unchanged.
    // A changed graph must be rejected and handled through setup().
    virtual bool update_values(const BlockOperator& blocks) = 0;
};

} // namespace cfdx::core
