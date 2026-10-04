#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/vector.h"

#include <string>

namespace cfdx::core {

// What apply() computes.
//
// Because the interface deliberately separates the Schur approximation from the
// outer block factorisation and from the pressure Krylov solver, it admits both
// kinds of action. Declaring which one an implementation provides is part of the
// contract rather than a hint: the two differ by a solve, so a caller that needs
// an inverse action and silently receives the operator has applied a
// mathematically different map with no way to detect it from the return value.
enum class SchurAction {
    // apply() returns S~ * rhs_p -- the approximated Schur operator itself.
    Operator,
    // apply() returns S~^{-1} * rhs_p -- one application of its inverse.
    InverseOperator,
};

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

    // The action apply() performs. Mandatory, so that a caller needing one
    // specific action can reject a mismatched implementation instead of
    // assuming it.
    virtual SchurAction action() const noexcept = 0;

    // Apply the configured approximation to a pressure right-hand side, in the
    // sense declared by action(). `out` receives the pressure-side result.
    virtual bool apply(const Vector& rhs_p, Vector& out) const = 0;

    // Refresh numerical values when the block sparsity patterns are unchanged.
    // A changed graph must be rejected and handled through setup().
    virtual bool update_values(const BlockOperator& blocks) = 0;
};

// True when the approximation provides exactly the requested action. A caller
// implementing a preconditioner needs an inverse action; a caller assembling a
// residual or an oracle comparison needs the operator.
inline bool provides_action(const SchurApproximation& approximation,
                            SchurAction required) noexcept {
    return approximation.action() == required;
}

} // namespace cfdx::core
