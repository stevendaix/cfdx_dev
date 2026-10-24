#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <limits>

namespace cfdx {
namespace core {

// M0.8-T04 — Solver outcome enumeration.
//
// This header is intentionally lightweight (C++14-compatible, no heavy
// dependencies) so that CUDA translation units can reference the solver
// outcome type without pulling in the full CG/GMRES/preconditioner headers,
// which themselves rely on C++17 features that older CUDA toolchains (11.5
// with GCC 11 libstdc++) cannot compile in device code.
enum class SolverStatus : std::uint8_t {
    CONVERGED = 0,
    MAX_ITER_REACHED,
    DIVERGED,
    NOT_APPLICABLE
};

inline const char* to_string(SolverStatus s) {
    switch (s) {
        case SolverStatus::CONVERGED:       return "converged";
        case SolverStatus::MAX_ITER_REACHED: return "max_iter_reached";
        case SolverStatus::DIVERGED:        return "diverged";
        case SolverStatus::NOT_APPLICABLE:  return "not_applicable";
        default:                             return "unknown";
    }
}

struct SolverDiagnosticSample {
    std::size_t iteration = 0;
    double true_residual = 0.0;
    double recursive_residual = 0.0;
    double preconditioned_dot = 0.0;
    double pAp = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
};

struct SolverResult {
    SolverStatus status = SolverStatus::NOT_APPLICABLE;
    std::size_t iterations = 0;
    double residual = 0.0;
    double residual_relative = 0.0;
    // Diagnostic evidence sampled during the iteration. This is intentionally
    // observational: it does not alter stopping criteria or Krylov updates.
    std::vector<SolverDiagnosticSample> diagnostics;
    // Residual-replacement evidence is observational and does not alter
    // convergence criteria. It is used to distinguish Krylov recurrence drift
    // from preconditioner/coarse-space limitations.
    std::size_t residual_replacements = 0;
    double max_true_recursive_gap = 0.0;
    double min_true_recursive_ratio = std::numeric_limits<double>::infinity();
    double max_true_recursive_ratio = 0.0;
    double max_true_residual = 0.0;
    double min_true_residual = std::numeric_limits<double>::infinity();
    std::size_t first_residual_replacement = 0;
    std::size_t last_residual_replacement = 0;
    // For solvers that relax the returned field after solving the linear
    // system, distinguish the linear candidate residual from the residual
    // of the actual returned field. NaN means not measured/not applicable.
    double residual_before_relaxation = std::numeric_limits<double>::quiet_NaN();
    double residual_after_relaxation = std::numeric_limits<double>::quiet_NaN();
    double residual_relative_before_relaxation = std::numeric_limits<double>::quiet_NaN();
    double residual_relative_after_relaxation = std::numeric_limits<double>::quiet_NaN();
};

}  // namespace core
}  // namespace cfdx