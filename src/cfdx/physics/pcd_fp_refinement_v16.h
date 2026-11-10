#pragma once
// V16: optional, bounded inner-solve refinement. An already converged
// *auxiliary* Fp may fail the independently evaluated projected-original
// residual because the two evaluations have different floating-point errors.
// A tighter inner request is a genuine numerical correction, never a looser
// acceptance threshold. False certificates and invalid iterates fail closed.
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/linalg/solver_status.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {
struct PcdFpRefinementV16 {
    bool certified = false;
    std::size_t extra_solves = 0;
    double last_requested_tolerance = 0.0;
};

template <class Solve, class Verify>
inline PcdFpRefinementV16 refine_pcd_fp_v16(
    cfdx::core::Vector& iterate,
    cfdx::core::SolverResult& result,
    double original_tolerance,
    std::size_t maximum_extra_solves,
    Solve&& solve_at_tolerance,
    Verify&& verify_against_original)
{
    if (!(original_tolerance > 0.0) || !std::isfinite(original_tolerance) ||
        maximum_extra_solves > 3 || !iterate.is_valid())
        throw std::invalid_argument("PCD V16: invalid refinement contract");
    PcdFpRefinementV16 report;
    report.last_requested_tolerance = original_tolerance;
    if (verify_against_original(result, iterate)) {
        report.certified = true;
        return report;
    }
    for (std::size_t attempt = 0; attempt < maximum_extra_solves; ++attempt) {
        // This is strictly tighter than the original tolerance, and must be
        // respected by the actual nested Krylov solver, not just the monitor.
        const double tighter = original_tolerance * std::pow(0.1, attempt + 1.0);
        if (!(tighter > 0.0) || tighter >= original_tolerance ||
            !std::isfinite(tighter)) break;
        const auto previous = iterate;
        const auto previous_result = result;
        auto candidate_result = solve_at_tolerance(iterate, tighter);
        ++report.extra_solves;
        report.last_requested_tolerance = tighter;
        if (!iterate.is_valid()) {
            iterate = previous;
            result = previous_result;
            break;
        }
        if (verify_against_original(candidate_result, iterate)) {
            result = candidate_result;
            report.certified = true;
            return report;
        }
        // Leave an unsuccessful attempt available as the warm start for the
        // next, even if the inner solver returned MAX_ITER_REACHED. The final
        // caller MUST reject unless certified==true.
        result = candidate_result;
    }
    return report;
}
}  // namespace cfdx::physics
