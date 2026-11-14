#pragma once

// V38 read-only forensic capture of the three ACTUAL PCD pressure matrices
// after their auxiliary boundary/gauge treatment, before any nested solves.
// Each CFP01 file carries ZERO placeholder vectors. They are NOT physical
// right-hand sides or solution guesses; the actual Fp actions are separately
// captured by CFDX_PCD_FP_CAPTURE_PREFIX. This facility does not modify
// a coefficient, select a numerical method or accept a physical state.

#include "cfdx/physics/pcd_fp_capture.h"

#include <atomic>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

inline std::size_t capture_pcd_pressure_operators_v38(
    const std::string& prefix,
    const cfdx::core::SparseMatrix& mp,
    const cfdx::core::SparseMatrix& kp,
    const cfdx::core::SparseMatrix& fp) {
    const std::size_t n = mp.n_rows();
    if (prefix.empty() || n == 0 || n > 25000 ||
        mp.n_cols() != n || kp.n_rows() != n || kp.n_cols() != n ||
        fp.n_rows() != n || fp.n_cols() != n ||
        !mp.is_consistent() || !kp.is_consistent() || !fp.is_consistent())
        throw std::invalid_argument("V38 PCD operator capture: invalid prefix or pressure matrices");

    // Strict process-wide bound; even on a long nonlinear run, at most the
    // first three complete operator triplets are written. Do not consume a
    // slot before validating the arguments.
    static std::atomic<std::size_t> samples{0};
    const std::size_t slot = samples.fetch_add(1, std::memory_order_relaxed) + 1;
    if (slot > 3) return 0;
    const std::string base = prefix + "_state" + std::to_string(slot);
    const cfdx::core::Vector dummy(n, 0.0);
    write_pcd_fp_capture(base + "_mp.cfp", mp, dummy, dummy);
    write_pcd_fp_capture(base + "_kp.cfp", kp, dummy, dummy);
    write_pcd_fp_capture(base + "_fp.cfp", fp, dummy, dummy);
    return slot;
}

}  // namespace cfdx::physics
