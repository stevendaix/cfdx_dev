#pragma once

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include <string>

namespace cfdx::core::linalg {

// ============================================
// Vertical Slice 1 — Hypre BoomerAMG Wrapper
// ============================================

//! AMG memory policy per v4 spec (A8-T08 / B2-T05)
enum class AMGMemoryPolicy {
    Low,     // Aggressive coarsening, minimized setup memory
    Balanced,  // Standard (HMIS coarsening)
    Fast    // Maximize convergence (strong threshold 0.5)
};

//! Hypre BoomerAMG wrapper implementing the preconditioner interface.
//!
//! Uses HYPRE's BoomerAMG library with configuration optimized for
//! CFD polyhedral meshes (coarsening HMIS, relaxation hybrid, etc.).
//! Memory policies control coarsening aggressiveness per v4 spec.
class HypreAMG {
public:
    HypreAMG();
    ~HypreAMG();

    //! Configure AMG with specified memory/performance policy
    void configure(AMGMemoryPolicy policy);

    //! Initialize AMG with matrix (setup phase — memory intensive)
    bool setup(const SparseMatrix& A, int num_partitions = 1);

    //! Apply preconditioner (y = M⁻¹ x) — called per CG iteration
    bool apply(const Vector& x, Vector& y) const;

    //! Measure memory usage (setup + solve peaks) for budget tracking
    struct AMGMemoryUsage {
        size_t setup_peak_bytes;
        size_t solve_peak_bytes;
        size_t temporary_bytes;
        int num_levels;
        int iterations;
    };
    AMGMemoryUsage measureMemoryUsage() const;

    //! Get estimated speedup vs no preconditioner (benchmark)
    double estimateSpeedup() const;

    //! Access underlying Hypre solver (for advanced tuning)
    void* hypreSolver();

private:
    AMGMemoryPolicy current_policy_;
    void* hypre_solver_;  // HYPRE_Solver (opaque pointer)
    bool is_initialized_;
};

} // namespace cfdx::core::linalg
