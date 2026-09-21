#include "hypre_amg.h"
#include <cstdio>

using namespace cfdx::core::linalg;

HypreAMG::HypreAMG() : current_policy_(AMGMemoryPolicy::Balanced),
                         hypre_solver_(nullptr), is_initialized_(false) {
    std::printf("[HypreAMG] Constructor — using HYPRE_BoomerAMG backend (stub)\n");
}

HypreAMG::~HypreAMG() {
    std::printf("[HypreAMG] Destructor\n");
}

void HypreAMG::configure(AMGMemoryPolicy policy) {
    current_policy_ = policy;
    std::printf("[HypreAMG] Policy set: %s\n",
                (policy == AMGMemoryPolicy::Low) ? "Low (aggressive coarsening)" :
                (policy == AMGMemoryPolicy::Balanced) ? "Balanced (HMIS)" :
                "Fast (strong threshold 0.5)");
}

bool HypreAMG::setup(const SparseMatrix& A, int num_partitions) {
    (void)A;
    (void)num_partitions;
    // HYPRE is optional and is not linked in this build. Never report a
    // successful setup for an inactive backend.
    is_initialized_ = false;
    return false;
}

bool HypreAMG::apply(const Vector& rhs, Vector& x) const {
    (void)rhs;
    (void)x;
    return false;
}

HypreAMG::AMGMemoryUsage HypreAMG::measureMemoryUsage() const {
    AMGMemoryUsage usage{};
    usage.setup_peak_bytes = 0;  // Would measure HYPRE memory
    usage.solve_peak_bytes = 0;
    usage.temporary_bytes = 0;
    usage.num_levels = 0;
    usage.iterations = 0
    return usage;
}

double HypreAMG::estimateSpeedup() const {
    return 0.0;  // No benchmark is available until HYPRE is actually configured
}

void* HypreAMG::hypreSolver() {
    return hypre_solver_;
}
