#pragma once

#include "cfdx/core/linalg/preconditioner.h"

#include <cstddef>
#include <memory>
#include <string>

namespace cfdx::core {

enum class AMGMemoryPolicy {
    Low,
    Balanced,
    Fast
};

// Dependency-free AMG recipe using the parts of the BoomerAMG method that fit
// CFDX's native algebra: strength-based aggregation, Galerkin coarse operators,
// smoothing, restriction/prolongation and a coarse solve. This is not HYPRE and
// does not claim binary or numerical equivalence with HYPRE BoomerAMG.
class NativeBoomerAMGPreconditioner final : public Preconditioner {
public:
    NativeBoomerAMGPreconditioner();
    ~NativeBoomerAMGPreconditioner() override;

    NativeBoomerAMGPreconditioner(const NativeBoomerAMGPreconditioner&) = delete;
    NativeBoomerAMGPreconditioner& operator=(const NativeBoomerAMGPreconditioner&) = delete;
    NativeBoomerAMGPreconditioner(NativeBoomerAMGPreconditioner&&) noexcept;
    NativeBoomerAMGPreconditioner& operator=(NativeBoomerAMGPreconditioner&&) noexcept;

    void configure(AMGMemoryPolicy policy);
    bool setup(const SparseMatrix& matrix) override;
    bool apply(const Vector& residual, Vector& correction) const override;
    const char* name() const override { return "NativeBoomerStyleAMG"; }

    bool is_ready() const noexcept;
    AMGMemoryPolicy memory_policy() const noexcept;
    std::size_t coarse_size() const noexcept;
    const std::string& last_error() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace cfdx::core

namespace cfdx::core::linalg {
using AMGMemoryPolicy = cfdx::core::AMGMemoryPolicy;
using HypreAMG [[deprecated("Use NativeBoomerAMGPreconditioner; this is not HYPRE")]] =
    cfdx::core::NativeBoomerAMGPreconditioner;
}
