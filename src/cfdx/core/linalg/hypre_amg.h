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

enum class NativeAMGMethod {
    BoomerStyle,
    SmoothedAggregation
};

// Shared dependency-free AMG lifecycle for CFDX's Boomer-style and smoothed-
// aggregation recipes. Neither method wraps or copies a third-party runtime.
class NativeAMGPreconditioner : public Preconditioner {
public:
    explicit NativeAMGPreconditioner(NativeAMGMethod method,
                                     bool constant_null_space = false);
    ~NativeAMGPreconditioner() override;

    NativeAMGPreconditioner(const NativeAMGPreconditioner&) = delete;
    NativeAMGPreconditioner& operator=(const NativeAMGPreconditioner&) = delete;
    NativeAMGPreconditioner(NativeAMGPreconditioner&&) noexcept;
    NativeAMGPreconditioner& operator=(NativeAMGPreconditioner&&) noexcept;

    void configure(AMGMemoryPolicy policy);
    bool setup(const SparseMatrix& matrix) override;
    bool update_values(const SparseMatrix& matrix) override;
    bool apply(const Vector& residual, Vector& correction) const override;
    const char* name() const override;

    bool is_ready() const noexcept;
    NativeAMGMethod method() const noexcept;
    bool uses_constant_null_space() const noexcept;
    AMGMemoryPolicy memory_policy() const noexcept;
    std::size_t coarse_size() const noexcept;
    std::size_t hierarchy_builds() const noexcept;
    std::size_t numeric_updates() const noexcept;
    const std::string& last_error() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class NativeBoomerAMGPreconditioner final : public NativeAMGPreconditioner {
public:
    explicit NativeBoomerAMGPreconditioner(bool constant_null_space = false)
        : NativeAMGPreconditioner(
              NativeAMGMethod::BoomerStyle, constant_null_space) {}
};

// Dependency-free smoothed-aggregation AMG following the transferable setup
// concepts of PETSc GAMG. This is not PETSc and does not copy PETSc source.
class NativeSmoothedAggregationAMGPreconditioner final
    : public NativeAMGPreconditioner {
public:
    explicit NativeSmoothedAggregationAMGPreconditioner(
        bool constant_null_space = false)
        : NativeAMGPreconditioner(
              NativeAMGMethod::SmoothedAggregation, constant_null_space) {}
};

} // namespace cfdx::core

namespace cfdx::core::linalg {
using AMGMemoryPolicy = cfdx::core::AMGMemoryPolicy;
using HypreAMG [[deprecated("Use NativeBoomerAMGPreconditioner; this is not HYPRE")]] =
    cfdx::core::NativeBoomerAMGPreconditioner;
}
