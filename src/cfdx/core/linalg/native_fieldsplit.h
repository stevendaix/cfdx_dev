#pragma once

#include "cfdx/core/linalg/block_preconditioner.h"

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::core {

enum class FieldSplitStrategy {
    Additive,
    DiagonalSchur,
    LowerSchur,
    UpperSchur,
    FullSchur
};

// Dependency-free two-field recipe matching the useful algebraic forms of a
// PETSc PCFIELDSPLIT setup. It reuses CFDX block and Schur preconditioners and
// intentionally does not expose PETSc names, options or runtime objects.
class NativeFieldSplitPreconditioner final : public Preconditioner {
public:
    NativeFieldSplitPreconditioner(std::vector<std::size_t> first_field,
                                   std::vector<std::size_t> second_field,
                                   FieldSplitStrategy strategy)
        : first_(std::move(first_field)),
          second_(std::move(second_field)),
          strategy_(strategy) {}

    bool setup(const SparseMatrix& matrix) override {
        implementation_.reset();
        std::unique_ptr<Preconditioner> next;
        if (strategy_ == FieldSplitStrategy::Additive) {
            std::vector<std::vector<std::size_t>> fields;
            fields.push_back(first_);
            fields.push_back(second_);
            next = std::make_unique<BlockDiagonalPreconditioner>(std::move(fields));
        } else {
            SchurFactorization factorization = SchurFactorization::Lower;
            if (strategy_ == FieldSplitStrategy::DiagonalSchur)
                factorization = SchurFactorization::Diagonal;
            else if (strategy_ == FieldSplitStrategy::UpperSchur)
                factorization = SchurFactorization::Upper;
            else if (strategy_ == FieldSplitStrategy::FullSchur)
                factorization = SchurFactorization::Full;
            next = std::make_unique<SchurComplementPreconditioner>(
                first_, second_, factorization);
        }
        if (!next->setup(matrix)) {
            last_error_ = "invalid field partition or singular field operator";
            return false;
        }
        implementation_ = std::move(next);
        last_error_.clear();
        return true;
    }

    bool apply(const Vector& residual, Vector& correction) const override {
        if (!implementation_) {
            last_error_ = "field split preconditioner has not been set up";
            return false;
        }
        if (!implementation_->apply(residual, correction)) {
            last_error_ = "field split application failed";
            return false;
        }
        last_error_.clear();
        return true;
    }

    const char* name() const override {
        switch (strategy_) {
            case FieldSplitStrategy::Additive: return "NativeFieldSplit-Additive";
            case FieldSplitStrategy::DiagonalSchur: return "NativeFieldSplit-DiagonalSchur";
            case FieldSplitStrategy::LowerSchur: return "NativeFieldSplit-LowerSchur";
            case FieldSplitStrategy::UpperSchur: return "NativeFieldSplit-UpperSchur";
            case FieldSplitStrategy::FullSchur: return "NativeFieldSplit-FullSchur";
        }
        return "NativeFieldSplit";
    }

    FieldSplitStrategy strategy() const noexcept { return strategy_; }
    const std::string& last_error() const noexcept { return last_error_; }

private:
    std::vector<std::size_t> first_;
    std::vector<std::size_t> second_;
    FieldSplitStrategy strategy_;
    std::unique_ptr<Preconditioner> implementation_;
    mutable std::string last_error_;
};

} // namespace cfdx::core
