#include "cfdx/core/linalg/hypre_amg.h"

#include "cfdx/core/linalg/amg_preconditioner.h"
#include "cfdx/core/linalg/linear_operator.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace cfdx::core {
namespace {

class OwnedSparseOperator final : public LinearOperatorBase {
public:
    explicit OwnedSparseOperator(const SparseMatrix& matrix) : matrix_(matrix) {}

    std::size_t rows() const noexcept override { return matrix_.n_rows(); }
    std::size_t cols() const noexcept override { return matrix_.n_cols(); }

    void apply(const Vector& input, Vector& output) const override {
        if (input.size() != matrix_.n_cols()) {
            throw std::invalid_argument("NativeBoomerAMG input size mismatch");
        }
        const auto values = matrix_.matvec(input);
        if (output.size() != values.size()) output.resize(values.size());
        for (std::size_t i = 0; i < values.size(); ++i) output(i) = values[i];
    }

    const SparseMatrix& matrix() const noexcept { return matrix_; }

private:
    SparseMatrix matrix_;
};

} // namespace

struct NativeBoomerAMGPreconditioner::Impl {
    AMGMemoryPolicy policy = AMGMemoryPolicy::Balanced;
    std::unique_ptr<OwnedSparseOperator> op;
    std::unique_ptr<MatrixFreeVcyclePreconditioner> amg;
    std::string error;
};

NativeBoomerAMGPreconditioner::NativeBoomerAMGPreconditioner()
    : impl_(std::make_unique<Impl>()) {}

NativeBoomerAMGPreconditioner::~NativeBoomerAMGPreconditioner() = default;
NativeBoomerAMGPreconditioner::NativeBoomerAMGPreconditioner(
    NativeBoomerAMGPreconditioner&&) noexcept = default;
NativeBoomerAMGPreconditioner& NativeBoomerAMGPreconditioner::operator=(
    NativeBoomerAMGPreconditioner&&) noexcept = default;

void NativeBoomerAMGPreconditioner::configure(AMGMemoryPolicy policy) {
    if (impl_->policy != policy) {
        impl_->amg.reset();
        impl_->op.reset();
    }
    impl_->policy = policy;
    impl_->error.clear();
}

bool NativeBoomerAMGPreconditioner::setup(const SparseMatrix& matrix) {
    impl_->amg.reset();
    impl_->op.reset();
    if (matrix.n_rows() == 0 || matrix.n_rows() != matrix.n_cols() ||
        !matrix.is_consistent()) {
        impl_->error = "AMG setup requires a non-empty, square, finalized CSR matrix";
        return false;
    }

    double omega = 0.7;
    std::size_t pre_sweeps = 4;
    std::size_t post_sweeps = 4;
    if (impl_->policy == AMGMemoryPolicy::Low) {
        omega = 0.65;
        pre_sweeps = 2;
        post_sweeps = 2;
    } else if (impl_->policy == AMGMemoryPolicy::Fast) {
        omega = 0.75;
        pre_sweeps = 6;
        post_sweeps = 6;
    }

    auto next_op = std::make_unique<OwnedSparseOperator>(matrix);
    auto next_amg = std::make_unique<MatrixFreeVcyclePreconditioner>(
        *next_op, omega, pre_sweeps, post_sweeps, 0.25,
        impl_->policy == AMGMemoryPolicy::Low ? 15 : 25);
    if (!next_amg->setup(next_op->matrix())) {
        impl_->error = "native AMG hierarchy construction failed";
        return false;
    }

    impl_->amg.reset();
    impl_->op = std::move(next_op);
    impl_->amg = std::move(next_amg);
    impl_->error.clear();
    return true;
}

bool NativeBoomerAMGPreconditioner::apply(
    const Vector& residual, Vector& correction) const {
    if (!impl_->amg) {
        impl_->error = "native AMG preconditioner has not been set up";
        return false;
    }
    if (!impl_->amg->apply(residual, correction)) {
        impl_->error = "native AMG V-cycle failed";
        return false;
    }
    impl_->error.clear();
    return true;
}

bool NativeBoomerAMGPreconditioner::is_ready() const noexcept {
    return static_cast<bool>(impl_->amg);
}

AMGMemoryPolicy NativeBoomerAMGPreconditioner::memory_policy() const noexcept {
    return impl_->policy;
}

std::size_t NativeBoomerAMGPreconditioner::coarse_size() const noexcept {
    return impl_->amg ? impl_->amg->coarse_size() : 0;
}

const std::string& NativeBoomerAMGPreconditioner::last_error() const noexcept {
    return impl_->error;
}

} // namespace cfdx::core
