#include "cfdx/core/linalg/hypre_amg.h"

#include "cfdx/core/linalg/amg_preconditioner.h"
#include "cfdx/core/linalg/linear_operator.h"

#include <cmath>
#include <limits>
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
    void update(const SparseMatrix& matrix) { matrix_ = matrix; }

private:
    SparseMatrix matrix_;
};

} // namespace

struct NativeAMGPreconditioner::Impl {
    explicit Impl(NativeAMGMethod selected_method, bool use_constant_null_space,
                 bool require_smoother_certificate,
                 bool algebraic_line_smoother, double prolongator_filter,
                 bool directional_semicoarsening)
        : method(selected_method),
          constant_null_space(use_constant_null_space),
          require_smoother_certificate(require_smoother_certificate),
          algebraic_line_smoother(algebraic_line_smoother),
          directional_semicoarsening(directional_semicoarsening),
          prolongator_filter(prolongator_filter) {}

    NativeAMGMethod method;
    bool constant_null_space = false;
    bool require_smoother_certificate = true;
    bool algebraic_line_smoother = false;
    bool directional_semicoarsening = false;
    double prolongator_filter = 0.0;
    std::size_t cf_interpolation_limit = 0;
    AMGMemoryPolicy policy = AMGMemoryPolicy::Balanced;
    std::unique_ptr<OwnedSparseOperator> op;
    std::unique_ptr<MatrixFreeVcyclePreconditioner> amg;
    std::string error;
    std::size_t hierarchy_builds = 0;
    std::size_t numeric_updates = 0;
};

NativeAMGPreconditioner::NativeAMGPreconditioner(
    NativeAMGMethod method, bool constant_null_space,
    bool require_smoother_certificate, bool algebraic_line_smoother,
    double prolongator_filter, bool directional_semicoarsening)
    : impl_(std::make_unique<Impl>(method, constant_null_space,
                                   require_smoother_certificate,
                                   algebraic_line_smoother, prolongator_filter,
                                   directional_semicoarsening)) {}

NativeAMGPreconditioner::~NativeAMGPreconditioner() = default;
NativeAMGPreconditioner::NativeAMGPreconditioner(
    NativeAMGPreconditioner&&) noexcept = default;
NativeAMGPreconditioner& NativeAMGPreconditioner::operator=(
    NativeAMGPreconditioner&&) noexcept = default;

void NativeAMGPreconditioner::configure(AMGMemoryPolicy policy) {
    if (impl_->policy != policy) {
        impl_->amg.reset();
        impl_->op.reset();
    }
    impl_->policy = policy;
    impl_->error.clear();
}

void NativeAMGPreconditioner::configure_cf_interpolation_limit(std::size_t max_entries) {
    if (max_entries == 1)
        throw std::invalid_argument("AMG C/F row cap must be 0 or >=2");
    if (impl_->cf_interpolation_limit != max_entries) {
        impl_->amg.reset();
        impl_->op.reset();
    }
    impl_->cf_interpolation_limit = max_entries;
    impl_->error.clear();
}

std::size_t NativeAMGPreconditioner::cf_interpolation_limit() const noexcept {
    return impl_->cf_interpolation_limit;
}

bool NativeAMGPreconditioner::setup(const SparseMatrix& matrix) {
    // Stage the next hierarchy before committing. An invalid new mesh or
    // failed factorization must not destroy a previously usable solver.
    if (matrix.n_rows() == 0 || matrix.n_rows() != matrix.n_cols() ||
        !matrix.is_consistent()) {
        impl_->error = "AMG setup requires a non-empty, square, finalized CSR matrix";
        return false;
    }

    double omega = 0.7;
    // Balanced native AMG is the qualification/default path. Use enough
    // Jacobi smoothing to make the V-cycle robust at the 4096/16384 scaling
    // points; the Fast policy remains explicit below.
    std::size_t pre_sweeps = 6;
    std::size_t post_sweeps = 6;
    if (impl_->policy == AMGMemoryPolicy::Low) {
        omega = 0.65;
        pre_sweeps = 2;
        post_sweeps = 2;
    } else if (impl_->policy == AMGMemoryPolicy::Fast) {
        omega = 0.75;
        pre_sweeps = 6;
        post_sweeps = 6;
    }

    if (impl_->directional_semicoarsening) {
        // A true block-line-GS sweep is stronger than the old Jacobi line
        // update; more sweeps can cost more than they save.
        pre_sweeps = 2;
        post_sweeps = 2;
    }
    auto next_op = std::make_unique<OwnedSparseOperator>(matrix);
    auto next_amg = std::make_unique<MatrixFreeVcyclePreconditioner>(
        *next_op, omega, pre_sweeps, post_sweeps, 0.25,
        impl_->policy == AMGMemoryPolicy::Low ? 15 : 25,
        impl_->directional_semicoarsening
            ? AMGInterpolationPolicy::DirectionalSemiCoarsening
            : impl_->method == NativeAMGMethod::SmoothedAggregation
                ? AMGInterpolationPolicy::SmoothedAggregation
                : AMGInterpolationPolicy::DirectCF,
        impl_->constant_null_space,
        impl_->require_smoother_certificate,
        impl_->directional_semicoarsening
            ? AMGSmootherPolicy::AlgebraicLineGaussSeidel
            : impl_->algebraic_line_smoother ? AMGSmootherPolicy::AlgebraicLines
                                            : AMGSmootherPolicy::DampedJacobi,
        impl_->prolongator_filter, impl_->cf_interpolation_limit);
    if (!next_amg->setup(next_op->matrix())) {
        impl_->error = "native AMG hierarchy construction failed";
        const auto& detail = next_amg->last_setup_error();
        if (!detail.empty()) impl_->error += ": " + detail;
        return false;
    }

    // Ensure the previous AMG (whose operator refers to the previous op)
    // is destroyed before its previous operator is released.
    impl_->amg = std::move(next_amg);
    impl_->op = std::move(next_op);
    ++impl_->hierarchy_builds;
    impl_->error.clear();
    return true;
}

bool NativeAMGPreconditioner::update_values(const SparseMatrix& matrix) {
    // update_values() is a numeric-only operation. Symbolic ownership remains
    // with the caller: an uninitialized hierarchy or a changed CSR graph must
    // be rejected rather than silently invoking setup().
    if (!impl_->op || !impl_->amg || matrix.n_rows() == 0 ||
        matrix.n_rows() != matrix.n_cols() || !matrix.is_consistent()) {
        impl_->error = "native AMG numeric update requires an existing compatible hierarchy";
        return false;
    }

    const SparseMatrix& current = impl_->op->matrix();
    if (current.n_rows() != matrix.n_rows() ||
        current.n_cols() != matrix.n_cols() ||
        current.nnz() != matrix.nnz() ||
        !current.is_consistent()) {
        impl_->error = "native AMG numeric update requires an unchanged CSR pattern";
        return false;
    }
    const auto* current_row = current.row_offsets_data();
    const auto* current_col = current.columns_data();
    const auto* new_row = matrix.row_offsets_data();
    const auto* new_col = matrix.columns_data();
    if (!std::equal(current_row, current_row + current.n_rows() + 1, new_row) ||
        !std::equal(current_col, current_col + current.nnz(), new_col)) {
        impl_->error = "native AMG numeric update requires an unchanged CSR pattern";
        return false;
    }

    if (!impl_->amg->update_values(matrix)) {
        impl_->error = "native AMG numeric update failed";
        return false;
    }
    impl_->op->update(matrix);
    ++impl_->numeric_updates;
    impl_->error.clear();
    return true;
}

bool NativeAMGPreconditioner::apply(
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

const char* NativeAMGPreconditioner::name() const {
    if (impl_->directional_semicoarsening)
        return "CFDXDirectionalSemiCoarseningLineGSAMG";
    if (impl_->algebraic_line_smoother)
        return impl_->method == NativeAMGMethod::SmoothedAggregation
            ? "CFDXAnisotropicSmoothedAggregationAMG"
            : "CFDXAnisotropicBoomerStyleAMG";
    return impl_->method == NativeAMGMethod::SmoothedAggregation
        ? "NativeSmoothedAggregationAMG"
        : "NativeBoomerStyleAMG";
}

bool NativeAMGPreconditioner::is_ready() const noexcept {
    return static_cast<bool>(impl_->amg);
}

NativeAMGMethod NativeAMGPreconditioner::method() const noexcept {
    return impl_->method;
}

bool NativeAMGPreconditioner::uses_constant_null_space() const noexcept {
    return impl_->constant_null_space;
}

AMGMemoryPolicy NativeAMGPreconditioner::memory_policy() const noexcept {
    return impl_->policy;
}

std::size_t NativeAMGPreconditioner::coarse_size() const noexcept {
    return impl_->amg ? impl_->amg->coarse_size() : 0;
}

std::vector<std::size_t> NativeAMGPreconditioner::hierarchy_level_sizes() const {
    return impl_->amg ? impl_->amg->hierarchy_level_sizes()
                      : std::vector<std::size_t>{};
}

std::vector<std::size_t> NativeAMGPreconditioner::hierarchy_level_nnz() const {
    return impl_->amg ? impl_->amg->hierarchy_level_nnz()
                      : std::vector<std::size_t>{};
}

double NativeAMGPreconditioner::operator_complexity() const {
    return impl_->amg ? impl_->amg->operator_complexity() : 0.0;
}

double NativeAMGPreconditioner::grid_complexity() const {
    return impl_->amg ? impl_->amg->grid_complexity() : 0.0;
}

std::size_t NativeAMGPreconditioner::finest_line_coverage() const {
    return impl_->amg ? impl_->amg->finest_line_coverage() : 0;
}

std::size_t NativeAMGPreconditioner::finest_directional_pairs() const {
    return impl_->amg ? impl_->amg->finest_directional_pairs() : 0;
}

double NativeAMGPreconditioner::finest_line_coupling_share() const {
    return impl_->amg ? impl_->amg->finest_line_coupling_share() : 0.0;
}

std::size_t NativeAMGPreconditioner::line_factor_fallback_levels() const {
    return impl_->amg ? impl_->amg->line_factor_fallback_levels() : 0;
}

std::size_t NativeAMGPreconditioner::numeric_line_fallbacks() const {
    return impl_->amg ? impl_->amg->numeric_line_fallbacks() : 0;
}

std::size_t NativeAMGPreconditioner::finest_longest_line() const {
    return impl_->amg ? impl_->amg->finest_longest_line() : 0;
}

auto NativeAMGPreconditioner::first_prolongation_nnz() const -> std::size_t {
    return impl_->amg ? impl_->amg->first_prolongation_nnz() : 0;
}

double NativeAMGPreconditioner::prolongation_row_sum_min() const {
    if (!impl_->amg) return std::numeric_limits<double>::quiet_NaN();
    return impl_->amg->prolongation_row_sum_min();
}

double NativeAMGPreconditioner::prolongation_row_sum_max() const {
    if (!impl_->amg) return std::numeric_limits<double>::quiet_NaN();
    return impl_->amg->prolongation_row_sum_max();
}

std::size_t NativeAMGPreconditioner::hierarchy_builds() const noexcept {
    return impl_->hierarchy_builds;
}

std::size_t NativeAMGPreconditioner::numeric_updates() const noexcept {
    return impl_->numeric_updates;
}

const std::string& NativeAMGPreconditioner::last_error() const noexcept {
    return impl_->error;
}

} // namespace cfdx::core
