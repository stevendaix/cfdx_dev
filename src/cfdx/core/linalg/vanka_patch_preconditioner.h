#pragma once

#include "cfdx/core/linalg/preconditioner.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::core {

// Algebraic, pressure-centered Vanka/overlapping Schwarz patches.
// Ordering: [Ux(0..N-1), Uy, Uz, p]. The pressure equation and every
// velocity row connected through D or G form a small local saddle block.
// Local inverses are pivoted; updates are multiplicative (latest residual),
// so the preconditioner must be paired with flexible nonsymmetric Krylov.
// No automatic insertion into the production COUPLED path is performed.
class VankaPressurePatchPreconditioner final : public Preconditioner {
public:
    explicit VankaPressurePatchPreconditioner(
        std::size_t n_cells, std::size_t max_patch = 24,
        double damping = 0.8, bool backward_sweep = true)
        : n_cells_(n_cells), max_patch_(max_patch), damping_(damping),
          backward_sweep_(backward_sweep) {}

    bool setup(const SparseMatrix& A) override {
        if (n_cells_ == 0 || A.n_rows() != 4*n_cells_ ||
            A.n_cols() != A.n_rows() || !A.is_consistent() ||
            max_patch_ < 2 || max_patch_ > 64 ||
            !(damping_ > 0.0 && damping_ <= 1.0) ||
            !std::isfinite(damping_)) {
            last_error_ = "invalid coupled CSR or Vanka parameters";
            return false;
        }
        const std::size_t nv = 3*n_cells_;
        const auto* ro = A.row_offsets_data();
        const auto* ci = A.columns_data();
        const auto* av = A.values_data();
        std::vector<std::vector<std::size_t>> pressure_to_u(n_cells_);
        // Both D and G are relevant on collocated meshes; the graph is
        // generally NOT the transpose of a compatible gradient graph.
        for (std::size_t i = 0; i < A.n_rows(); ++i)
            for (std::size_t k = ro[i]; k < ro[i+1]; ++k) {
                const std::size_t j = ci[k];
                if (i >= nv && j < nv && av[k] != 0.0)
                    pressure_to_u[i-nv].push_back(j);
                if (i < nv && j >= nv && av[k] != 0.0)
                    pressure_to_u[j-nv].push_back(i);
            }
        std::vector<Patch> next;
        next.reserve(n_cells_);
        std::vector<bool> covered(A.n_rows(), false);
        std::vector<double> uncovered_inverse_diag(A.n_rows(),0.0);
        std::size_t largest = 0;
        for (std::size_t p = 0; p < n_cells_; ++p) {
            Patch patch;
            patch.indices = std::move(pressure_to_u[p]);
            std::sort(patch.indices.begin(), patch.indices.end());
            patch.indices.erase(std::unique(patch.indices.begin(), patch.indices.end()),
                                patch.indices.end());
            patch.indices.push_back(nv+p);
            const std::size_t m = patch.indices.size();
            if (m > max_patch_) {
                last_error_ = "pressure Vanka patch exceeds configured limit";
                return false;
            }
            patch.lu.assign(m*m, 0.0);
            patch.pivots.assign(m, 0);
            for (std::size_t r = 0; r < m; ++r) {
                const std::size_t i = patch.indices[r];
                for (std::size_t k = ro[i]; k < ro[i+1]; ++k) {
                    auto it = std::lower_bound(patch.indices.begin(),
                                               patch.indices.end(),ci[k]);
                    if (it != patch.indices.end() && *it == ci[k]) {
                        const auto c = static_cast<std::size_t>(it-patch.indices.begin());
                        patch.lu[r*m+c] += av[k];
                    }
                }
            }
            if (!factor(patch)) {
                last_error_ = "Vanka patch is singular or ill-conditioned (pressure index " +
                              std::to_string(p) + ")";
                return false;
            }
            for (const auto row : patch.indices) covered[row] = true;
            largest = std::max(largest,m);
            next.push_back(std::move(patch));
        }
        // Vanka patches do not necessarily touch all momentum unknowns.
        // E.g. transverse velocity dofs on a 2D/1D flow have no gradient
        // link to pressure. Without a complement solve, the preconditioner
        // is rank deficient and FGMRES can catastrophically stagnate.
        for (std::size_t i=0;i<A.n_rows();++i) {
            if (covered[i]) continue;
            double diag=0.0,scale=0.0;
            for(std::size_t k=ro[i];k<ro[i+1];++k) {
                if(ci[k]==i) diag+=av[k];
                scale=std::max(scale,std::abs(av[k]));
            }
            if(!std::isfinite(diag) || !std::isfinite(scale) ||
               std::abs(diag)<=128.0*std::numeric_limits<double>::epsilon()*scale) {
                last_error_ = "uncovered velocity variable has no safe diagonal";
                return false;
            }
            uncovered_inverse_diag[i]=1.0/diag;
        }
        // Commit after all blocks have factored: failing setup() preserves
        // the old preconditioner and its original coefficient matrix.
        A_ = A;
        patches_ = std::move(next);
        uncovered_inverse_diag_ = std::move(uncovered_inverse_diag);
        longest_patch_ = largest;
        ready_ = true;
        last_error_.clear();
        return true;
    }

    bool update_values(const SparseMatrix& A) override {
        if (!ready_ || A.n_rows() != A_.n_rows() || A.nnz() != A_.nnz() ||
            !A.is_consistent()) {
            last_error_ = "Vanka numeric update requires existing CSR graph";
            return false;
        }
        for (std::size_t i=0;i<=A.n_rows();++i)
            if (A.row_offsets_data()[i] != A_.row_offsets_data()[i]) {
                last_error_ = "Vanka input row pattern changed";
                return false;
            }
        for (std::size_t k=0;k<A.nnz();++k)
            if (A.columns_data()[k] != A_.columns_data()[k]) {
                last_error_ = "Vanka input column pattern changed";
                return false;
            }
        // A graph-preserving coefficient change can still render an
        // individual pressure block singular: stage all pivots first.
        return setup(A);
    }

    bool apply(const Vector& rhs, Vector& correction) const override {
        if (!ready_ || rhs.size() != A_.n_rows() ||
            correction.size() != rhs.size() || !rhs.is_valid()) return false;
        const Vector* source = &rhs;
        Vector saved;
        if (&rhs == &correction) { saved = rhs; source = &saved; }
        correction.fill(0.0);
        const auto* ro = A_.row_offsets_data();
        const auto* ci = A_.columns_data();
        const auto* av = A_.values_data();
        std::vector<double> workspace(longest_patch_, 0.0);
        auto sweep = [&](bool reverse) {
            for (std::size_t step = 0; step < patches_.size(); ++step) {
                const Patch& patch = patches_[reverse ? patches_.size()-1-step : step];
                const std::size_t m=patch.indices.size();
                for(std::size_t local=0;local<m;++local) {
                    const std::size_t i=patch.indices[local];
                    double residual=(*source)(i);
                    for(std::size_t k=ro[i];k<ro[i+1];++k)
                        residual-=av[k]*correction(ci[k]);
                    workspace[local]=residual;
                }
                if (!solve(patch,workspace)) return false;
                for (std::size_t i=0;i<m;++i) {
                    correction(patch.indices[i]) += damping_ * workspace[i];
                    if (!std::isfinite(correction(patch.indices[i]))) return false;
                }
            }
            return true;
        };
        if (!sweep(false) || (backward_sweep_ && !sweep(true))) return false;
        for (std::size_t i=0;i<correction.size();++i) {
            if (uncovered_inverse_diag_[i] == 0.0) continue;
            double residual=(*source)(i);
            for (std::size_t k=ro[i];k<ro[i+1];++k)
                residual-=av[k]*correction(ci[k]);
            correction(i)+=uncovered_inverse_diag_[i]*residual;
            if(!std::isfinite(correction(i))) return false;
        }
        return correction.is_valid();
    }

    const char* name() const override { return "NativePressurePatchVanka"; }
    const std::string& last_error() const noexcept { return last_error_; }
    std::size_t patch_count() const noexcept { return patches_.size(); }
    std::size_t largest_patch() const noexcept { return longest_patch_; }

private:
    struct Patch {
        std::vector<std::size_t> indices;
        std::vector<std::size_t> pivots;
        std::vector<double> lu;
    };
    static bool factor(Patch& patch) {
        const std::size_t n=patch.indices.size();
        for(std::size_t k=0;k<n;++k) {
            std::size_t best=k;
            double best_value=std::abs(patch.lu[k*n+k]);
            for(std::size_t i=k+1;i<n;++i) {
                const double value=std::abs(patch.lu[i*n+k]);
                if (value>best_value) {best_value=value;best=i;}
            }
            if (!std::isfinite(best_value)) return false;
            // Relative to the surviving row norm, never an absolute
            // 1e-12 pivot that rejects physically small pressure units.
            double row_norm=0.0;
            for (std::size_t j=k;j<n;++j)
                row_norm=std::max(row_norm,std::abs(patch.lu[best*n+j]));
            if (!(best_value>128.0*std::numeric_limits<double>::epsilon()*row_norm))
                return false;
            patch.pivots[k]=best;
            if (best!=k) for(std::size_t j=0;j<n;++j)
                std::swap(patch.lu[k*n+j],patch.lu[best*n+j]);
            const double pivot=patch.lu[k*n+k];
            for(std::size_t i=k+1;i<n;++i) {
                patch.lu[i*n+k]/=pivot;
                for(std::size_t j=k+1;j<n;++j)
                    patch.lu[i*n+j]-=patch.lu[i*n+k]*patch.lu[k*n+j];
            }
        }
        return true;
    }
    static bool solve(const Patch& patch,std::vector<double>& rhs) {
        const std::size_t n=patch.indices.size();
        for(std::size_t k=0;k<n;++k)
            if(patch.pivots[k]!=k) std::swap(rhs[k],rhs[patch.pivots[k]]);
        for(std::size_t i=0;i<n;++i)
            for(std::size_t k=0;k<i;++k) rhs[i]-=patch.lu[i*n+k]*rhs[k];
        for(std::size_t ii=0;ii<n;++ii) {
            const std::size_t i=n-1-ii;
            for(std::size_t k=i+1;k<n;++k) rhs[i]-=patch.lu[i*n+k]*rhs[k];
            rhs[i]/=patch.lu[i*n+i];
            if(!std::isfinite(rhs[i])) return false;
        }
        return true;
    }
    std::size_t n_cells_;
    std::size_t max_patch_;
    double damping_;
    bool backward_sweep_;
    std::size_t longest_patch_=0;
    SparseMatrix A_;
    std::vector<Patch> patches_;
    std::vector<double> uncovered_inverse_diag_;
    bool ready_=false;
    std::string last_error_;
};

} // namespace cfdx::core