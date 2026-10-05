#pragma once

#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "cfdx/core/linalg/pcd_schur.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
#include "cfdx/core/linalg/simplerc_schur.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace cfdx::core {

// PETSc/HYPRE-inspired block reduction controls. PETSc PCFIELDSPLIT exposes
// Schur factorization variants (diag/lower/upper/full), while hypre MGR uses
// user-supplied block information and block-diagonal reductions. CFDX keeps
// the first implementation deliberately small and explicit: the field layout
// remains [Ux, Uy, Uz, p], while the reduction strategy is configurable.
enum class CoupledSchurFactorization {
    Diagonal,
    Lower,
    Upper,
    Full
};

enum class CoupledSchurVelocityApproximation {
    Block,
    Diagonal
};

enum class CoupledSchurApproximationModel {
    BlockLocal,
    PCD,
    LSC,
    BFBT,
    SIMPLE,
    SIMPLEC
};

struct CoupledBlockSchurOptions {
    CoupledSchurFactorization factorization =
        CoupledSchurFactorization::Full;
    CoupledSchurVelocityApproximation velocity_approximation =
        CoupledSchurVelocityApproximation::Block;
    CoupledSchurApproximationModel schur_approximation =
        CoupledSchurApproximationModel::BlockLocal;
    // PETSc flips the sign of the Schur solve for the diagonal saddle-point
    // form by default. CFDX does not assume that sign convention: it is an
    // explicit opt-in because CFDX pressure matrices may already be sign
    // normalized by the discretization.
    bool diagonal_schur_sign_flip = false;
};

// Block preconditioner for the pressure-based coupled ordering
// [Ux, Uy, Uz, p], using a cell-local 3x3 velocity block and an AMG solve
// of the explicitly assembled pressure Schur approximation
//
//   S~ = C - D M_b^{-1} G,
//
// where M_b is the block-diagonal collection of cell-local 3x3 velocity
// blocks -- NOT the inverse of the full velocity operator A_uu. The exact
// Schur complement S = C - D A_uu^{-1} G is not computed; that would
// require an inner solve for A_uu at every apply(). S~ is a valid and
// cheap Schur *approximation* in the SIMPLE/PCD family, and this class
// is intentionally named CoupledBlockSchurAMGPreconditioner to make the
// block-local approximation part of the type name (tracked in #481).
//
// The pressure AMG is a genuine preconditioner for S~, not a scalar Jacobi
// replacement. The class is intentionally serial/native until MPI/GPU block
// AMG is qualified separately.
class CoupledBlockSchurAMGPreconditioner final : public Preconditioner {
public:
    explicit CoupledBlockSchurAMGPreconditioner(
        std::size_t n_cells,
        CoupledBlockSchurOptions options = {})
        : n_cells_(n_cells), nv_(3 * n_cells),
          options_(options), pressure_amg_(false) {}

    void set_pcd_schur(std::unique_ptr<PcdSchurApproximation> pcd) {
        pcd_schur_ = std::move(pcd);
    }

    void set_algebraic_schur(std::unique_ptr<SchurApproximation> schur) {
        algebraic_schur_ = std::move(schur);
    }

    void set_simpler_schur(std::unique_ptr<SimplerSchurApproximation> schur) {
        simpler_schur_ = std::move(schur);
    }

    bool setup(const SparseMatrix& A) override {
        reset();
        if (n_cells_ == 0 || A.n_rows() != A.n_cols() ||
            A.n_rows() != 4 * n_cells_)
            return fail("expected a square 4N coupled matrix");

        if (!prepare_numeric_state(A)) return false;
        if (options_.schur_approximation == CoupledSchurApproximationModel::PCD) {
            if (!pcd_schur_)
                return fail("PCD Schur approximation was not configured");
            Auu_ = extract_block(A, 0, 0, 3 * n_cells_, 3 * n_cells_);
            G_ = extract_block(A, 0, 1, 3 * n_cells_, n_cells_);
            D_ = extract_block(A, 1, 0, n_cells_, 3 * n_cells_);
            C_ = extract_block(A, 1, 1, n_cells_, n_cells_);
            pcd_blocks_ =
                std::make_unique<BlockOperator>(Auu_, G_, D_, C_);
            if (!pcd_schur_->setup(*pcd_blocks_))
                return fail("PCD Schur approximation setup failed");
            pcd_ready_ = true;
            ready_ = true;
            return true;
        }
        if (options_.schur_approximation == CoupledSchurApproximationModel::SIMPLE ||
            options_.schur_approximation == CoupledSchurApproximationModel::SIMPLEC) {
            if (!simpler_schur_)
                return fail("SIMPLE/SIMPLEC Schur approximation was not configured");
            Auu_ = extract_block(A, 0, 0, 3 * n_cells_, 3 * n_cells_);
            G_ = extract_block(A, 0, 1, 3 * n_cells_, n_cells_);
            D_ = extract_block(A, 1, 0, n_cells_, 3 * n_cells_);
            C_ = extract_block(A, 1, 1, n_cells_, n_cells_);
            simpler_blocks_ = std::make_unique<BlockOperator>(Auu_, G_, D_, C_);
            if (!simpler_schur_->setup(*simpler_blocks_))
                return fail("SIMPLE/SIMPLEC Schur approximation setup failed");
            schur_ = simpler_schur_->assembled_operator();
            if (schur_.n_rows() != n_cells_ || schur_.nnz() == 0)
                return fail("SIMPLE/SIMPLEC Schur assembly failed");
            if (!pressure_amg_.setup(schur_))
                return fail(std::string("pressure SIMPLE/SIMPLEC AMG setup failed: ") +
                            pressure_amg_.last_error());
            simpler_ready_ = true;
            ready_ = true;
            return true;
        }
        if (options_.schur_approximation == CoupledSchurApproximationModel::LSC ||
            options_.schur_approximation == CoupledSchurApproximationModel::BFBT) {
            if (!algebraic_schur_)
                return fail("LSC/BFBt Schur approximation was not configured");
            Auu_ = extract_block(A, 0, 0, 3 * n_cells_, 3 * n_cells_);
            G_ = extract_block(A, 0, 1, 3 * n_cells_, n_cells_);
            D_ = extract_block(A, 1, 0, n_cells_, 3 * n_cells_);
            C_ = extract_block(A, 1, 1, n_cells_, n_cells_);
            algebraic_blocks_ = std::make_unique<BlockOperator>(Auu_, G_, D_, C_);
            if (!algebraic_schur_->setup(*algebraic_blocks_))
                return fail("LSC/BFBt Schur approximation setup failed");
            algebraic_ready_ = true;
            ready_ = true;
            return true;
        }
        if (!pressure_amg_.setup(schur_)) {
            return fail(std::string("pressure Schur AMG setup failed: ") +
                        pressure_amg_.last_error());
        }
        ready_ = true;
        return true;
    }

    bool update_values(const SparseMatrix& A) override {
        // Rebuild only the numeric Schur data and velocity inverses here.
        // In particular, do not call setup() on a temporary AMG object:
        // values-only updates must reuse the existing symbolic hierarchy.
        // Numeric updates never own the symbolic rebuild. The caller must
        // call setup() explicitly for an uninitialized hierarchy or a new
        // matrix size.
        if (!ready_ || A.n_rows() != A.n_cols() ||
            A.n_rows() != 4 * n_cells_ || !A.is_consistent())
            return fail_update("Schur numeric update requires an existing compatible hierarchy");

        // Reject a changed input CSR graph before rebuilding any temporary
        // velocity inverses or Schur matrix. The candidate Schur comparison
        // below remains necessary because numerical cancellation can also alter
        // the resulting Schur graph even when the input graph is unchanged.
        if (!same_input_pattern(A))
            return fail_update("coupled CSR pattern changed; explicit setup() required");

        CoupledBlockSchurAMGPreconditioner candidate(n_cells_, options_);
        if (!candidate.prepare_numeric_state(A)) {
            last_error_ = candidate.last_error_;
            return false;
        }

        if (options_.schur_approximation == CoupledSchurApproximationModel::SIMPLE ||
            options_.schur_approximation == CoupledSchurApproximationModel::SIMPLEC) {
            if (!simpler_ready_ || !simpler_schur_ || !simpler_blocks_)
                return fail_update("SIMPLE/SIMPLEC Schur approximation is not initialized");
            const SparseMatrix new_Auu = extract_block(A, 0, 0, 3 * n_cells_, 3 * n_cells_);
            const SparseMatrix new_G = extract_block(A, 0, 1, 3 * n_cells_, n_cells_);
            const SparseMatrix new_D = extract_block(A, 1, 0, n_cells_, 3 * n_cells_);
            const SparseMatrix new_C = extract_block(A, 1, 1, n_cells_, n_cells_);
            if (!same_pattern(Auu_, new_Auu) || !same_pattern(G_, new_G) ||
                !same_pattern(D_, new_D) || !same_pattern(C_, new_C))
                return fail_update("SIMPLE/SIMPLEC coupled block graph changed; explicit setup() required");
            std::copy(new_Auu.values_data(), new_Auu.values_data() + new_Auu.nnz(), Auu_.values_data());
            std::copy(new_G.values_data(), new_G.values_data() + new_G.nnz(), G_.values_data());
            std::copy(new_D.values_data(), new_D.values_data() + new_D.nnz(), D_.values_data());
            std::copy(new_C.values_data(), new_C.values_data() + new_C.nnz(), C_.values_data());
            simpler_blocks_ = std::make_unique<BlockOperator>(Auu_, G_, D_, C_);
            if (!simpler_schur_->update_values(*simpler_blocks_))
                return fail_update("SIMPLE/SIMPLEC Schur numeric update failed");
            const SparseMatrix new_schur = simpler_schur_->assembled_operator();
            if (new_schur.n_rows() != n_cells_ || new_schur.nnz() == 0)
                return fail_update("SIMPLE/SIMPLEC Schur assembly failed during numeric update");
            if (!same_pattern(schur_, new_schur))
                return fail_update("SIMPLE/SIMPLEC Schur CSR pattern changed; explicit setup() required");
            if (!pressure_amg_.update_values(new_schur)) {
                last_error_ = pressure_amg_.last_error();
                return false;
            }
            schur_ = new_schur;
            velocity_inv_ = std::move(candidate.velocity_inv_);
            velocity_inv_diag_ = std::move(candidate.velocity_inv_diag_);
            row_ = std::move(candidate.row_);
            col_ = std::move(candidate.col_);
            val_ = std::move(candidate.val_);
            return true;
        }

        if (!same_pattern(schur_, candidate.schur_)) {
            // A changed Schur graph invalidates the symbolic hierarchy.
            return fail_update("Schur CSR pattern changed; explicit setup() required");
        }

        if (options_.schur_approximation == CoupledSchurApproximationModel::PCD) {
            if (!pcd_ready_ || !pcd_schur_)
                return fail_update("PCD Schur approximation is not initialized");
            const SparseMatrix new_Auu =
                extract_block(A, 0, 0, 3 * n_cells_, 3 * n_cells_);
            const SparseMatrix new_G =
                extract_block(A, 0, 1, 3 * n_cells_, n_cells_);
            const SparseMatrix new_D =
                extract_block(A, 1, 0, n_cells_, 3 * n_cells_);
            const SparseMatrix new_C =
                extract_block(A, 1, 1, n_cells_, n_cells_);
            if (!same_pattern(Auu_, new_Auu) ||
                !same_pattern(G_, new_G) ||
                !same_pattern(D_, new_D) ||
                !same_pattern(C_, new_C))
                return fail_update("PCD coupled block graph changed; explicit setup() required");

            std::copy(new_Auu.values_data(),
                      new_Auu.values_data() + new_Auu.nnz(), Auu_.values_data());
            std::copy(new_G.values_data(),
                      new_G.values_data() + new_G.nnz(), G_.values_data());
            std::copy(new_D.values_data(),
                      new_D.values_data() + new_D.nnz(), D_.values_data());
            std::copy(new_C.values_data(),
                      new_C.values_data() + new_C.nnz(), C_.values_data());

            if (!pcd_blocks_ || !pcd_schur_->update_values(*pcd_blocks_))
                return fail_update("PCD Schur numeric update failed");

            velocity_inv_ = std::move(candidate.velocity_inv_);
            velocity_inv_diag_ = std::move(candidate.velocity_inv_diag_);
            row_ = std::move(candidate.row_);
            col_ = std::move(candidate.col_);
            val_ = std::move(candidate.val_);
            schur_ = std::move(candidate.schur_);
            return true;
        }

        if (options_.schur_approximation == CoupledSchurApproximationModel::LSC ||
            options_.schur_approximation == CoupledSchurApproximationModel::BFBT) {
            if (!algebraic_ready_ || !algebraic_schur_ || !algebraic_blocks_)
                return fail_update("LSC/BFBt Schur approximation is not initialized");
            const SparseMatrix new_Auu = extract_block(A, 0, 0, 3 * n_cells_, 3 * n_cells_);
            const SparseMatrix new_G = extract_block(A, 0, 1, 3 * n_cells_, n_cells_);
            const SparseMatrix new_D = extract_block(A, 1, 0, n_cells_, 3 * n_cells_);
            const SparseMatrix new_C = extract_block(A, 1, 1, n_cells_, n_cells_);
            if (!same_pattern(Auu_, new_Auu) || !same_pattern(G_, new_G) ||
                !same_pattern(D_, new_D) || !same_pattern(C_, new_C))
                return fail_update("LSC/BFBt coupled block graph changed; explicit setup() required");
            std::copy(new_Auu.values_data(), new_Auu.values_data() + new_Auu.nnz(), Auu_.values_data());
            std::copy(new_G.values_data(), new_G.values_data() + new_G.nnz(), G_.values_data());
            std::copy(new_D.values_data(), new_D.values_data() + new_D.nnz(), D_.values_data());
            std::copy(new_C.values_data(), new_C.values_data() + new_C.nnz(), C_.values_data());
            algebraic_blocks_ = std::make_unique<BlockOperator>(Auu_, G_, D_, C_);
            if (!algebraic_schur_->update_values(*algebraic_blocks_))
                return fail_update("LSC/BFBt Schur numeric update failed");
            velocity_inv_ = std::move(candidate.velocity_inv_);
            velocity_inv_diag_ = std::move(candidate.velocity_inv_diag_);
            row_ = std::move(candidate.row_);
            col_ = std::move(candidate.col_);
            val_ = std::move(candidate.val_);
            schur_ = std::move(candidate.schur_);
            return true;
        }

        if (!pressure_amg_.update_values(candidate.schur_)) {
            last_error_ = pressure_amg_.last_error();
            return false;
        }

        velocity_inv_ = std::move(candidate.velocity_inv_);
        velocity_inv_diag_ = std::move(candidate.velocity_inv_diag_);
        row_ = std::move(candidate.row_);
        col_ = std::move(candidate.col_);
        val_ = std::move(candidate.val_);
        schur_ = std::move(candidate.schur_);
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (!ready_ || r.size() != 4 * n_cells_ || z.size() != r.size())
            return false;

        Vector ru(nv_), rp(n_cells_), y(nv_), schur_rhs(n_cells_);
        for (std::size_t i = 0; i < nv_; ++i) ru(i) = r(i);
        for (std::size_t c = 0; c < n_cells_; ++c) rp(c) = r(nv_ + c);

        const auto solve_velocity = [&](const Vector& rhs_u, const Vector& pressure,
                                         Vector& out_u) {
            for (std::size_t c = 0; c < n_cells_; ++c) {
                double in[3] = {
                    rhs_u(c), rhs_u(n_cells_ + c), rhs_u(2 * n_cells_ + c)};
                for (std::size_t rcomp = 0; rcomp < 3; ++rcomp) {
                    const std::size_t gr = rcomp * n_cells_ + c;
                    for (std::size_t k = row_[gr]; k < row_[gr + 1]; ++k) {
                        if (col_[k] >= nv_)
                            in[rcomp] -= val_[k] * pressure(col_[k] - nv_);
                    }
                }
                for (std::size_t rcomp = 0; rcomp < 3; ++rcomp) {
                    double out = 0.0;
                    for (std::size_t q = 0; q < 3; ++q)
                        out += velocity_inv_[c][3 * rcomp + q] * in[q];
                    out_u(rcomp * n_cells_ + c) = out;
                }
            }
        };

        const auto solve_pressure = [&](const Vector& pressure_rhs, Vector& pressure) {
            if (options_.schur_approximation == CoupledSchurApproximationModel::PCD) {
                if (!pcd_ready_ || !pcd_schur_ ||
                    !pcd_schur_->apply(pressure_rhs, pressure))
                    return false;
            } else if (options_.schur_approximation == CoupledSchurApproximationModel::LSC ||
                       options_.schur_approximation == CoupledSchurApproximationModel::BFBT) {
                if (!algebraic_ready_ || !algebraic_schur_ ||
                    !algebraic_schur_->apply(pressure_rhs, pressure))
                    return false;
            } else if (!pressure_amg_.apply(pressure_rhs, pressure)) return false;
            if (options_.diagonal_schur_sign_flip) {
                for (std::size_t c = 0; c < n_cells_; ++c)
                    pressure(c) = -pressure(c);
            }
            return true;
        };

        Vector zp(n_cells_, 0.0);

        // PCFIELDSPLIT-style factorization choices:
        //   diagonal: M^-1 and S^-1 independently
        //   lower:    M^-1 followed by the Schur correction
        //   upper:    Schur solve followed by velocity back-substitution
        //   full:     lower + upper triangular factors (the default).
        if (options_.factorization == CoupledSchurFactorization::Diagonal) {
            apply_velocity_inverse(ru, y);
            Vector pressure_rhs(n_cells_);
            for (std::size_t c = 0; c < n_cells_; ++c)
                pressure_rhs(c) = rp(c);
            if (!solve_pressure(pressure_rhs, zp)) return false;
            for (std::size_t i = 0; i < nv_; ++i) z(i) = y(i);
        } else if (options_.factorization == CoupledSchurFactorization::Lower) {
            apply_velocity_inverse(ru, y);
            Vector pressure_rhs(n_cells_);
            for (std::size_t c = 0; c < n_cells_; ++c) {
                double value = rp(c);
                const std::size_t pr = nv_ + c;
                for (std::size_t k = row_[pr]; k < row_[pr + 1]; ++k)
                    if (col_[k] < nv_) value -= val_[k] * y(col_[k]);
                pressure_rhs(c) = value;
            }
            if (!solve_pressure(pressure_rhs, zp)) return false;
            for (std::size_t i = 0; i < nv_; ++i) z(i) = y(i);
        } else if (options_.factorization == CoupledSchurFactorization::Upper) {
            Vector pressure_rhs(n_cells_);
            for (std::size_t c = 0; c < n_cells_; ++c)
                pressure_rhs(c) = rp(c);
            if (!solve_pressure(pressure_rhs, zp)) return false;
            solve_velocity(ru, zp, y);
            for (std::size_t i = 0; i < nv_; ++i) z(i) = y(i);
        } else {
            apply_velocity_inverse(ru, y);
            Vector pressure_rhs(n_cells_);
            for (std::size_t c = 0; c < n_cells_; ++c) {
                double value = rp(c);
                const std::size_t pr = nv_ + c;
                for (std::size_t k = row_[pr]; k < row_[pr + 1]; ++k)
                    if (col_[k] < nv_) value -= val_[k] * y(col_[k]);
                pressure_rhs(c) = value;
            }
            if (!solve_pressure(pressure_rhs, zp)) return false;
            // Full block-LDU inverse: after the Schur solve, apply the
            // upper triangular correction u <- M^-1(r_u - G p).
            solve_velocity(ru, zp, y);
            for (std::size_t i = 0; i < nv_; ++i) z(i) = y(i);
        }
        for (std::size_t c = 0; c < n_cells_; ++c) z(nv_ + c) = zp(c);
        return z.is_valid();
    }

    const char* name() const override {
        return "CoupledBlockSchur-AMG";
    }

    bool is_ready() const noexcept { return ready_; }
    const SparseMatrix& schur_matrix() const noexcept { return schur_; }
    std::size_t pressure_coarse_size() const noexcept {
        return pressure_amg_.coarse_size();
    }
    std::size_t pressure_hierarchy_builds() const noexcept {
        return pressure_amg_.hierarchy_builds();
    }
    std::size_t pressure_numeric_updates() const noexcept {
        return pressure_amg_.numeric_updates();
    }
    CoupledSchurFactorization factorization() const noexcept {
        return options_.factorization;
    }
    CoupledSchurVelocityApproximation velocity_approximation() const noexcept {
        return options_.velocity_approximation;
    }
    const std::string& last_error() const noexcept { return last_error_; }

private:
    static SparseMatrix extract_block(const SparseMatrix& A,
                                      std::size_t row_block,
                                      std::size_t col_block,
                                      std::size_t row_size,
                                      std::size_t col_size) {
        SparseMatrix block(row_size, col_size);
        // For the 4N ordering [Ux,Uy,Uz,p], velocity occupies 3N rows/cols.
        const std::size_t velocity_size = A.n_rows() * 3 / 4;
        const std::size_t roffset = row_block == 0 ? 0 : velocity_size;
        const std::size_t coffset = col_block == 0 ? 0 : velocity_size;
        for (std::size_t r = 0; r < row_size; ++r) {
            const std::size_t gr = roffset + r;
            for (std::uint32_t k = A.row_offsets_data()[gr];
                 k < A.row_offsets_data()[gr + 1]; ++k) {
                const std::size_t gc = A.columns_data()[k];
                if (gc < coffset || gc >= coffset + col_size) continue;
                block.push_back(r, gc - coffset, A.values_data()[k]);
            }
        }
        block.finalize();
        return block;
    }

    static bool invert3x3(const std::array<double, 9>& a,
                          std::array<double, 9>& inv) {
        const double det =
            a[0] * (a[4] * a[8] - a[5] * a[7]) -
            a[1] * (a[3] * a[8] - a[5] * a[6]) +
            a[2] * (a[3] * a[7] - a[4] * a[6]);
        double scale = 0.0;
        for (double x : a) scale = std::max(scale, std::abs(x));
        if (!std::isfinite(det) || scale == 0.0 ||
            std::abs(det) <= 128.0 * std::numeric_limits<double>::epsilon() *
                                 scale * scale * scale)
            return false;
        inv[0] = (a[4]*a[8]-a[5]*a[7])/det;
        inv[1] = (a[2]*a[7]-a[1]*a[8])/det;
        inv[2] = (a[1]*a[5]-a[2]*a[4])/det;
        inv[3] = (a[5]*a[6]-a[3]*a[8])/det;
        inv[4] = (a[0]*a[8]-a[2]*a[6])/det;
        inv[5] = (a[2]*a[3]-a[0]*a[5])/det;
        inv[6] = (a[3]*a[7]-a[4]*a[6])/det;
        inv[7] = (a[1]*a[6]-a[0]*a[7])/det;
        inv[8] = (a[0]*a[4]-a[1]*a[3])/det;
        for (double x : inv) if (!std::isfinite(x)) return false;
        return true;
    }

    void apply_velocity_inverse(const Vector& r, Vector& z) const {
        for (std::size_t c = 0; c < n_cells_; ++c) {
            const double in[3] = {r(c), r(n_cells_ + c), r(2*n_cells_ + c)};
            for (std::size_t i = 0; i < 3; ++i) {
                double out = 0.0;
                for (std::size_t j = 0; j < 3; ++j)
                    out += velocity_inv_[c][3*i+j] * in[j];
                z(i*n_cells_ + c) = out;
            }
        }
    }

    bool prepare_numeric_state(const SparseMatrix& A) {
        if (A.n_rows() != A.n_cols() || A.n_rows() != 4 * n_cells_)
            return fail("expected a square 4N coupled matrix");

        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        row_.assign(row, row + A.n_rows() + 1);
        col_.assign(col, col + A.nnz());
        val_.assign(val, val + A.nnz());

        velocity_inv_.assign(n_cells_, {});
        velocity_inv_diag_.assign(n_cells_, std::array<double, 3>{});
        for (std::size_t c = 0; c < n_cells_; ++c) {
            std::array<double, 9> block{};
            for (std::size_t r = 0; r < 3; ++r) {
                const std::size_t gr = r * n_cells_ + c;
                for (std::size_t k = row[gr]; k < row[gr + 1]; ++k) {
                    const std::size_t gc = col[k];
                    for (std::size_t q = 0; q < 3; ++q) {
                        if (gc == q * n_cells_ + c) {
                            block[3 * r + q] += val[k];
                            break;
                        }
                    }
                }
            }
            if (!invert3x3(block, velocity_inv_[c]))
                return fail("singular or invalid cell-local velocity block");
            for (std::size_t q = 0; q < 3; ++q) {
                const double d = block[3 * q + q];
                if (!std::isfinite(d) || d == 0.0)
                    return fail("singular or invalid velocity diagonal block");
                velocity_inv_diag_[c][q] = 1.0 / d;
            }
        }
        return build_schur(A);
    }

    bool build_schur(const SparseMatrix& A) {
        // Build a sparse representation of G grouped by velocity cell:
        // G_{(component,cell),q}.  This keeps the Schur assembly proportional
        // to the actual coupling graph instead of scanning every pressure
        // column for every D entry.
        struct GEntry {
            std::size_t pressure;
            std::array<double, 3> g{};
        };
        std::vector<std::vector<GEntry>> g_by_cell(n_cells_);
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();

        for (std::size_t cell = 0; cell < n_cells_; ++cell) {
            for (std::size_t component = 0; component < 3; ++component) {
                const std::size_t gr = component * n_cells_ + cell;
                for (std::size_t k = row[gr]; k < row[gr + 1]; ++k) {
                    if (col[k] < nv_) continue;
                    const std::size_t q = col[k] - nv_;
                    auto it = std::find_if(
                        g_by_cell[cell].begin(), g_by_cell[cell].end(),
                        [q](const GEntry& e) { return e.pressure == q; });
                    if (it == g_by_cell[cell].end()) {
                        GEntry e;
                        e.pressure = q;
                        e.g[component] = val[k];
                        g_by_cell[cell].push_back(e);
                    } else {
                        it->g[component] += val[k];
                    }
                }
            }
        }

        std::vector<std::vector<std::pair<std::size_t, double>>> rows(n_cells_);
        for (std::size_t c = 0; c < n_cells_; ++c) {
            const std::size_t pr = nv_ + c;

            auto add = [&](std::size_t pressure, double value) {
                for (auto& e : rows[c]) {
                    if (e.first == pressure) {
                        e.second += value;
                        return;
                    }
                }
                rows[c].push_back({pressure, value});
            };

            // C block.
            for (std::size_t k = row[pr]; k < row[pr + 1]; ++k) {
                if (col[k] >= nv_) {
                    add(col[k] - nv_, val[k]);
                }
            }

            // -D M_b^{-1} G.
            for (std::size_t k = row[pr]; k < row[pr + 1]; ++k) {
                const std::size_t u = col[k];
                if (u >= nv_) continue;

                const std::size_t cell = u % n_cells_;
                const std::size_t component = u / n_cells_;
                const double d = val[k];

                for (const auto& g : g_by_cell[cell]) {
                    double m_inv_g_component = 0.0;
                    for (std::size_t j = 0; j < 3; ++j) {
                        m_inv_g_component +=
                            options_.velocity_approximation ==
                                    CoupledSchurVelocityApproximation::Block
                                ? velocity_inv_[cell][3 * component + j] * g.g[j]
                                : (component == j
                                       ? velocity_inv_diag_[cell][component] * g.g[j]
                                       : 0.0);
                    }
                    add(g.pressure, -d * m_inv_g_component);
                }
            }
        }

        return finalize_schur(rows);
    }

    bool finalize_schur(const std::vector<std::vector<std::pair<std::size_t,double>>>& rows) {
        schur_ = SparseMatrix(n_cells_, n_cells_);
        for (std::size_t r = 0; r < n_cells_; ++r) {
            auto entries = rows[r];
            std::sort(entries.begin(), entries.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });
            for (const auto& [c,v] : entries) {
                if (!std::isfinite(v)) return fail("non-finite Schur coefficient");
                schur_.push_back(r,c,v);
            }
        }
        schur_.finalize();
        return schur_.n_rows() == n_cells_ && schur_.nnz() > 0;
    }

    bool same_input_pattern(const SparseMatrix& A) const {
        if (A.n_rows() != row_.size() - 1 ||
            A.nnz() != col_.size())
            return false;
        const auto* input_row = A.row_offsets_data();
        const auto* input_col = A.columns_data();
        for (std::size_t i = 0; i < row_.size(); ++i)
            if (input_row[i] != row_[i]) return false;
        for (std::size_t i = 0; i < col_.size(); ++i)
            if (input_col[i] != col_[i]) return false;
        return true;
    }

    static bool same_pattern(const SparseMatrix& a, const SparseMatrix& b) {
        if (a.n_rows()!=b.n_rows() || a.n_cols()!=b.n_cols() || a.nnz()!=b.nnz())
            return false;
        const auto* ar=a.row_offsets_data(); const auto* br=b.row_offsets_data();
        for (std::size_t i=0;i<=a.n_rows();++i) if (ar[i]!=br[i]) return false;
        const auto* ac=a.columns_data(); const auto* bc=b.columns_data();
        for (std::size_t i=0;i<a.nnz();++i) if (ac[i]!=bc[i]) return false;
        return true;
    }

    bool fail(const std::string& message) {
        last_error_=message; ready_=false; return false;
    }

    bool fail_update(const std::string& message) {
        // A rejected numeric update must not invalidate the still-valid old
        // hierarchy. The caller can continue using it or explicitly call setup()
        // when a symbolic rebuild is intended.
        last_error_=message;
        return false;
    }

    void reset() {
        ready_=false; pcd_ready_=false; pcd_blocks_.reset(); algebraic_ready_=false; algebraic_blocks_.reset(); simpler_ready_=false; simpler_blocks_.reset(); last_error_.clear(); velocity_inv_.clear();
        velocity_inv_diag_.clear();
        row_.clear(); col_.clear(); val_.clear(); schur_=SparseMatrix();
    }

    std::size_t n_cells_{0}, nv_{0};
    CoupledBlockSchurOptions options_{};
    std::vector<std::array<double,9>> velocity_inv_;
    std::vector<std::array<double,3>> velocity_inv_diag_;
    std::vector<std::size_t> row_;
    std::vector<std::uint32_t> col_;
    std::vector<double> val_;
    SparseMatrix schur_;
    SparseMatrix Auu_;
    SparseMatrix G_;
    SparseMatrix D_;
    SparseMatrix C_;
    NativeBoomerAMGPreconditioner pressure_amg_;
    std::unique_ptr<PcdSchurApproximation> pcd_schur_;
    std::unique_ptr<BlockOperator> pcd_blocks_;
    std::unique_ptr<SchurApproximation> algebraic_schur_;
    std::unique_ptr<BlockOperator> algebraic_blocks_;
    std::unique_ptr<SimplerSchurApproximation> simpler_schur_;
    std::unique_ptr<BlockOperator> simpler_blocks_;
    bool pcd_ready_{false};
    bool algebraic_ready_{false};
    bool simpler_ready_{false};
    bool ready_{false};
    std::string last_error_;
};

} // namespace cfdx::core
