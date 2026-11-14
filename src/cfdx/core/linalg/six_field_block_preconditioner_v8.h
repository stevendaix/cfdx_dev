#pragma once

#include "cfdx/core/linalg/preconditioner.h"
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx::core {

// A reusable 6-field JFNK preconditioning adapter for variable-major
// [Ux, Uy, Uz, p, k, omega] unknowns, N cells per field. The first 4N
// component may be any independently prepared native COUPLED 4N Schur
// (BlockLocal / PCD / LSC / BFBT / SIMPLE / SIMPLEC, AMG, Braess-Sarazin,
// MGR, Schwarz). The last N+N components may be native AMG or other scalar
// preconditioners. Off-diagonal SST-momentum couplings remain in the true
// matrix-free F, and are deliberately omitted ONLY from this block-diagonal
// approximation. This is NOT by itself a monolithic production FV residual.
//
// The adapter borrows the three preconditioners. Its setup() receives an
// approximate assembled 6N CSR and calls each block's native setup() using
// a separately extracted CSR, preserving the block preconditioner contract.
class SixFieldBlockPreconditionerV8 final : public Preconditioner {
public:
    SixFieldBlockPreconditionerV8(
        std::size_t n_cells, Preconditioner& pressure_velocity,
        Preconditioner& k, Preconditioner& omega)
        : nc_(n_cells), up_(&pressure_velocity), k_(&k), omega_(&omega) {
        if (nc_ == 0 || nc_ > std::numeric_limits<std::size_t>::max() / 6)
            throw std::invalid_argument("V8 six-field block: invalid cell count");
        if (up_ == k_ || up_ == omega_ || k_ == omega_)
            throw std::invalid_argument("V8 six-field block: shared mutable block preconditioner");
    }
    bool setup(const SparseMatrix& A) override {
        ready_=false;
        if (!A.is_consistent() || A.n_rows()!=6*nc_ || A.n_cols()!=6*nc_)
            return false;
        const auto slice=[&](std::size_t offset,std::size_t dimension){
            SparseMatrix block(dimension,dimension);
            for(std::size_t row=0;row<dimension;++row) {
                const std::size_t global=row+offset;
                for(std::size_t at=A.row_offsets_data()[global];
                    at<A.row_offsets_data()[global+1];++at) {
                    const auto col=A.columns_data()[at];
                    if(col>=offset && col<offset+dimension)
                        block.push_back(row,col-offset,A.values_data()[at]);
                }
            }
            block.finalize();
            return block;
        };
        const auto velocity_pressure=slice(0,4*nc_);
        const auto k_matrix=slice(4*nc_,nc_);
        const auto omega_matrix=slice(5*nc_,nc_);
        ready_=up_->setup(velocity_pressure) &&
               k_->setup(k_matrix) && omega_->setup(omega_matrix);
        return ready_;
    }
    bool apply(const Vector& r,Vector& z) const override {
        if (!ready_ || r.size()!=6*nc_ || z.size()!=6*nc_ || !r.is_valid())
            return false;
        Vector rup(4*nc_,0.), rk(nc_,0.), rw(nc_,0.);
        for (std::size_t i=0;i<4*nc_;++i) rup(i)=r(i);
        for (std::size_t i=0;i<nc_;++i) {
            rk(i)=r(4*nc_+i);
            rw(i)=r(5*nc_+i);
        }
        Vector zup(4*nc_,0.), zk(nc_,0.), zw(nc_,0.);
        if(!up_->apply(rup,zup) || !k_->apply(rk,zk) ||
           !omega_->apply(rw,zw) ||
           !zup.is_valid() || !zk.is_valid() || !zw.is_valid())
            return false;
        // Publish atomically: a rejected scalar block must not partially
        // mutate the user-visible Krylov preconditioner image.
        for (std::size_t i=0;i<4*nc_;++i) z(i)=zup(i);
        for (std::size_t i=0;i<nc_;++i) {
            z(4*nc_+i)=zk(i);
            z(5*nc_+i)=zw(i);
        }
        return true;
    }
    const char* name() const override { return "SixFieldBlockV8(UP+K+Omega)"; }
    bool ready() const noexcept { return ready_; }
private:
    std::size_t nc_;
    Preconditioner* up_;
    Preconditioner* k_;
    Preconditioner* omega_;
    bool ready_=false;
};

} // namespace cfdx::core
