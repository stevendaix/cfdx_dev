#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::physics {

// Frozen right preconditioner for the SST logarithmic Newton variables.
// Each cell stores [dFk/dlog(k), dFk/dlog(omega);
//                  dFw/dlog(k), dFw/dlog(omega)] AFTER the same fixed row
// scaling used for the Newton residual. F itself is never modified.
// This is a block-Jacobi *approximation*: off-cell diffusion, gradient,
// viscosity, wall G and wall-target derivatives remain in matrix-free Jv.
class SSTBlockPreconditionerV23 {
public:
    using Block = std::array<double,4>; // row major: kk,kw,wk,ww

    explicit SSTBlockPreconditionerV23(const std::vector<Block>& raw)
        : n_(raw.size()), raw_blocks_(raw) {
        if (n_==0) throw std::invalid_argument("V23 SST block: empty Jacobian");
        inverse_.reserve(n_);
        for (const auto& b:raw) {
            for (const double v:b)
                if (!std::isfinite(v))
                    throw std::invalid_argument("V23 SST block: nonfinite Jacobian");
            const double scale=std::max({std::abs(b[0]),std::abs(b[1]),
                                         std::abs(b[2]),std::abs(b[3])});
            if (!(scale>0.))
                throw std::invalid_argument("V23 SST block: empty local operator");
            const double a=b[0]/scale,c=b[1]/scale,d=b[2]/scale,e=b[3]/scale;
            const double det=a*e-c*d;
            // Reject tiny determinants relative to the product of row norms.
            // A singular local factor must NEVER silently replace the global
            // physical or matrix-free Jacobian. Use stable diagonal fallback.
            const double row_product=(std::abs(a)+std::abs(c)) *
                                     (std::abs(d)+std::abs(e));
            if (!(std::abs(det)>1e-8*row_product)) {
                ++diagonal_fallbacks_;
                // Identity is an invertible, bounded fallback; unlike an
                // invented inverse it can be used for both zero diagonals.
                inverse_.push_back({1.,0.,0.,1.});
                continue;
            }
            const double multiplier=(1.0/scale)/det;
            const Block inv{e*multiplier,-c*multiplier,
                            -d*multiplier,a*multiplier};
            bool valid=true;
            for (const double v:inv) valid &= std::isfinite(v);
            if (!valid) {
                ++diagonal_fallbacks_;
                inverse_.push_back({1.,0.,0.,1.});
            } else {
                if (b[1]!=0. || b[2]!=0.) ++coupled_blocks_;
                inverse_.push_back(inv);
            }
        }
    }

    std::vector<double> apply(const std::vector<double>& rhs) const {
        if (rhs.size()!=2*n_)
            throw std::invalid_argument("V23 SST block: invalid residual vector dimension");
        std::vector<double> direction(rhs.size());
        for (std::size_t i=0;i<n_;++i) {
            const auto& m=inverse_[i];
            const double rk=rhs[i],rw=rhs[n_+i];
            if (!std::isfinite(rk) || !std::isfinite(rw))
                throw std::invalid_argument("V23 SST block: nonfinite rhs");
            direction[i]=m[0]*rk+m[1]*rw;
            direction[n_+i]=m[2]*rk+m[3]*rw;
            if (!std::isfinite(direction[i]) || !std::isfinite(direction[n_+i]))
                throw std::domain_error("V23 SST block: inverse overflow");
        }
        return direction;
    }

    const std::vector<Block>& raw_blocks() const noexcept {return raw_blocks_;}
    const std::vector<Block>& inverse_blocks() const noexcept {return inverse_;}
    std::size_t cells() const noexcept {return n_;}
    std::size_t coupled_blocks() const noexcept {return coupled_blocks_;}
    std::size_t diagonal_fallbacks() const noexcept {return diagonal_fallbacks_;}

private:
    std::size_t n_ = 0;
    std::size_t coupled_blocks_ = 0;
    std::size_t diagonal_fallbacks_ = 0;
    std::vector<Block> inverse_;
    std::vector<Block> raw_blocks_;
};

} // namespace cfdx::physics
