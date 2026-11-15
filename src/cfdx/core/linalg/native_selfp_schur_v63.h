#pragma once
// V63 RESEARCH ONLY. Additional, explicitly named seventh Schur option.
// Nothing here changes PCD/LSC/BFBt/BlockLocal/SIMPLE/SIMPLEC.
// S~ = C - D diag(Auu)^{-1} G is built from the *actual* 4N CSR including
// its pressure gauge, then ILU(0) is applied to Auu and S~ in a full
// block-LDU preconditioner. This is NOT the true Schur inverse and never
// claims nonlinear or BFS physical qualification.

#include "cfdx/core/linalg/advanced_preconditioners.h"
#include "cfdx/core/linalg/coupled_block_extraction.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::core {

class NativeSelfpILUSchurV63 final : public Preconditioner {
public:
    explicit NativeSelfpILUSchurV63(std::size_t cells) : n_(cells), nu_(3*cells) {}

    bool setup(const SparseMatrix& matrix) override {
        clear();
        if (!n_ || n_ > std::numeric_limits<std::size_t>::max()/4 ||
            !matrix.is_consistent() || matrix.n_rows()!=4*n_ ||
            matrix.n_cols()!=4*n_)
            return fail("invalid 4N physical coupled CSR");
        const auto u=extract_coupled_block(matrix,0,0,nu_,nu_);
        const auto g=extract_coupled_block(matrix,0,1,nu_,n_);
        const auto d=extract_coupled_block(matrix,1,0,n_,nu_);
        const auto c=extract_coupled_block(matrix,1,1,n_,n_);
        std::vector<double> inv(nu_,0.0);
        for (std::size_t row=0;row<nu_;++row) {
            bool diagonal=false;
            double rowscale=0.0,value=0.0;
            for(std::size_t k=u.row_offsets_data()[row];
                k<u.row_offsets_data()[row+1];++k) {
                rowscale=std::max(rowscale,std::abs(u.values_data()[k]));
                if(u.columns_data()[k]==row) {
                    diagonal=true;value=u.values_data()[k];
                }
            }
            if(!diagonal || !std::isfinite(value) ||
               std::abs(value)<=128.*std::numeric_limits<double>::epsilon()*
                                        std::max(rowscale,1e-300))
                return fail("singular or absent velocity diagonal; no regularization");
            inv[row]=1./value;
            if(!std::isfinite(inv[row]))return fail("nonfinite inverse diagonal");
        }
        SparseMatrix schur(n_,n_);
        std::size_t terms=0;
        for(std::size_t row=0;row<n_;++row) {
            std::map<std::size_t,long double> entries;
            for(std::size_t k=c.row_offsets_data()[row];
                k<c.row_offsets_data()[row+1];++k)
                entries[c.columns_data()[k]]+=c.values_data()[k];
            for(std::size_t k=d.row_offsets_data()[row];
                k<d.row_offsets_data()[row+1];++k) {
                const auto velocity=d.columns_data()[k];
                const long double scale=static_cast<long double>(d.values_data()[k])*inv[velocity];
                for(std::size_t j=g.row_offsets_data()[velocity];
                    j<g.row_offsets_data()[velocity+1];++j) {
                    entries[g.columns_data()[j]]-=scale*g.values_data()[j];
                    ++terms;
                }
            }
            // The graph and gauge come from the real A. No artificially
            // added pin, sign flip or zero-mode suppression.
            for(const auto& item:entries) {
                if(!std::isfinite(item.second) ||
                    std::abs(item.second)>std::numeric_limits<double>::max())
                    return fail("nonfinite Schur coefficient");
                if(item.second!=0.0L)
                    schur.push_back(row,item.first,static_cast<double>(item.second));
            }
        }
        schur.finalize();
        ILU0Preconditioner velocity_ilu,pressure_ilu;
        if(!velocity_ilu.setup(u))return fail("native ILU0(Auu) pivot failure");
        if(!pressure_ilu.setup(schur))return fail("native ILU0(selfp) pivot failure");
        // Publish only after both factorizations succeed.
        u_=u;g_=g;d_=d;schur_=std::move(schur);
        velocity_ilu_=std::move(velocity_ilu);
        pressure_ilu_=std::move(pressure_ilu);
        triple_terms_=terms;
        ready_=true;
        return true;
    }

    bool apply(const Vector& rhs,Vector& result) const override {
        if(!ready_ || rhs.size()!=4*n_ || result.size()!=4*n_ || !rhs.is_valid())
            return false;
        Vector ru(nu_),rp(n_),u0(nu_),pressure_rhs(n_),p(n_),gu(nu_),u(nu_);
        for(std::size_t i=0;i<nu_;++i)ru(i)=rhs(i);
        for(std::size_t i=0;i<n_;++i)rp(i)=rhs(nu_+i);
        if(!velocity_ilu_.apply(ru,u0))return false;
        d_.matvec_into(u0,pressure_rhs);
        for(std::size_t i=0;i<n_;++i)pressure_rhs(i)=rp(i)-pressure_rhs(i);
        if(!pressure_ilu_.apply(pressure_rhs,p))return false;
        g_.matvec_into(p,gu);
        for(std::size_t i=0;i<nu_;++i)ru(i)-=gu(i);
        if(!velocity_ilu_.apply(ru,u))return false;
        if(!u.is_valid() || !p.is_valid())return false;
        for(std::size_t i=0;i<nu_;++i)result(i)=u(i);
        for(std::size_t i=0;i<n_;++i)result(nu_+i)=p(i);
        return result.is_valid();
    }

    const char* name() const override {return "NativeSelfp-ILU0-V63[experimental]";}
    const SparseMatrix& schur_matrix() const noexcept {return schur_;}
    std::size_t triple_terms() const noexcept {return triple_terms_;}
    bool ready() const noexcept {return ready_;}
    const std::string& last_error() const noexcept {return error_;}

private:
    bool fail(const char* message) {clear();error_=message;return false;}
    void clear() {ready_=false;triple_terms_=0;u_=SparseMatrix{};
                  g_=SparseMatrix{};d_=SparseMatrix{};schur_=SparseMatrix{};}
    std::size_t n_=0,nu_=0,triple_terms_=0;
    bool ready_=false;
    std::string error_;
    SparseMatrix u_,g_,d_,schur_;
    ILU0Preconditioner velocity_ilu_,pressure_ilu_;
};
} // namespace cfdx::core
