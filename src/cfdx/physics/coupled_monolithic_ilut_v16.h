#pragma once
// EXPERIMENTAL small-system fallback for diagnosing failed Schur/PCD outer
// FGMRES.  Pure C++ ILUT with bounded fill; the original monolithic saddle
// matrix and independent physical residual are NEVER modified.  This is not
// a scalable replacement for AMG/MGR or a default production preconditioner.
#include "cfdx/core/linalg/preconditioner.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <utility>
#include <vector>

namespace cfdx::physics {

class CoupledMonolithicIlutV16 final : public cfdx::core::Preconditioner {
public:
    // A bounded amount of fill in each row keeps memory predictable.  Entries
    // whose absolute value is insignificant relative to the original row
    // magnitude are discarded, except for the mandatory diagonal.
    explicit CoupledMonolithicIlutV16(double drop_tolerance=1e-5,
                                       std::size_t max_fill_factor=30)
      : drop_tolerance_(drop_tolerance), max_fill_factor_(max_fill_factor) {}

    bool setup(const cfdx::core::SparseMatrix& A) override {
        lower_.clear();upper_.clear();diagonal_.clear();regularized_pivots_=0;
        if (!A.is_consistent() || A.n_rows()!=A.n_cols() || A.n_rows()==0 ||
            A.n_rows()>4096 || !(drop_tolerance_>=0.) ||
            !std::isfinite(drop_tolerance_) || max_fill_factor_<1 ||
            max_fill_factor_>100) return false;
        const auto n=A.n_rows();
        lower_.resize(n);upper_.resize(n);diagonal_.resize(n);
        for (std::size_t i=0;i<n;++i) {
            std::map<std::size_t,long double> work;
            long double original_norm=0;
            const std::size_t row_entries=A.row_offsets_data()[i+1]-A.row_offsets_data()[i];
            for (std::size_t k=A.row_offsets_data()[i]; k<A.row_offsets_data()[i+1];++k) {
                const auto col=A.columns_data()[k];const auto value=A.values_data()[k];
                if (!std::isfinite(value) || col>=n) return false;
                work[col]+=value;
                original_norm=std::max(original_norm,std::abs(static_cast<long double>(value)));
            }
            if (!(original_norm>0)) return false;
            const long double drop=static_cast<long double>(drop_tolerance_)*original_norm;
            for (auto it=work.begin();it!=work.end() && it->first<i;++it) {
                const auto j=it->first;
                const long double wj=it->second;
                if (std::abs(wj)<=drop) {it->second=0.;continue;}
                const long double mult=wj/diagonal_[j];
                if (!std::isfinite(mult)) return false;
                it->second=mult;
                for (const auto& [col,u]:upper_[j]) {
                    const long double updated=work[col]-mult*u;
                    if (!std::isfinite(updated)) return false;
                    work[col]=updated;
                }
            }
            const auto diag_it=work.find(i);
            long double diag=diag_it!=work.end()?diag_it->second:0.0L;
            const long double pivot_floor = original_norm*1e-13L;
            if (std::abs(diag)<pivot_floor) {
                diag=diag<0 ? -pivot_floor : pivot_floor;
                ++regularized_pivots_;
            }
            if (!std::isfinite(diag) || diag==0.) return false;
            diagonal_[i]=static_cast<double>(diag);
            if (!std::isfinite(diagonal_[i]) || diagonal_[i]==0.) return false;
            using Entry=std::pair<std::size_t,double>;
            std::vector<Entry> lower, upper;
            for (const auto& [col,value]:work) {
                if (col==i) continue;
                const bool keep= col<i ? value!=0. : std::abs(value)>drop;
                if (!keep) continue;
                if (!std::isfinite(value) ||
                    std::abs(value)>std::numeric_limits<double>::max()) return false;
                (col<i ? lower:upper).emplace_back(col,static_cast<double>(value));
            }
            const auto per_row_limit = std::max<std::size_t>(8,row_entries*max_fill_factor_);
            const auto trim=[&](std::vector<Entry>& rows) {
                if (rows.size()>per_row_limit) {
                    std::nth_element(rows.begin(),rows.begin()+per_row_limit,rows.end(),
                        [](const Entry& a,const Entry& b) {
                            return std::abs(a.second)>std::abs(b.second);
                        });
                    rows.resize(per_row_limit);
                }
                std::sort(rows.begin(),rows.end(),
                    [](const Entry& a,const Entry& b){return a.first<b.first;});
            };
            trim(lower);trim(upper);
            lower_[i]=std::move(lower);upper_[i]=std::move(upper);
        }
        return true;
    }

    bool apply(const cfdx::core::Vector& rhs,cfdx::core::Vector& result) const override {
        const auto n=diagonal_.size();
        if (n==0 || rhs.size()!=n || result.size()!=n) return false;
        std::vector<double> intermediate(n,0.);
        for (std::size_t i=0;i<n;++i) {
            long double value=rhs(i);
            if (!std::isfinite(value)) return false;
            for (const auto& [j,l]:lower_[i]) value-=static_cast<long double>(l)*intermediate[j];
            if (!std::isfinite(value) ||
                std::abs(value)>std::numeric_limits<double>::max()) return false;
            intermediate[i]=static_cast<double>(value);
        }
        for (std::size_t ii=n;ii-->0;) {
            long double value=intermediate[ii];
            for (const auto& [j,u]:upper_[ii]) value-=static_cast<long double>(u)*result(j);
            value/=diagonal_[ii];
            if (!std::isfinite(value) ||
                std::abs(value)>std::numeric_limits<double>::max()) return false;
            result(ii)=static_cast<double>(value);
        }
        return true;
    }
    const char* name() const override { return "V16-Monolithic-ILUT-Experimental"; }
    std::size_t regularized_pivots() const noexcept {return regularized_pivots_;}
private:
    double drop_tolerance_;
    std::size_t max_fill_factor_;
    std::size_t regularized_pivots_=0;
    std::vector<double> diagonal_;
    std::vector<std::vector<std::pair<std::size_t,double>>> lower_,upper_;
};
} // namespace cfdx::physics
