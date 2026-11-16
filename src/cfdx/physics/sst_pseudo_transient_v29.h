#pragma once
// CFDX V29: local diagonal pseudo-time continuation for STEADY SST.
// This modifies only the frozen systems used to produce a CANDIDATE. The
// accepted state must pass the unshifted, reassembled physical FV residual.
// A_ii/CFL is algebraically rho*V/dtau_i with
// dtau_i = CFL*rho*V/A_ii. This is a diagonal time-scale approximation, not
// the full convective/spectral-radius local time-step of SU2 or Fluent.
#include "cfdx/physics/finite_volume_transport.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct SSTV29PseudoTimeState {
    double cfl=0.5;
    static constexpr double min_cfl=1./128.;
    static constexpr double max_cfl=16.;
    unsigned selected_steps=0, rejected_steps=0;
    void selected(double old_merit, double new_merit) {
        if (!(old_merit>0.) || !std::isfinite(old_merit) ||
            !(new_merit>=0.) || !std::isfinite(new_merit) ||
            !(new_merit<old_merit))
            throw std::invalid_argument("V29 pseudo-time requires actual physical contraction");
        // Change CFL only AFTER a candidate wins the final FV selection.
        const double fractional_gain=(old_merit-new_merit)/old_merit;
        cfl=std::min(max_cfl,cfl*(fractional_gain>0.1 ? 2.0 : 1.25));
        ++selected_steps;
    }
    void rejected() noexcept {
        cfl=std::max(min_cfl,0.5*cfl);
        ++rejected_steps;
    }
};

struct SSTV29PseudoTimeAudit {
    std::size_t shifted_rows=0, skipped_wall_rows=0;
    double added_diagonal_sum=0.;
};

inline SSTV29PseudoTimeAudit sst_add_local_pseudo_time_v29(
    ScalarEquation& equation, const cfdx::core::Vector& previous,
    double cfl, const std::vector<unsigned char>& fixed_owner_mask={})
{
    const std::size_t n=previous.size();
    if (!n || equation.rhs.size()!=n || equation.matrix.n_rows()!=n ||
        equation.matrix.n_cols()!=n || equation.diagonal.size()!=n ||
        !(cfl>=SSTV29PseudoTimeState::min_cfl &&
          cfl<=SSTV29PseudoTimeState::max_cfl) ||
        !std::isfinite(cfl) ||
        (!fixed_owner_mask.empty() && fixed_owner_mask.size()!=n))
        throw std::invalid_argument("V29 pseudo-time: malformed system/CFL/mask");
    const auto* ro=equation.matrix.row_offsets_data();
    const auto* col=equation.matrix.columns_data();
    std::vector<std::size_t> diagonal_entry(n);
    std::vector<double> mass(n,0.);
    SSTV29PseudoTimeAudit audit;
    // Check entire update BEFORE mutating the equation (transactional).
    for (std::size_t i=0;i<n;++i) {
        if (!std::isfinite(previous(i)))
            throw std::invalid_argument("V29 pseudo-time: nonfinite previous value");
        if (!fixed_owner_mask.empty() && fixed_owner_mask[i]) {
            ++audit.skipped_wall_rows;
            continue;
        }
        bool found=false;
        double diagonal=0.;
        for (auto j=ro[i];j<ro[i+1];++j)
            if (col[j]==i) {found=true;diagonal=equation.matrix.values_data()[j];diagonal_entry[i]=j;break;}
        if (!found || !(diagonal>0.) || !std::isfinite(diagonal) ||
            !(equation.diagonal[i]>0.) || !std::isfinite(equation.diagonal[i]))
            throw std::invalid_argument("V29 pseudo-time: no positive row diagonal");
        mass[i]=diagonal/cfl;
        if (!(mass[i]>0.) || !std::isfinite(mass[i]) ||
            !std::isfinite(mass[i]*previous(i)) ||
            !std::isfinite(equation.rhs(i)+mass[i]*previous(i)) ||
            !std::isfinite(diagonal+mass[i]))
            throw std::invalid_argument("V29 pseudo-time: unsafe pseudo-mass");
        ++audit.shifted_rows;
        audit.added_diagonal_sum+=mass[i];
    }
    if (!std::isfinite(audit.added_diagonal_sum))
        throw std::invalid_argument("V29 pseudo-time: coefficient overflow");
    auto staged=equation;
    auto* a=staged.matrix.values_data();
    for (std::size_t i=0;i<n;++i) {
        if (!mass[i]) continue;
        a[diagonal_entry[i]]+=mass[i];
        staged.diagonal[i]+=mass[i];
        staged.rhs(i)+=mass[i]*previous(i);
    }
    equation=std::move(staged);
    return audit;
}
} // namespace cfdx::physics
