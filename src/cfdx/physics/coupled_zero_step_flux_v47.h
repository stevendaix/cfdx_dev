#pragma once
// Diagnostic only. Different Rhie-Chow inverse momentum diagonals can map
// *identical U,p fields* to different face fluxes between Picard iterations.
// A true nonlinear line search must not compare such unmatched operators as
// if F(x_previous) were a fixed, unique function evaluation.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct CoupledZeroStepFluxAuditV47 {
    double previous_continuity = 0.0;
    double current_diagonal_zero_step_continuity = 0.0;
    double divergence_jump_l1 = 0.0;
    double relative_face_flux_change = 0.0;
    double maximum_face_flux_change = 0.0;
    std::size_t worst_cell = 0;
};

template<typename Owner,typename Neighbour,typename OldFlux,typename NewFlux>
CoupledZeroStepFluxAuditV47 audit_coupled_zero_step_flux_v47(
    std::size_t cells, std::size_t faces, Owner owner, Neighbour neighbour,
    OldFlux prior, NewFlux remapped) {
    if (cells == 0 || faces == 0)
        throw std::invalid_argument("V47 zero-step flux: empty connectivity");
    std::vector<double> div_old(cells,0.0), div_new(cells,0.0);
    double base=0.0, newbase=0.0, diff=0.0;
    CoupledZeroStepFluxAuditV47 result;
    for (std::size_t f=0;f<faces;++f) {
        const std::size_t o=owner(f);
        const auto n=neighbour(f);
        const double po=prior(f), pn=remapped(f);
        if (o>=cells || n< -1 || (n>=0 && static_cast<std::size_t>(n)>=cells) ||
            !std::isfinite(po) || !std::isfinite(pn))
            throw std::invalid_argument("V47 zero-step flux: corrupt face");
        div_old[o]+=po;
        div_new[o]+=pn;
        if (n>=0) {
            div_old[static_cast<std::size_t>(n)]-=po;
            div_new[static_cast<std::size_t>(n)]-=pn;
        }
        base+=std::abs(po); newbase+=std::abs(pn);
        diff+=std::abs(po-pn);
        result.maximum_face_flux_change=std::max(result.maximum_face_flux_change,std::abs(po-pn));
    }
    double jump=0.0, largest=0.0;
    for (std::size_t i=0;i<cells;++i) {
        const double delta=std::abs(div_new[i]-div_old[i]);
        result.previous_continuity+=std::abs(div_old[i]);
        result.current_diagonal_zero_step_continuity+=std::abs(div_new[i]);
        jump+=delta;
        if (delta>largest) { largest=delta;result.worst_cell=i; }
    }
    result.previous_continuity /= std::max(base,1e-300);
    result.current_diagonal_zero_step_continuity /= std::max(newbase,1e-300);
    result.divergence_jump_l1=jump/std::max(base,1e-300);
    result.relative_face_flux_change=diff/std::max(base,1e-300);
    for (const double x: {result.previous_continuity,
         result.current_diagonal_zero_step_continuity,
         result.divergence_jump_l1,result.relative_face_flux_change})
        if (!std::isfinite(x))
            throw std::invalid_argument("V47 zero-step flux: nonfinite result");
    return result;
}
} // namespace cfdx::physics
