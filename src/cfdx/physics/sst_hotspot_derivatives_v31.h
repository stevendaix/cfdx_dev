#pragma once
// Diagnostic V31: differentiate the ACTUAL reassembled finite-volume budgets.
// No approximate SST formula is used by this helper. Calls to the evaluator
// must be side-effect-free and use the same wall/face refresh as production.
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>
#include <algorithm>

namespace cfdx::physics {
struct SSTV31RowComponents {
    // A*x-b, separated in the sign convention of the real FV budget.
    std::array<double, 6> parts{}; // convection, diffusion, production, destruction, cross, residual
    double closure_error=0.;
    bool constrained=false;
};
struct SSTV31Derivative {
    std::size_t cell=0;
    double perturbation=0.;
    std::array<double,6> d_dlog{};
    double maximum_relative_closure=0.;
};
using SSTV31BudgetEvaluator = std::function<SSTV31RowComponents(const std::vector<double>&)>;
inline void sst_v31_validate_components(const SSTV31RowComponents& v) {
    if (v.constrained) throw std::invalid_argument("V31 hotspot: constrained row is not a transport PDE");
    double sum=0.;
    for (std::size_t i=0;i<5;++i) {
        if (!std::isfinite(v.parts[i])) throw std::invalid_argument("V31 hotspot: nonfinite source budget");
        sum+=v.parts[i];
    }
    if (!std::isfinite(v.parts[5]) || !std::isfinite(v.closure_error) ||
        std::abs(sum-v.parts[5])>1e-8*std::max({1.,std::abs(sum),std::abs(v.parts[5])}))
        throw std::invalid_argument("V31 hotspot: FV row budget does not close");
}
inline SSTV31Derivative sst_v31_log_derivative(
    const std::vector<double>& state, std::size_t variable,
    double lower_bound, double h, const SSTV31BudgetEvaluator& evaluator)
{
    if (state.empty() || variable>=state.size() ||
        !(lower_bound>0.) || !(h>1e-8 && h<0.1) ||
        !std::isfinite(state[variable]) || !(state[variable]>lower_bound))
        throw std::invalid_argument("V31 hotspot: invalid perturbation");
    auto p=state,m=state;
    p[variable]=state[variable]*std::exp(h);
    m[variable]=state[variable]*std::exp(-h);
    if (!(m[variable]>lower_bound) || !std::isfinite(p[variable]))
        throw std::invalid_argument("V31 hotspot: step would cross positivity bound");
    const auto plus=evaluator(p),minus=evaluator(m);
    sst_v31_validate_components(plus);
    sst_v31_validate_components(minus);
    SSTV31Derivative out;
    out.cell=variable;
    out.perturbation=h;
    for (std::size_t i=0;i<out.d_dlog.size();++i) {
        out.d_dlog[i]=(plus.parts[i]-minus.parts[i])/(2.*h);
        if (!std::isfinite(out.d_dlog[i]))
            throw std::invalid_argument("V31 hotspot: nonfinite derivative");
    }
    double sum=0.;
    for (std::size_t i=0;i<5;++i) sum+=out.d_dlog[i];
    out.maximum_relative_closure=std::abs(sum-out.d_dlog[5])/
        std::max({1.,std::abs(sum),std::abs(out.d_dlog[5])});
    if (out.maximum_relative_closure>1e-8)
        throw std::invalid_argument("V31 hotspot: differentiated FV budget does not close");
    return out;
}
} // namespace cfdx::physics
