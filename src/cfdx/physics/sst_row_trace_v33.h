#pragma once
// V33 research/diagnostic only. No production turbulence model or solver
// coefficient is modified by this module. All indices are zero-based.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

inline std::vector<std::size_t> sst_v33_trace_cells(
    const char* specification, std::size_t cell_count)
{
    if (specification == nullptr) return {};
    // Each row audit reassembles FV equations and scans faces. Never run this
    // exhaustive tracer on the 20,540-cell BFS; use a separate top-K collector.
    if (cell_count == 0 || cell_count > 4096)
        throw std::invalid_argument("V33 SST trace: mesh exceeds diagnostic limit");
    const std::string input(specification);
    if (input.empty())
        throw std::invalid_argument("V33 SST trace: empty cell list");
    std::vector<std::size_t> cells;
    std::size_t pos = 0;
    while (pos < input.size()) {
        const auto end = input.find(',',pos);
        const auto stop = end == std::string::npos ? input.size() : end;
        if (stop == pos || cells.size() == 8)
            throw std::invalid_argument("V33 SST trace: invalid cell count/list");
        std::size_t index = 0;
        for (auto i = pos; i < stop; ++i) {
            const auto ch = input[i];
            if (ch < '0' || ch > '9')
                throw std::invalid_argument("V33 SST trace: non-digit cell ID");
            const auto digit = static_cast<std::size_t>(ch - '0');
            if (index > (std::numeric_limits<std::size_t>::max()-digit)/10)
                throw std::invalid_argument("V33 SST trace: overflowing cell ID");
            index = index*10+digit;
        }
        if (index >= cell_count ||
            std::find(cells.begin(),cells.end(),index)!=cells.end())
            throw std::invalid_argument("V33 SST trace: missing or duplicate cell");
        cells.push_back(index);
        if (end == std::string::npos) break;
        pos = end+1;
        if (pos == input.size())
            throw std::invalid_argument("V33 SST trace: trailing comma");
    }
    return cells;
}

struct SSTV33SourceBranches {
    bool omega_viscosity_branch = false;
    bool production_limited = false;
    bool wall_production = false;
    bool k_on_floor = false;
    bool omega_on_floor = false;
    bool close_to_viscosity_switch = false;
    bool close_to_production_switch = false;
    double unconstrained_production = 0.0;
    double limited_production = 0.0;
    double turbulent_viscosity = 0.0;
    // Frozen F1,F2,strain,gradients and wall G: source-only Jacobian
    // PER UNIT VOLUME. Multiply by the cell volume to compare with
    // the FV-integrated discrete cell budgets. Neither is the Jacobian
    // of the complete FV residual (faces, blending and wall targets vary).
    bool frozen_source_derivative_valid = false;
    double d_rk_dlog_k_frozen_source = 0.0;
    double d_rk_dlog_omega_frozen_source = 0.0;
};

// Branch audit of the *existing* SST source algebra, NOT an alternative
// production model or an approximation of the fully reassembled Jacobian.
// F1/F2, strain, wall G and gradients are supplied from the current FV state.
inline SSTV33SourceBranches sst_v33_source_branches(
    double k, double omega, double strain, double f2,
    double wall_g_kinematic, double density, double a1,
    double beta_star, double production_limiter, double k_min,
    double omega_min)
{
    if (!std::isfinite(k) || !std::isfinite(omega) ||
        !std::isfinite(strain) || !std::isfinite(f2) ||
        !std::isfinite(wall_g_kinematic) || !std::isfinite(density) ||
        !std::isfinite(a1) || !std::isfinite(beta_star) ||
        !std::isfinite(production_limiter) || !std::isfinite(k_min) ||
        !std::isfinite(omega_min) || k_min<=0 || omega_min<=0 ||
        strain<0 || f2<0 || f2>1 || density<=0 || a1<=0 ||
        beta_star<=0 || production_limiter<=0)
        throw std::invalid_argument("V33 SST branch: invalid source state");
    SSTV33SourceBranches out;
    out.k_on_floor = k <= k_min;
    out.omega_on_floor = omega <= omega_min;
    const auto ki=std::max(k,k_min), wi=std::max(omega,omega_min);
    const auto omega_denom=a1*wi;
    const auto strain_denom=strain*f2;
    const auto denominator=std::max(omega_denom,strain_denom);
    if (!(denominator>0.) || !std::isfinite(denominator))
        throw std::invalid_argument("V33 SST branch: invalid viscosity branch");
    out.omega_viscosity_branch=omega_denom >= strain_denom;
    out.close_to_viscosity_switch=
        std::abs(omega_denom-strain_denom) <=
            1e-6*std::max({omega_denom,strain_denom,1e-30});
    out.turbulent_viscosity=a1*ki/denominator;
    out.wall_production=wall_g_kinematic>=0;
    const auto raw_per_mass=out.wall_production
        ? wall_g_kinematic : out.turbulent_viscosity*strain*strain;
    out.unconstrained_production=density*raw_per_mass;
    const auto ceiling=production_limiter*beta_star*density*ki*wi;
    out.production_limited=out.unconstrained_production>ceiling;
    out.limited_production=std::min(out.unconstrained_production,ceiling);
    out.close_to_production_switch=
        std::abs(out.unconstrained_production-ceiling) <=
            1e-6*std::max({std::abs(out.unconstrained_production),ceiling,1e-30});
    // At branch switches the classical central derivative need not exist.
    // Do not supply a fabricated preconditioner coefficient there.
    out.frozen_source_derivative_valid = !out.k_on_floor &&
        !out.omega_on_floor && !out.close_to_viscosity_switch &&
        !out.close_to_production_switch;
    if (out.frozen_source_derivative_valid) {
        const double destruction=density*beta_star*ki*wi;
        const double dpk_dlog_k=out.production_limited
            ? out.limited_production
            : (out.wall_production ? 0.0 : out.limited_production);
        const double dpk_dlog_omega=out.production_limited
            ? out.limited_production
            : (out.wall_production ? 0.0 :
                (out.omega_viscosity_branch ? -out.limited_production : 0.0));
        out.d_rk_dlog_k_frozen_source=destruction-dpk_dlog_k;
        out.d_rk_dlog_omega_frozen_source=destruction-dpk_dlog_omega;
    }
    if (!std::isfinite(out.d_rk_dlog_k_frozen_source) ||
        !std::isfinite(out.d_rk_dlog_omega_frozen_source))
        throw std::invalid_argument("V33 SST branch: derivative overflow");
    if (!std::isfinite(out.turbulent_viscosity) ||
        !std::isfinite(out.unconstrained_production) ||
        !std::isfinite(out.limited_production))
        throw std::invalid_argument("V33 SST branch: source overflow");
    return out;
}

inline bool sst_v33_branch_changed(const SSTV33SourceBranches& a,
                                   const SSTV33SourceBranches& b)
{
    return a.omega_viscosity_branch != b.omega_viscosity_branch ||
           a.production_limited != b.production_limited ||
           a.wall_production != b.wall_production ||
           a.k_on_floor != b.k_on_floor ||
           a.omega_on_floor != b.omega_on_floor;
}

} // namespace cfdx::physics
