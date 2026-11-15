@@ -0,0 +1,100 @@
#pragma once
// CFDX V89 - EXPERIMENTAL, opt-in, BFS SIMPLE/SST only.
// A trial is recomputed with a smaller nonlinear step when rejected;
// NO clipping of p, U, Rhie-Chow face flux or SST and NO change to final gates.
// Pressure decisions use ranges of pressure DIFFERENCES and are invariant to
// independent additive gauge constants applied to successive pressure fields.
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
namespace cfdx::physics {
struct BfsPressureTransversePreV89 {
    bool may_continue_sst=false;
    double transverse_increment_over_inlet=0.0;
    double pressure_increment_span_over_qin=0.0;
    std::string reason;
};
inline BfsPressureTransversePreV89 check_bfs_pressure_transverse_pre_v89(
    double delta_u_y_linf, double delta_pressure_span, double inlet_speed,
    double density, bool baseline_state,
    double transverse_limit=0.15, double pressure_limit_q=16.0)
{
    if (!(std::isfinite(delta_u_y_linf) && delta_u_y_linf>=0.0 &&
          std::isfinite(delta_pressure_span) && delta_pressure_span>=0.0 &&
          std::isfinite(inlet_speed) && inlet_speed>0.0 &&
          std::isfinite(density) && density>0.0 &&
          std::isfinite(transverse_limit) && transverse_limit>0.0 &&
          std::isfinite(pressure_limit_q) && pressure_limit_q>0.0))
        throw std::invalid_argument("V89 pre-SST pressure/transverse trust: invalid input");
    const double q=0.5*density*inlet_speed*inlet_speed;
    BfsPressureTransversePreV89 result;
    result.transverse_increment_over_inlet=delta_u_y_linf/inlet_speed;
    result.pressure_increment_span_over_qin=delta_pressure_span/q;
    if (!std::isfinite(result.transverse_increment_over_inlet) ||
        !std::isfinite(result.pressure_increment_span_over_qin))
        throw std::invalid_argument("V89 pre-SST pressure/transverse trust: overflow");
    // Startup is a baseline, NOT an accredited physical result; on an
    // initial uniform U/p state there is no meaningful prior pressure scale.
    if (baseline_state) {
        result.may_continue_sst=true;
        result.reason="startup_baseline_not_certified";
    } else if (result.pressure_increment_span_over_qin>pressure_limit_q) {
        result.reason="gauge_invariant_pressure_increment";
    } else if (result.transverse_increment_over_inlet>transverse_limit) {
        result.reason="transverse_momentum_increment";
    } else {
        result.may_continue_sst=true;
        result.reason="pressure_transverse_preflight_not_certified";
    }
    return result;
}
struct BfsPressureTransversePostV89 {
    bool may_commit=false;
    double pressure_span_over_qin=0.0;
    double transverse_momentum_relative=0.0;
    std::string reason;
};
inline BfsPressureTransversePostV89 check_bfs_pressure_transverse_post_v89(
    double pressure_span, double baseline_pressure_span, bool have_baseline,
    double transverse_momentum_relative, double previous_transverse_relative,
    bool have_previous, bool fully_reassembled, bool sst_bounded,
    double continuity, double continuity_tolerance, double inlet_speed,
    double density, double max_span_vs_initial=3.0,
    double max_transverse_residual_growth=1.5)
{
    if (!(std::isfinite(pressure_span) && pressure_span>=0.0 &&
          std::isfinite(baseline_pressure_span) && baseline_pressure_span>=0.0 &&
          std::isfinite(transverse_momentum_relative) && transverse_momentum_relative>=0.0 &&
          std::isfinite(previous_transverse_relative) && previous_transverse_relative>=0.0 &&
          std::isfinite(continuity) && continuity>=0.0 &&
          std::isfinite(continuity_tolerance) && continuity_tolerance>0.0 &&
          std::isfinite(inlet_speed) && inlet_speed>0.0 &&
          std::isfinite(density) && density>0.0 &&
          std::isfinite(max_span_vs_initial) && max_span_vs_initial>=1.0 &&
          std::isfinite(max_transverse_residual_growth) && max_transverse_residual_growth>=1.0))
        throw std::invalid_argument("V89 full physical trust: invalid input");
    const double q=0.5*density*inlet_speed*inlet_speed;
    BfsPressureTransversePostV89 result;
    result.pressure_span_over_qin=pressure_span/q;
    result.transverse_momentum_relative=transverse_momentum_relative;
    if (!std::isfinite(result.pressure_span_over_qin))
        throw std::invalid_argument("V89 full physical trust: overflow");
    if (!fully_reassembled || !sst_bounded) {
        result.reason="missing_full_U_p_k_omega_certificate";
    } else if (continuity>continuity_tolerance) {
        result.reason="rhie_chow_conservation_not_closed";
    } else if (have_baseline &&
        pressure_span>max_span_vs_initial*std::max(baseline_pressure_span,q)) {
        result.reason="cumulative_gauge_invariant_pressure_span";
    } else if (have_previous &&
        transverse_momentum_relative >
            max_transverse_residual_growth*std::max(previous_transverse_relative,1e-12)) {
        result.reason="transverse_momentum_residual_growth";
    } else {
        result.may_commit=true;
        result.reason="nonlinear_progress_only_not_converged";
    }
    return result;
}
} // namespace cfdx::physics
