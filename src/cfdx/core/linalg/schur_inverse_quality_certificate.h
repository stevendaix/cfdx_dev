#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::core {

// Diagnostic only. On the EXACT same assembled, already gauge-fixed blocks,
// measure the quality of an approximate inverse B~=S^-1 by evaluating
//
//    e = (C - D Auu^-1 G) (B~ r) - r.
//
// No dense Schur matrix is formed. The cost is one sufficiently accurate
// velocity solve per pressure right-hand side. This permits apples-to-apples
// comparisons of PCD/LSC/BFBt/SIMPLE on captured systems with nonzero flux.
// It certifies algebraic inverse-action quality, NOT nonlinear conservation.
struct SchurInverseQualityControls {
    double relative_tolerance = 1e-7;
    double absolute_tolerance = 1e-10;
    // Precision of the independent Auu inverse used to form S*z.
    double velocity_relative_tolerance = 1e-9;
    double velocity_absolute_tolerance = 1e-12;
};

struct SchurInverseQualityCertificate {
    bool complete = false;
    bool passed = false;
    std::size_t probes_checked = 0;
    std::size_t worst_probe = 0;
    double worst_schur_defect = 0.0;
    double worst_relative_defect = 0.0;
    double worst_tolerance_ratio = 0.0;
    double worst_velocity_solve_defect = 0.0;
};

using SchurReferenceVelocitySolve =
    std::function<bool(const Vector& rhs, Vector& solution)>;

inline SchurInverseQualityCertificate certify_schur_inverse_action(
    const BlockOperator& blocks,
    const SchurApproximation& approximation,
    const SchurReferenceVelocitySolve& solve_auu,
    const std::vector<Vector>& pressure_rhs,
    const SchurInverseQualityControls& controls = {})
{
    blocks.validate();
    const auto valid_tolerance = [](double value) {
        return std::isfinite(value) && value >= 0.0;
    };
    if (blocks.velocity_size() == 0 || blocks.pressure_size() == 0 ||
        !solve_auu || pressure_rhs.empty() ||
        approximation.action() != SchurAction::InverseOperator ||
        !valid_tolerance(controls.relative_tolerance) ||
        !valid_tolerance(controls.absolute_tolerance) ||
        !valid_tolerance(controls.velocity_relative_tolerance) ||
        !valid_tolerance(controls.velocity_absolute_tolerance) ||
        (controls.relative_tolerance == 0.0 && controls.absolute_tolerance == 0.0) ||
        (controls.velocity_relative_tolerance == 0.0 &&
         controls.velocity_absolute_tolerance == 0.0))
        throw std::invalid_argument("Schur inverse quality: invalid inputs or action");

    SchurInverseQualityCertificate report;
    const std::size_t nu = blocks.velocity_size();
    const std::size_t np = blocks.pressure_size();
    const auto norm = [](const Vector& value) {
        double n = 0.0;
        for (std::size_t i = 0; i < value.size(); ++i)
            n = std::hypot(n, value(i));
        return n;
    };
    for (std::size_t sample = 0; sample < pressure_rhs.size(); ++sample) {
        const Vector& rhs = pressure_rhs[sample];
        if (rhs.size() != np || !rhs.is_valid()) return report;
        Vector p(np, 0.0);
        if (!approximation.apply(rhs, p) || p.size() != np || !p.is_valid())
            return report;

        Vector gp(nu, 0.0), auu_inv_gp(nu, 0.0);
        blocks.G().matvec_into(p, gp);
        if (!gp.is_valid() || !solve_auu(gp, auu_inv_gp) ||
            auu_inv_gp.size() != nu || !auu_inv_gp.is_valid()) return report;
        Vector solved_gp(nu, 0.0);
        blocks.Auu().matvec_into(auu_inv_gp, solved_gp);
        if (!solved_gp.is_valid()) return report;
        for (std::size_t i = 0; i < nu; ++i)
            solved_gp(i) -= gp(i);
        const double velocity_defect = norm(solved_gp);
        const double velocity_bound = controls.velocity_absolute_tolerance +
            controls.velocity_relative_tolerance * norm(gp);
        if (!std::isfinite(velocity_defect) || velocity_defect > velocity_bound)
            return report; // independent Schur oracle is not precise enough

        Vector cp(np, 0.0), dy(np, 0.0);
        blocks.C().matvec_into(p, cp);
        blocks.D().matvec_into(auu_inv_gp, dy);
        if (!cp.is_valid() || !dy.is_valid()) return report;
        Vector defect(np, 0.0);
        for (std::size_t i = 0; i < np; ++i)
            defect(i) = cp(i) - dy(i) - rhs(i);
        if (!defect.is_valid()) return report;
        const double absolute = norm(defect);
        const double rhs_norm = norm(rhs);
        const double relative = absolute / std::max(rhs_norm, controls.absolute_tolerance);
        const double bound = controls.absolute_tolerance +
            controls.relative_tolerance * rhs_norm;
        const double ratio = absolute / bound;
        if (!std::isfinite(ratio) || !std::isfinite(relative)) return report;
        if (ratio > report.worst_tolerance_ratio) {
            report.worst_probe = sample;
            report.worst_tolerance_ratio = ratio;
        }
        report.worst_schur_defect = std::max(report.worst_schur_defect, absolute);
        report.worst_relative_defect = std::max(report.worst_relative_defect, relative);
        report.worst_velocity_solve_defect =
            std::max(report.worst_velocity_solve_defect, velocity_defect);
        ++report.probes_checked;
    }
    report.complete = report.probes_checked == pressure_rhs.size();
    report.passed = report.complete && report.worst_tolerance_ratio <= 1.0;
    return report;
}

// Deterministic reference RHS: canonical axes plus a nontrivial alternating
// vector. This is only a probe generator, not a proof on the full space unless
// the entire pressure basis is used and the actions are linear and fixed.
inline std::vector<Vector> schur_inverse_quality_basis(
    std::size_t np, std::size_t max_basis = 0)
{
    if (np == 0) throw std::invalid_argument("Schur probe basis: empty pressure space");
    const std::size_t count = max_basis == 0 ? np : std::min(np, max_basis);
    std::vector<Vector> probes;
    probes.reserve(count + 1);
    for (std::size_t sample = 0; sample < count; ++sample) {
        const std::size_t col = count == 1 ? 0 : sample * (np - 1) / (count - 1);
        Vector rhs(np, 0.0);
        rhs(col) = 1.0;
        probes.push_back(std::move(rhs));
    }
    Vector alternating(np, 0.0);
    for (std::size_t i = 0; i < np; ++i)
        alternating(i) = (i % 2 == 0) ? 1.0 : -1.0;
    probes.push_back(std::move(alternating));
    return probes;
}

} // namespace cfdx::core
