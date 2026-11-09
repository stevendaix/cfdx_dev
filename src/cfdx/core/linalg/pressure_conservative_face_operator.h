#pragma once

// Independently evaluates the implicit FV pressure-correction operator from
// EXACT face owner/neighbour and coefficient arrays. Unlike the rounded CSR,
// zero-gradient interior rows carry NO invented diagonal row-sum leakage.
// This is not the full Rhie-Chow operator: deferred pressure gradients,
// pressure relaxation and momentum updates remain separate requirements.

#include "pressure_compensated_refinement.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::core {

struct PressureFaceRefinementReport {
    bool valid = false;
    bool implicit_face_converged = false;
    bool original_csr_converged = false;
    bool physical_fv_flux_certified = false;
    std::size_t accepted_refinements = 0;
    double initial_face_relative = std::numeric_limits<double>::infinity();
    double final_face_relative = std::numeric_limits<double>::infinity();
    double final_csr_relative = std::numeric_limits<double>::infinity();
};

// This certificate checks the actual owner/neighbour divergence of a base
// face flux after adding the high/low *implicit* pressure response.  Unlike
// A*x=b, it detects mismatched RHS/face ordering and a non-conservative source.
// It is deliberately not a certificate of deferred/non-orthogonal RC fluxes.
struct PressureFaceFluxBridgeReport {
    bool valid = false;
    bool rhs_matches_base_flux = false;
    bool implicit_linearized_continuity_converged = false;
    bool complete_rhie_chow_certified = false;
    double rhs_base_flux_relative = std::numeric_limits<double>::infinity();
    double implicit_flux_continuity_relative = std::numeric_limits<double>::infinity();
    double implicit_operator_residual_relative = std::numeric_limits<double>::infinity();
    double maximum_corrected_face_flux = 0.0;
};

// Decomposes the complete RC face residual into its already-cancelled
// implicit part and the pressure-gradient response absent from that solve.
// Spatial indices refer to native, unfiltered face/cell indexing.
struct PressureFullFluxGapReport {
    bool valid = false;
    double full_continuity_relative = std::numeric_limits<double>::infinity();
    double implicit_continuity_relative = std::numeric_limits<double>::infinity();
    double deferred_divergence_relative = std::numeric_limits<double>::infinity();
    double deferred_face_l2_relative = std::numeric_limits<double>::infinity();
    double maximum_abs_deferred_face = 0.0;
    std::size_t worst_deferred_face = 0;
    double maximum_abs_cell_divergence = 0.0;
    std::size_t worst_divergence_cell = 0;
    double maximum_abs_deferred_cell_divergence = 0.0;
    std::size_t worst_deferred_cell = 0;
};

class ConservativePressureFaceOperator {
public:
    struct Face {
        std::size_t owner = 0;
        std::int64_t neighbour = -1;
        double coefficient = 0.0; // positive integrated rho*V/A_P*orth_face/d
    };

    explicit ConservativePressureFaceOperator(std::size_t cells,
                                              std::vector<Face> faces)
        : cells_(cells), faces_(std::move(faces)) {
        if (!cells_ || faces_.empty())
            throw std::invalid_argument("pressure face operator: empty mesh");
        std::size_t active = 0;
        for (const Face& f : faces_) {
            if (f.owner >= cells_ || f.neighbour < -1 ||
                f.neighbour >= static_cast<std::int64_t>(cells_) ||
                f.neighbour == static_cast<std::int64_t>(f.owner) ||
                !std::isfinite(f.coefficient) || f.coefficient < 0.0)
                throw std::invalid_argument("pressure face operator: invalid face");
            if (f.coefficient > 0.0) ++active;
        }
        if (!active)
            throw std::invalid_argument("pressure face operator: no active pressure faces");
    }

    std::size_t size() const noexcept { return cells_; }
    std::size_t active_face_count() const noexcept {
        return static_cast<std::size_t>(std::count_if(
            faces_.begin(), faces_.end(),
            [](const Face& f) { return f.coefficient > 0.0; }));
    }

    // The face array must follow the same ordering as the supplied base flux.
    // The base flux must be the authoritative Rhie-Chow flux BEFORE the
    // pressure increment, with the identical frozen predictor, BC and rAU.
    // Internal pressure contributions appear once on each face, with the
    // opposite signs used in the two adjacent cell balances.  The caller
    // separately checks the complete deferred flux and velocity update.
    PressureFaceFluxBridgeReport audit_implicit_flux_bridge(
        const Vector& rhs, const Vector& high, const Vector& low,
        const std::vector<double>& base_face_flux, double tolerance) const {
        if (rhs.size() != cells_ || high.size() != cells_ ||
            low.size() != cells_ || base_face_flux.size() != faces_.size() ||
            !(tolerance > 0.0) || !std::isfinite(tolerance) ||
            std::numeric_limits<long double>::digits <=
                std::numeric_limits<double>::digits)
            throw std::invalid_argument("pressure face flux bridge: invalid inputs");
        std::vector<long double> before(cells_, 0.0L), after(cells_, 0.0L);
        long double norm_rhs = 0.0L;
        for (std::size_t c = 0; c < cells_; ++c) {
            if (!std::isfinite(rhs(c)) || !std::isfinite(high(c)) ||
                !std::isfinite(low(c)))
                throw std::domain_error("pressure face flux bridge: nonfinite cell input");
            norm_rhs += static_cast<long double>(rhs(c)) * rhs(c);
        }
        if (!(norm_rhs > 0.0L))
            throw std::invalid_argument("pressure face flux bridge: nonzero RHS required");
        PressureFaceFluxBridgeReport result;
        for (std::size_t face = 0; face < faces_.size(); ++face) {
            const auto& f = faces_[face];
            const double original = base_face_flux[face];
            if (!std::isfinite(original))
                throw std::domain_error("pressure face flux bridge: nonfinite base flux");
            const long double p_owner =
                static_cast<long double>(high(f.owner)) + low(f.owner);
            const long double p_neighbour = f.neighbour < 0 ? 0.0L :
                static_cast<long double>(high(static_cast<std::size_t>(f.neighbour))) +
                low(static_cast<std::size_t>(f.neighbour));
            const long double response =
                static_cast<long double>(f.coefficient) * (p_owner - p_neighbour);
            const long double corrected = static_cast<long double>(original) + response;
            if (!std::isfinite(corrected))
                throw std::domain_error("pressure face flux bridge: nonfinite corrected flux");
            result.maximum_corrected_face_flux = std::max(
                result.maximum_corrected_face_flux,
                static_cast<double>(std::abs(corrected)));
            before[f.owner] += original;
            after[f.owner] += corrected;
            if (f.neighbour >= 0) {
                const auto n = static_cast<std::size_t>(f.neighbour);
                before[n] -= original;
                after[n] -= corrected;
            }
        }
        long double b_error2 = 0.0L, phi_error2 = 0.0L;
        for (std::size_t c = 0; c < cells_; ++c) {
            const long double mismatch = before[c] + rhs(c);
            b_error2 += mismatch * mismatch;
            phi_error2 += after[c] * after[c];
        }
        result.rhs_base_flux_relative =
            static_cast<double>(std::sqrt(b_error2 / norm_rhs));
        result.implicit_flux_continuity_relative =
            static_cast<double>(std::sqrt(phi_error2 / norm_rhs));
        result.implicit_operator_residual_relative =
            residual(rhs, high, low).second;
        if (!std::isfinite(result.rhs_base_flux_relative) ||
            !std::isfinite(result.implicit_flux_continuity_relative))
            throw std::domain_error("pressure face flux bridge: invalid flux norm");
        result.valid = true;
        result.rhs_matches_base_flux = result.rhs_base_flux_relative <= tolerance;
        result.implicit_linearized_continuity_converged =
            result.rhs_matches_base_flux &&
            result.implicit_flux_continuity_relative <= tolerance &&
            result.implicit_operator_residual_relative <= tolerance;
        return result;
    }

    // Test the reconstructed COMPLETE RC flux against the implicit face
    // response. This is diagnostic, never a nonlinear or physical certificate.
    PressureFullFluxGapReport diagnose_full_flux_gap(
        const Vector& rhs, const Vector& high, const Vector& low,
        const std::vector<double>& base_face_flux,
        const std::vector<long double>& full_face_flux) const {
        if (rhs.size() != cells_ || high.size() != cells_ || low.size() != cells_ ||
            base_face_flux.size() != faces_.size() ||
            full_face_flux.size() != faces_.size())
            throw std::invalid_argument("pressure full flux gap: size mismatch");
        std::vector<long double> implicit(cells_,0.0L), full(cells_,0.0L),
            missing(cells_,0.0L);
        long double rhs2=0.0L, face_gap2=0.0L;
        for (std::size_t c=0;c<cells_;++c) {
            if (!std::isfinite(rhs(c)) || !std::isfinite(high(c)) ||
                !std::isfinite(low(c)))
                throw std::domain_error("pressure full flux gap: nonfinite cell");
            rhs2+=static_cast<long double>(rhs(c))*rhs(c);
        }
        if (!(rhs2>0.0L))
            throw std::invalid_argument("pressure full flux gap: zero RHS");
        PressureFullFluxGapReport result;
        for (std::size_t f=0;f<faces_.size();++f) {
            const auto& entry=faces_[f];
            if (!std::isfinite(base_face_flux[f]) ||
                !std::isfinite(full_face_flux[f]))
                throw std::domain_error("pressure full flux gap: nonfinite face");
            const long double po=static_cast<long double>(high(entry.owner))+low(entry.owner);
            const long double pn=entry.neighbour<0 ? 0.0L :
                static_cast<long double>(high(static_cast<std::size_t>(entry.neighbour)))+
                low(static_cast<std::size_t>(entry.neighbour));
            const long double implicit_flux=static_cast<long double>(base_face_flux[f])+
                static_cast<long double>(entry.coefficient)*(po-pn);
            const long double missing_flux=full_face_flux[f]-implicit_flux;
            if (!std::isfinite(implicit_flux) || !std::isfinite(missing_flux))
                throw std::domain_error("pressure full flux gap: nonfinite response");
            face_gap2+=missing_flux*missing_flux;
            if (std::abs(missing_flux)>result.maximum_abs_deferred_face) {
                result.maximum_abs_deferred_face=static_cast<double>(std::abs(missing_flux));
                result.worst_deferred_face=f;
            }
            implicit[entry.owner]+=implicit_flux;
            full[entry.owner]+=full_face_flux[f];
            missing[entry.owner]+=missing_flux;
            if (entry.neighbour>=0) {
                const auto n=static_cast<std::size_t>(entry.neighbour);
                implicit[n]-=implicit_flux;
                full[n]-=full_face_flux[f];
                missing[n]-=missing_flux;
            }
        }
        long double implicit2=0.0L,full2=0.0L,missing2=0.0L;
        for (std::size_t c=0;c<cells_;++c) {
            implicit2+=implicit[c]*implicit[c];
            full2+=full[c]*full[c];
            missing2+=missing[c]*missing[c];
            if (std::abs(full[c])>result.maximum_abs_cell_divergence) {
                result.maximum_abs_cell_divergence=static_cast<double>(std::abs(full[c]));
                result.worst_divergence_cell=c;
            }
            if (std::abs(missing[c])>result.maximum_abs_deferred_cell_divergence) {
                result.maximum_abs_deferred_cell_divergence=static_cast<double>(std::abs(missing[c]));
                result.worst_deferred_cell=c;
            }
        }
        result.implicit_continuity_relative=static_cast<double>(std::sqrt(implicit2/rhs2));
        result.full_continuity_relative=static_cast<double>(std::sqrt(full2/rhs2));
        result.deferred_divergence_relative=static_cast<double>(std::sqrt(missing2/rhs2));
        result.deferred_face_l2_relative=static_cast<double>(std::sqrt(face_gap2/rhs2));
        result.valid=std::isfinite(result.implicit_continuity_relative) &&
                     std::isfinite(result.full_continuity_relative) &&
                     std::isfinite(result.deferred_divergence_relative) &&
                     std::isfinite(result.deferred_face_l2_relative);
        if (!result.valid)
            throw std::domain_error("pressure full flux gap: invalid norm");
        return result;
    }

    // Returns b-A_face*(high+low), not a row-by-row rounded CSR residual.
    // Every internal face is written equal/opposite and every Dirichlet face
    // is written once to its owner. No fake boundary rows are constructed.
    std::pair<Vector, double> residual(const Vector& rhs,
                                       const Vector& high,
                                       const Vector& low) const {
        if (rhs.size() != cells_ || high.size() != cells_ ||
            low.size() != cells_ ||
            std::numeric_limits<long double>::digits <=
                std::numeric_limits<double>::digits)
            throw std::invalid_argument("pressure face operator: incompatible fields or precision");
        std::vector<long double> divergence(cells_, 0.0L);
        long double rhs_norm2 = 0.0L;
        for (std::size_t i = 0; i < cells_; ++i) {
            if (!std::isfinite(rhs(i)) || !std::isfinite(high(i)) ||
                !std::isfinite(low(i)))
                throw std::domain_error("pressure face operator: nonfinite field");
            rhs_norm2 += static_cast<long double>(rhs(i)) * rhs(i);
        }
        if (rhs_norm2 <= 0.0L)
            throw std::invalid_argument("pressure face operator: nonzero RHS required");
        for (const auto& f : faces_) {
            if (f.coefficient == 0.0) continue;
            const auto o = f.owner;
            const long double po = static_cast<long double>(high(o)) + low(o);
            const long double pn = f.neighbour < 0 ? 0.0L :
                static_cast<long double>(high(static_cast<std::size_t>(f.neighbour))) +
                low(static_cast<std::size_t>(f.neighbour));
            const long double flux = static_cast<long double>(f.coefficient) * (po - pn);
            divergence[o] += flux;
            if (f.neighbour >= 0)
                divergence[static_cast<std::size_t>(f.neighbour)] -= flux;
        }
        Vector result(cells_, 0.0);
        long double norm2 = 0.0L;
        for (std::size_t i = 0; i < cells_; ++i) {
            const long double r = static_cast<long double>(rhs(i)) - divergence[i];
            if (!std::isfinite(r) ||
                std::abs(r) > std::numeric_limits<double>::max())
                throw std::domain_error("pressure face operator: residual overflow");
            result(i) = static_cast<double>(r);
            norm2 += r * r;
        }
        return {std::move(result),
                static_cast<double>(std::sqrt(norm2 / rhs_norm2))};
    }

private:
    std::size_t cells_;
    std::vector<Face> faces_;
};

template<class IncrementSolver>
PressureFaceRefinementReport refine_pressure_on_implicit_faces(
    const ConservativePressureFaceOperator& faces, const SparseMatrix& csr,
    const Vector& rhs, Vector& high, Vector& low,
    double original_tolerance, std::size_t max_passes,
    IncrementSolver&& solve_increment) {
    if (faces.size() != csr.n_rows() || csr.n_rows() != csr.n_cols() ||
        rhs.size() != faces.size() || !(original_tolerance > 0.0) ||
        !std::isfinite(original_tolerance) || max_passes > 8)
        throw std::invalid_argument("pressure face refinement: invalid input");
    auto [r, relative] = faces.residual(rhs, high, low);
    PressureFaceRefinementReport report;
    report.valid = true;
    report.initial_face_relative = relative;
    report.final_face_relative = relative;
    for (std::size_t pass = 0; pass < max_passes &&
         relative > original_tolerance; ++pass) {
        Vector delta(csr.n_rows(), 0.0);
        if (!solve_increment(csr, r, delta)) break;
        Vector h = high, l = low;
        bool finite = true;
        for (std::size_t i = 0; i < csr.n_rows(); ++i)
            finite = compensated_add(h(i), l(i), delta(i)) && finite;
        if (!finite) break;
        auto [proposed, error] = faces.residual(rhs, h, l);
        if (!(error < relative)) break;
        high = std::move(h);
        low = std::move(l);
        r = std::move(proposed);
        relative = error;
        ++report.accepted_refinements;
        report.final_face_relative = relative;
    }
    report.implicit_face_converged = relative <= original_tolerance;
    report.final_csr_relative =
        compensated_pressure_residual(csr, rhs, high, low).second;
    report.original_csr_converged =
        report.final_csr_relative <= original_tolerance;
    // No acceptance based only on this implicit-face correction.
    return report;
}

} // namespace cfdx::core
