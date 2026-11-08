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
        return std::count_if(faces_.begin(), faces_.end(),
                             [](const Face& f) { return f.coefficient > 0.0; });
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
