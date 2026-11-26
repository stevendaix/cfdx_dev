#pragma once

// Vrai solveur matriciel couplé 2x2 k/omega.
//
// Résout le système bloc :
//
//   [ A_kk   C_kw ] [ k ]   [ b_k + C_kw*w0 ]
//   [ C_wk   A_ww ] [ w ] = [ b_w + C_wk*k0 ]
//
// C_kw=-d(S_k)/d(omega) and C_wk=-d(S_omega)/d(k), in derivatives
// with respect to the dimensional fields (not logarithmic variables).
// The affine RHS shift makes the tangent pass through the current iterate,
// preserving the original segregated FV residual at the linearization point.
// The assembled ScalarEquation RHS is used so boundary and relaxation terms
// are not silently dropped.
//
// par itération de Richardson préconditionnée par bloc-Jacobi.
//
// Contrat de linéarisation : les coefficients croisés sont les dérivées
// négatives des sources explicites par rapport aux champs dimensionnels.
// Les contributions de destruction implicites restent dans A_kk/A_ww.
// Le décalage affine du RHS est indispensable : sans lui, l'ajout des blocs
// croisés change le point fixe des deux équations scalaires.
// Critère de convergence normalisé par équation :
//      R_k* = ||r_k|| / (||b_k|| + eps), R_w* = ||r_w|| / (||b_w|| + eps)
//   3. Pas de coupure de positivité forcée à l'intérieur du solveur linéaire
//      (la positivité est gérée par le code appelant via enforce_turbulence_bounds).
//   4. Fallback si le déterminant d'un bloc 2x2 est trop petit : on utilise
//      uniquement la diagonale (Jacobi simple) pour cette cellule.
//
// Aucun changement à finite_volume_transport.h n'est nécessaire : A_kk et
// A_ww sont assemblés par assemble_scalar_equation (existant), et leurs
// matvec sont utilisés directement.

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <algorithm>

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/turbulence_transport.h"

namespace cfdx::physics {

// Convert a relative finite-difference perturbation of a source density
// into the negative derivative coefficient used by the assembled FV matrix.
// The cell volume is essential: ScalarEquation stores integrated source terms.
inline double sst_negative_source_derivative_coefficient(
    double source_at_base, double source_at_perturbed,
    double relative_step, double base_variable, double cell_volume)
{
    if (!std::isfinite(source_at_base) || !std::isfinite(source_at_perturbed) ||
        !(relative_step > 0.0) || !std::isfinite(relative_step) ||
        !(base_variable > 0.0) || !std::isfinite(base_variable) ||
        !(cell_volume > 0.0) || !std::isfinite(cell_volume))
        throw std::invalid_argument("SST coupled 2x2: invalid source derivative data");
    const double coefficient =
        -(source_at_perturbed - source_at_base) /
        (relative_step * base_variable) * cell_volume;
    if (!std::isfinite(coefficient))
        throw std::domain_error("SST coupled 2x2: source derivative overflow");
    return coefficient;
}

// Résultat du solveur couplé 2x2 matriciel.
struct Coupled2x2SolveResult {
    std::size_t iterations = 0;
    double res_k_norm = 0.0;   // résidu normalisé de k
    double res_w_norm = 0.0;   // résidu normalisé de omega
    bool converged = false;
};

// Résout le système couplé 2x2 par Richardson + bloc-Jacobi.
inline Coupled2x2SolveResult solve_coupled_2x2_matrix(
    const ScalarEquation& eqk,
    const ScalarEquation& eqw,
    std::vector<double>& k,
    std::vector<double>& w,
    const std::vector<double>& d_kw,  // -d(S_k)/d(omega), dimensional cross coefficient
    const std::vector<double>& d_wk,  // -d(S_omega)/d(k), dimensional cross coefficient
    double relaxation = 0.7,
    std::size_t max_iterations = 10,
    double tolerance = 1e-4,
    double k_min = 1e-12,
    double omega_min = 1e-12)
{
    const std::size_t n = k.size();
    Coupled2x2SolveResult result;
    if (n == 0 || w.size() != n || d_kw.size() != n || d_wk.size() != n ||
        eqk.rhs.size() != n || eqw.rhs.size() != n ||
        eqk.matrix.n_rows() != n || eqk.matrix.n_cols() != n ||
        eqw.matrix.n_rows() != n || eqw.matrix.n_cols() != n ||
        eqk.diagonal.size() != n || eqw.diagonal.size() != n ||
        !(relaxation > 0.0 && relaxation <= 1.0) ||
        !std::isfinite(relaxation) || max_iterations == 0 ||
        !(tolerance > 0.0) || !std::isfinite(tolerance) ||
        !(k_min > 0.0) || !(omega_min > 0.0))
        throw std::invalid_argument("SST coupled 2x2: invalid system or controls");

    // Start from the complete assembled scalar RHS, including boundary
    // contributions. The affine tangent shift preserves the base residual:
    // b_tangent - A*x0 - C*y0 == b_FV - A*x0.
    std::vector<double> b_k(n), b_w(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(eqk.rhs(i)) || !std::isfinite(eqw.rhs(i)) ||
            !std::isfinite(k[i]) || !std::isfinite(w[i]) ||
            !std::isfinite(d_kw[i]) || !std::isfinite(d_wk[i]) ||
            !(k[i] > 0.0) || !(w[i] > 0.0))
            throw std::invalid_argument("SST coupled 2x2: non-finite or non-positive state");
        b_k[i] = eqk.rhs(i) + d_kw[i] * w[i];
        b_w[i] = eqw.rhs(i) + d_wk[i] * k[i];
        if (!std::isfinite(b_k[i]) || !std::isfinite(b_w[i]))
            throw std::invalid_argument("SST coupled 2x2: affine RHS overflow");
    }

    const auto& diag_kk = eqk.diagonal;
    const auto& diag_ww = eqw.diagonal;

    // Préconditionneur bloc-Jacobi : inverse [diag_kk  d_kw; d_wk  diag_ww]
    // par cellule. Si le déterminant est trop petit, on utilise uniquement
    // la diagonale (Jacobi simple) pour cette cellule.
    std::vector<double> inv_det(n, 0.0);
    std::vector<double> inv_a_kk(n, 0.0);
    std::vector<double> inv_a_ww(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double a = diag_kk[i];
        const double d = diag_ww[i];
        const double b = d_kw[i];
        const double c = d_wk[i];
        if (!std::isfinite(a) || !std::isfinite(d) ||
            !(a > 0.0) || !(d > 0.0))
            throw std::invalid_argument("SST coupled 2x2: non-positive scalar diagonal");
        const long double ad = static_cast<long double>(a) * d;
        const long double bc = static_cast<long double>(b) * c;
        const long double det = ad - bc;
        const long double det_scale = std::max(std::abs(ad),std::abs(bc));
        const long double det_tolerance =
            64.0L * std::numeric_limits<double>::epsilon() * det_scale;
        if (std::isfinite(det) && std::abs(det) > det_tolerance) {
            const double inverse = static_cast<double>(1.0L / det);
            if (std::isfinite(inverse)) {
                inv_det[i] = inverse;
                continue;
            }
        }
        // Scale-independent fallback: both scalar diagonals have already
        // been validated positive, so no dimensioned absolute cutoff is used.
        inv_a_kk[i] = 1.0 / a;
        inv_a_ww[i] = 1.0 / d;
        if (!std::isfinite(inv_a_kk[i]) || !std::isfinite(inv_a_ww[i]))
            throw std::domain_error("SST coupled 2x2: scalar diagonal inverse overflow");
    }

    // Normes des côtés explicites pour la normalisation des résidus.
    double norm_b_k = 0.0, norm_b_w = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        norm_b_k = std::max(norm_b_k, std::abs(b_k[i]));
        norm_b_w = std::max(norm_b_w, std::abs(b_w[i]));
    }
    norm_b_k = std::max(norm_b_k, 1e-12);
    norm_b_w = std::max(norm_b_w, 1e-12);

    // Vecteurs de travail.
    std::vector<double> res_k(n), res_w(n);
    std::vector<double> dk(n), dw(n);
    double final_res_k = 0.0, final_res_w = 0.0;

    for (std::size_t it = 0; it < max_iterations; ++it) {
        // Matvec : A_kk * k et A_ww * w.
        const auto akk = eqk.matrix.matvec(k);
        const auto aww = eqw.matrix.matvec(w);

        // Résidus : r = b - A*x - cross*other.
        // Attention : A_kk contient déjà spk (destruction de k) sur sa diagonale,
        // A_ww contient déjà spw (destruction de omega) sur sa diagonale.
        // d_kw et d_wk contiennent uniquement les termes croisés hors-diagonale.
        double max_res_k = 0.0, max_res_w = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            res_k[i] = b_k[i] - akk[i] - d_kw[i] * w[i];
            res_w[i] = b_w[i] - d_wk[i] * k[i] - aww[i];
            max_res_k = std::max(max_res_k, std::abs(res_k[i]));
            max_res_w = std::max(max_res_w, std::abs(res_w[i]));
        }

        const double res_k_norm = max_res_k / norm_b_k;
        const double res_w_norm = max_res_w / norm_b_w;
        final_res_k = res_k_norm;
        final_res_w = res_w_norm;

        if (res_k_norm < tolerance && res_w_norm < tolerance) {
            result.iterations = it;
            result.res_k_norm = res_k_norm;
            result.res_w_norm = res_w_norm;
            result.converged = true;
            return result;
        }

        // Préconditionneur bloc-Jacobi : M^{-1} * r.
        for (std::size_t i = 0; i < n; ++i) {
            const double idet = inv_det[i];
            const double rk = res_k[i];
            const double rw = res_w[i];

            if (idet != 0.0) {
                // Inverse complet du bloc 2x2.
                const double a = diag_kk[i];
                const double d = diag_ww[i];
                const double b = d_kw[i];
                const double c = d_wk[i];
                dk[i] = idet * (d * rk - b * rw);
                dw[i] = idet * (-c * rk + a * rw);
            } else {
                // Jacobi simple (diagonale uniquement).
                dk[i] = inv_a_kk[i] * rk;
                dw[i] = inv_a_ww[i] * rw;
            }

            // Bound corrections in dimensionless relative units; k and
            // omega have different physical dimensions and cannot share an
            // absolute increment cap.
            const double max_delta_k = 10.0 * std::max(k[i], k_min);
            const double max_delta_w = 10.0 * std::max(w[i], omega_min);
            dk[i] = std::clamp(dk[i], -max_delta_k, max_delta_k);
            dw[i] = std::clamp(dw[i], -max_delta_w, max_delta_w);
        }

        // Mise à jour Richardson avec sous-relaxation + floor de positivité.
        // Le floor empêche k/omega de devenir négatifs durant les itérations
        // intérieures (ce qui casserait le matvec suivant via log/racine dans
        // les sources SST).
        for (std::size_t i = 0; i < n; ++i) {
            k[i] += relaxation * dk[i];
            w[i] += relaxation * dw[i];
            if (!(k[i] > k_min)) k[i] = k_min;
            if (!(w[i] > omega_min)) w[i] = omega_min;
        }
    }

    result.iterations = max_iterations;
    result.res_k_norm = final_res_k;
    result.res_w_norm = final_res_w;
    result.converged = false;
    return result;
}

} // namespace cfdx::physics
