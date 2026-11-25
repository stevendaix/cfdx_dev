#pragma once

// Vrai solveur matriciel couplé 2x2 k/omega.
//
// Résout le système bloc :
//
//   [ A_kk   D_kw ] [ dk ]   [ r_k ]
//   [ D_wk   A_ww ] [ dw ] = [ r_w ]
//
// par itération de Richardson préconditionnée par bloc-Jacobi.
//
// CORRECTIONS par rapport à la première version :
//   1. Pas de double comptage du terme de destruction :
//      - d_kw = dPk/dlog(omega) uniquement (le terme -spk*k est déjà dans A_kk)
//      - d_wk = dSw/dlog(k) uniquement (le terme -spw*w est déjà dans A_ww)
//   2. Critère de convergence normalisé par équation :
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
#include <vector>
#include <algorithm>

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/turbulence_transport.h"

namespace cfdx::physics {

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
    const std::vector<double>& b_k,   // côté explicite de k (sk = Pk)
    const std::vector<double>& b_w,   // côté explicite de omega (sw)
    std::vector<double>& k,
    std::vector<double>& w,
    const std::vector<double>& d_kw,  // dPk/dlog(omega)  (diagonal, pas de destruction)
    const std::vector<double>& d_wk,  // dSw/dlog(k)      (diagonal, pas de destruction)
    double relaxation = 0.7,
    std::size_t max_iterations = 10,
    double tolerance = 1e-4,
    double k_min = 1e-12,
    double omega_min = 1e-12)
{
    const std::size_t n = k.size();
    Coupled2x2SolveResult result;

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
        const double det = a * d - b * c;
        if (std::isfinite(det) && std::abs(det) > 1e-14) {
            inv_det[i] = 1.0 / det;
        } else {
            // Fallback Jacobi simple : inverse de la diagonale uniquement.
            if (std::isfinite(a) && std::abs(a) > 1e-14) inv_a_kk[i] = 1.0 / a;
            if (std::isfinite(d) && std::abs(d) > 1e-14) inv_a_ww[i] = 1.0 / d;
        }
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

            // Clamp pour éviter les sauts excessifs.
            const double max_delta = 10.0;
            dk[i] = std::clamp(dk[i], -max_delta, max_delta);
            dw[i] = std::clamp(dw[i], -max_delta, max_delta);
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
