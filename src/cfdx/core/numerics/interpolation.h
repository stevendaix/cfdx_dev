// M0.6-T01 — Cell→Face interpolation
//
// Spécification CFDX v0.7 §29 :
//   L'interpolation cellule → face supporte :
//     linear
//     upwind
//     limited
//
//   Elle est générique et indépendante de la physique.
//
// Schémas (pour un champ de dimension D, composante par composante) :
//   - linear   : 0.5 * (owner + neighbour) pour les faces internes ;
//                owner pour les faces de frontière.
//   - upwind   : owner pour toutes les faces (suppose le flux positif).
//   - limited  : linéaire borné par [min(owner,neighbour), max(owner,neighbour)].
//
// Le résultat est un Field<double, FACE> de même dimension que le champ source.

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/ownership.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <string>
#include <algorithm>
#include <stdexcept>

namespace cfdx {
namespace core {

enum class InterpScheme : std::uint8_t {
    LINEAR = 0,
    UPWIND,
    LIMITED
};

enum class LimiterType : std::uint8_t {
    NONE = 0,
    MINMOD,
    VANLEER,
    SUPERBEE,
    VAN_ALBADA
};

inline const char* to_string(InterpScheme s) {
    switch (s) {
        case InterpScheme::LINEAR:  return "linear";
        case InterpScheme::UPWIND:  return "upwind";
        case InterpScheme::LIMITED: return "limited";
        default:                    return "unknown";
    }
}

inline const char* to_string(LimiterType l) {
    switch (l) {
        case LimiterType::NONE:      return "none";
        case LimiterType::MINMOD:    return "minmod";
        case LimiterType::VANLEER:   return "vanleer";
        case LimiterType::SUPERBEE:  return "superbee";
        case LimiterType::VAN_ALBADA: return "vanalbada";
        default:                     return "unknown";
    }
}

inline InterpScheme interp_scheme_from_string(const std::string& s) {
    if (s == "linear")  return InterpScheme::LINEAR;
    if (s == "upwind")  return InterpScheme::UPWIND;
    if (s == "limited") return InterpScheme::LIMITED;
    throw std::runtime_error("interp_scheme_from_string: unknown scheme '" + s + "'");
}

inline LimiterType limiter_type_from_string(const std::string& s) {
    if (s == "none")      return LimiterType::NONE;
    if (s == "minmod")    return LimiterType::MINMOD;
    if (s == "vanleer")   return LimiterType::VANLEER;
    if (s == "superbee")  return LimiterType::SUPERBEE;
    if (s == "vanalbada") return LimiterType::VAN_ALBADA;
    throw std::runtime_error("limiter_type_from_string: unknown limiter '" + s + "'");
}

// Applique un limiteur TVD réel (M0.6-T02 correction).
// v_owner   : valeur dans la cellule propriétaire
// v_upwind  : valeur dans la cellule amont (upwind)
// v_extrap  : valeur extrapolée linéairement au centre de la face (v_owner + grad · d)
// limiter    : MINMOD, VANLEER, SUPERBEE, VAN_ALBADA
inline double apply_limiter_tvd(double v_owner, double v_upwind, double v_extrap, LimiterType limiter) {
    double denominator = v_extrap - v_owner;
    double r = 1.0;    
    if (std::abs(denominator) > 1e-15) {
        r = (v_upwind - v_owner) / denominator;
    }

    double psi = 1.0;
    switch (limiter) {
        case LimiterType::MINMOD:
            psi = std::max(0.0, std::min(1.0, r));
            break;
        case LimiterType::VANLEER:
            psi = (r + std::abs(r)) / (1.0 + std::abs(r) + 1e-15);
            break;
        case LimiterType::SUPERBEE:
            psi = std::max(0.0, std::max(std::min(1.0, 2.0 * r), std::min(2.0, r)));
            break;
        case LimiterType::VAN_ALBADA:
            psi = (r * r + r) / (r * r + 1.0 + 1e-15);
            break;
        default:
            psi = 1.0;
    }
    
    // Valeur limitée
    double v_limited = v_owner + psi * (v_extrap - v_owner);
    
    // Bornage strict (garantie TVD) : la valeur à la face ne peut pas dépasser les extrêmes locaux
    double v_min = std::min(v_owner, v_upwind);
    double v_max = std::max(v_owner, v_upwind);
    
    return std::max(v_min, std::min(v_max, v_limited));
}

// Interpole un champ cellulaire (Field<double, CELL>) vers un champ de faces.
// Le champ source peut être scalaire (dim=1) ou vecteur (dim=3) — le résultat
// hérite de la même dimension.
//
// Pour UPWIND : nécessite un champ de flux (phi_f = U_f · Sf_f) pour déterminer le signe.
// Si flux >= 0 : valeur owner ; si flux < 0 : valeur neighbour.
inline Field<double, Location::FACE> interpolate_cell_to_face(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh,
    const GeometryCache& geometry,
    InterpScheme scheme,
    const Field<double, Location::FACE>* face_flux = nullptr,
    LimiterType limiter_type = LimiterType::NONE,
    const Field<double, Location::CELL>* cell_gradient = nullptr)
{
    const std::size_t n_faces = mesh.n_faces();
    const std::size_t n_cells = mesh.n_cells();
    const std::size_t dim = cell_field.dimension();

    if (!is_valid(geometry, mesh))
        throw std::invalid_argument("interpolate_cell_to_face: invalid geometry cache");
    if (cell_field.size() != n_cells)
        throw std::runtime_error("interpolate_cell_to_face: field size != n_cells");
    if (dim == 0)
        throw std::runtime_error("interpolate_cell_to_face: dimension must be >= 1");
    if (face_flux && (face_flux->size() != n_faces || face_flux->dimension() != 1))
        throw std::runtime_error("interpolate_cell_to_face: face_flux must be scalar with n_faces values");
    if (scheme == InterpScheme::UPWIND && !face_flux)
        throw std::runtime_error("interpolate_cell_to_face: UPWIND requires face_flux");
    if (scheme == InterpScheme::LIMITED) {
        if (!face_flux)
            throw std::runtime_error("interpolate_cell_to_face: LIMITED requires face_flux");
        if (!cell_gradient)
            throw std::runtime_error("interpolate_cell_to_face: LIMITED requires cell_gradient");
        if (cell_field.dimension() != 1 || cell_gradient->dimension() != 3 ||
            cell_gradient->size() != n_cells)
            throw std::runtime_error(
                "interpolate_cell_to_face: LIMITED currently supports scalar fields with a 3-component gradient");
    }

    Field<double, Location::FACE> result(
        n_faces, cell_field.name(), cell_field.metadata().unit, dim);

    const FaceOwnership& own = mesh.ownership();
    const double* flux = face_flux ? face_flux->component_data(0) : nullptr;

    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t owner = own.owner(f);
        if (owner >= n_cells)
            throw std::runtime_error("interpolate_cell_to_face: owner index out of range");

        const std::int64_t nb_i = own.neighbour(f);
        const bool internal = nb_i >= 0;
        if (internal && static_cast<std::size_t>(nb_i) >= n_cells)
            throw std::runtime_error("interpolate_cell_to_face: neighbour index out of range");
        const std::size_t neighbour = internal ? static_cast<std::size_t>(nb_i) : owner;

        for (std::size_t d = 0; d < dim; ++d) {
            const double vo = cell_field.component_data(d)[owner];
            const double vn = cell_field.component_data(d)[neighbour];
            double v = vo;

            if (!internal) {
                v = vo;
            } else if (scheme == InterpScheme::LINEAR) {
                v = 0.5 * (vo + vn);
            } else if (scheme == InterpScheme::UPWIND) {
                v = (flux[f] >= 0.0) ? vo : vn;
            } else if (scheme == InterpScheme::LIMITED) {
                // Reconstruct from the upwind cell to the actual face centre.
                // The limiter is applied to the directional ratio along the
                // owner-neighbour line; the final value is strictly bounded by
                // the two adjacent cell values.
                const bool owner_upwind = flux[f] >= 0.0;
                const std::size_t up = owner_upwind ? owner : neighbour;
                const std::size_t down = owner_upwind ? neighbour : owner;
                const double v_up = cell_field.component_data(0)[up];
                const double v_down = cell_field.component_data(0)[down];

                const double* gx = cell_gradient->component_data(0);
                const double* gy = cell_gradient->component_data(1);
                const double* gz = cell_gradient->component_data(2);
                const Vec3 dface = geometry.face_centres[f] - geometry.cell_centres[up];
                const double delta_extrap =
                    gx[up] * dface.x + gy[up] * dface.y + gz[up] * dface.z;

                const double delta_neighbour = v_down - v_up;
                const double r = (std::abs(delta_extrap) > 1e-14)
                    ? delta_neighbour / delta_extrap
                    : 0.0;

                double psi = 1.0;
                switch (limiter_type) {
                    case LimiterType::NONE:
                        psi = 1.0;
                        break;
                    case LimiterType::MINMOD:
                        psi = std::max(0.0, std::min(1.0, r));
                        break;
                    case LimiterType::VANLEER:
                        psi = (r + std::abs(r)) / (1.0 + std::abs(r) + 1e-15);
                        break;
                    case LimiterType::SUPERBEE:
                        psi = std::max(0.0,
                            std::max(std::min(1.0, 2.0 * r),
                                     std::min(2.0, r)));
                        break;
                    case LimiterType::VAN_ALBADA:
                        psi = (r * r + r) / (r * r + 1.0 + 1e-15);
                        break;
                }

                const double reconstructed = v_up + psi * delta_extrap;
                v = std::max(std::min(v_up, v_down),
                             std::min(std::max(v_up, v_down), reconstructed));
            }
            result.component_data(d)[f] = v;
        }
    }
    return result;
}

inline Field<double, Location::FACE> interpolate_cell_to_face(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh,
    InterpScheme scheme,
    const Field<double, Location::FACE>* face_flux = nullptr,
    LimiterType limiter_type = LimiterType::NONE,
    const Field<double, Location::CELL>* cell_gradient = nullptr)
{
    const GeometryCache geometry = make_geometry_cache(mesh);
    return interpolate_cell_to_face(
        cell_field, mesh, geometry, scheme, face_flux, limiter_type, cell_gradient);
}

}  // namespace core
}  // namespace cfdx
