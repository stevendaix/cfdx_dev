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
#include <cstddef>
#include <algorithm>
#include <stdexcept>

namespace cfdx {
namespace core {

enum class InterpScheme : std::uint8_t {
    LINEAR = 0,
    UPWIND,
    LIMITED
};

inline const char* to_string(InterpScheme s) {
    switch (s) {
        case InterpScheme::LINEAR:  return "linear";
        case InterpScheme::UPWIND:  return "upwind";
        case InterpScheme::LIMITED: return "limited";
        default:                    return "unknown";
    }
}

inline InterpScheme interp_scheme_from_string(const std::string& s) {
    if (s == "linear")  return InterpScheme::LINEAR;
    if (s == "upwind")  return InterpScheme::UPWIND;
    if (s == "limited") return InterpScheme::LIMITED;
    throw std::runtime_error("interp_scheme_from_string: unknown scheme '" + s + "'");
}

// Interpole un champ cellulaire (Field<double, CELL>) vers un champ de faces.
// Le champ source peut être scalaire (dim=1) ou vecteur (dim=3) — le résultat
// hérite de la même dimension.
inline Field<double, Location::FACE> interpolate_cell_to_face(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh,
    InterpScheme scheme)
{
    const std::size_t n_faces = mesh.n_faces();
    const std::size_t dim = cell_field.dimension();
    if (dim == 0) throw std::runtime_error("interpolate_cell_to_face: dimension must be >= 1");

    Field<double, Location::FACE> result(
        n_faces, cell_field.name(), cell_field.metadata().unit, dim);

    const FaceOwnership& own = mesh.ownership();

    const double* src = cell_field.data();
    double* dst = result.data();

    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t owner = own.owner(f);
        if (owner >= mesh.n_cells()) {
            throw std::runtime_error("interpolate_cell_to_face: owner index out of range");
        }

        const bool is_internal = (own.neighbour(f) >= 0);
        const std::size_t neighbour = is_internal
            ? static_cast<std::size_t>(own.neighbour(f))
            : owner;

        if (is_internal && neighbour >= mesh.n_cells()) {
            throw std::runtime_error("interpolate_cell_to_face: neighbour index out of range");
        }

        for (std::size_t d = 0; d < dim; ++d) {
            const double v_owner = src[owner * dim + d];
            const double v_neigh = src[neighbour * dim + d];
            double v;

            switch (scheme) {
                case InterpScheme::UPWIND:
                    v = v_owner;
                    break;
                case InterpScheme::LINEAR:
                    v = 0.5 * (v_owner + v_neigh);
                    break;
                case InterpScheme::LIMITED: {
                    const double lo = std::min(v_owner, v_neigh);
                    const double hi = std::max(v_owner, v_neigh);
                    const double linear = 0.5 * (v_owner + v_neigh);
                    v = std::max(lo, std::min(hi, linear));
                    break;
                }
                default:
                    throw std::runtime_error("interpolate_cell_to_face: unknown scheme");
            }
            dst[f * dim + d] = v;
        }
    }

    return result;
}

}  // namespace core
}  // namespace cfdx