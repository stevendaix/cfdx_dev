// M0.1-T05 — Boundary patches
//
// Spécification CFDX v0.7 §14 :
//   /mesh/boundary/patches/
//
// Chaque patch contient :
//   name
//   type
//   face_ids
//   metadata
//
// Le concept de patch reste indépendant des conditions physiques (§14).

#pragma once

#include "index_types.h"
#include "ownership.h"
#include <vector>
#include <string>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <map>

namespace cfdx {
namespace core {

class Mesh;  // Forward declaration

enum class PatchType : std::uint8_t {
    WALL = 0,
    INLET,
    OUTLET,
    SYMMETRY,
    PERIODIC,
    INTERFACE,
    EMPTY,
    UNKNOWN
};

inline const char* to_string(PatchType t) {
    switch (t) {
        case PatchType::WALL:      return "wall";
        case PatchType::INLET:     return "inlet";
        case PatchType::OUTLET:    return "outlet";
        case PatchType::SYMMETRY:  return "symmetry";
        case PatchType::PERIODIC:  return "periodic";
        case PatchType::INTERFACE: return "interface";
        case PatchType::EMPTY:     return "empty";
        default:                   return "unknown";
    }
}

inline PatchType patch_type_from_string(const std::string& s) {
    if (s == "wall")      return PatchType::WALL;
    if (s == "inlet")     return PatchType::INLET;
    if (s == "outlet")    return PatchType::OUTLET;
    if (s == "symmetry")  return PatchType::SYMMETRY;
    if (s == "periodic")  return PatchType::PERIODIC;
    if (s == "interface") return PatchType::INTERFACE;
    if (s == "empty")     return PatchType::EMPTY;
    return PatchType::UNKNOWN;
}

struct Patch {
    std::string name;
    PatchType type = PatchType::UNKNOWN;
    std::vector<FaceIndex> face_ids;
    std::map<std::string, std::string> metadata;

    std::size_t size() const noexcept { return face_ids.size(); }
    bool empty() const noexcept { return face_ids.empty(); }
};

class BoundaryPatches {
public:
    BoundaryPatches() = default;

    // --- Accès ---

    std::size_t n_patches() const noexcept { return patches_.size(); }

    const Patch& patch(std::size_t i) const {
        check_index(i);
        return patches_[i];
    }

    Patch& patch(std::size_t i) {
        check_index(i);
        return patches_[i];
    }

    // Trouve un patch par son nom. Retourne n_patches() si absent.
    std::size_t find(const std::string& name) const {
        for (std::size_t i = 0; i < patches_.size(); ++i) {
            if (patches_[i].name == name) return i;
        }
        return patches_.size();
    }

    bool has_patch(const std::string& name) const {
        return find(name) < patches_.size();
    }

    // --- Modification ---

    std::size_t add_patch(const Patch& p) {
        if (has_patch(p.name)) {
            throw std::runtime_error("BoundaryPatches: duplicate patch name '" + p.name + "'");
        }
        patches_.push_back(p);
        return patches_.size() - 1;
    }

    std::size_t add_patch(const std::string& name, PatchType type) {
        Patch p;
        p.name = name;
        p.type = type;
        return add_patch(p);
    }

    void remove_patch(std::size_t i) {
        check_index(i);
        patches_.erase(patches_.begin() + static_cast<std::ptrdiff_t>(i));
    }

    void clear() { patches_.clear(); }

    // --- Validation ---

    // Vérifie que les face_ids sont valides et non chevauchées.
    bool is_consistent(std::size_t n_faces) const {
        std::vector<bool> seen(n_faces, false);
        for (const auto& p : patches_) {
            for (const auto fid : p.face_ids) {
                if (fid >= static_cast<FaceIndex>(n_faces)) return false;
                if (seen[fid]) return false;  // face dans deux patches
                seen[fid] = true;
            }
        }
        return true;
    }

    // Vérifie la cohérence complète des patches avec la topologie du maillage (§19).
    // Invariant : union(all patch faces) == all faces where neighbour == -1
    //             intersection(patches) = ∅
    //             internal faces (neighbour >= 0) must not be in any patch
    bool is_consistent_with_mesh(const Mesh& m) const;

    // Nombre total de faces de frontière (somme des patches).
    std::size_t n_boundary_faces() const {
        std::size_t total = 0;
        for (const auto& p : patches_) total += p.size();
        return total;
    }

    // --- Accès bulk ---

    const std::vector<Patch>& patches() const noexcept { return patches_; }

private:
    void check_index(std::size_t i) const {
        if (i >= n_patches()) {
            throw std::out_of_range("BoundaryPatches: patch index out of range");
        }
    }

    std::vector<Patch> patches_;
};

}  // namespace core
}  // namespace cfdx
