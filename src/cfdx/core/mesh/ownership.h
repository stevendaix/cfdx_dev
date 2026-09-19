// M0.1-T03 — Owner / Neighbour storage
//
// Spécification CFDX v0.7 §12 :
//   /mesh/faces/owner        uint32
//   /mesh/faces/neighbour    int32
//
// Une face interne possède :
//   owner >= 0, neighbour >= 0
// Une face frontière possède :
//   neighbour = -1
//
// Cette représentation reprend le principe général de polyMesh (§12).

#pragma once

#include <vector>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace cfdx {
namespace core {

class FaceOwnership {
public:
    using Owner = std::uint32_t;
    using Neighbour = std::int32_t;
    static constexpr Neighbour BOUNDARY = -1;

    FaceOwnership() = default;

    explicit FaceOwnership(std::size_t n_faces) {
        owner_.resize(n_faces, 0);
        neighbour_.resize(n_faces, BOUNDARY);
    }

    // --- Accès ---

    Owner owner(std::size_t i) const {
        check_index(i);
        return owner_[i];
    }

    Neighbour neighbour(std::size_t i) const {
        check_index(i);
        return neighbour_[i];
    }

    void set_owner(std::size_t i, Owner o) {
        check_index(i);
        owner_[i] = o;
    }

    void set_neighbour(std::size_t i, Neighbour n) {
        check_index(i);
        neighbour_[i] = n;
    }

    // --- Dimensions ---

    std::size_t size() const noexcept { return owner_.size(); }
    bool empty() const noexcept { return owner_.empty(); }

    void resize(std::size_t n) {
        owner_.resize(n, 0);
        neighbour_.resize(n, BOUNDARY);
    }

    void clear() {
        owner_.clear();
        neighbour_.clear();
    }

    // --- Accès bulk ---

    const Owner* owner_data() const noexcept { return owner_.data(); }
    const Neighbour* neighbour_data() const noexcept { return neighbour_.data(); }

    Owner* owner_data() noexcept { return owner_.data(); }
    Neighbour* neighbour_data() noexcept { return neighbour_.data(); }

    // --- Validation ---

    // Vérifie la cohérence owner/neighbour (§19) :
    //   - owner >= 0
    //   - si neighbour >= 0, owner != neighbour (pas de face auto-connectée)
    //   - si neighbour == -1, c'est une face de frontière
    bool is_consistent(std::size_t n_cells) const noexcept {
        for (std::size_t i = 0; i < size(); ++i) {
            if (owner_[i] >= static_cast<Owner>(n_cells)) return false;
            if (neighbour_[i] >= 0) {
                if (owner_[i] == static_cast<Owner>(neighbour_[i])) return false;
                if (static_cast<std::size_t>(neighbour_[i]) >= n_cells) return false;
            } else if (neighbour_[i] != BOUNDARY) {
                return false;  // valeur d'invalidité inconnue
            }
        }
        return true;
    }

    // Nombre de faces internes (avec voisin >= 0).
    std::size_t n_internal_faces() const noexcept {
        std::size_t count = 0;
        for (const auto n : neighbour_) {
            if (n >= 0) ++count;
        }
        return count;
    }

    // Nombre de faces de frontière (neighbour == -1).
    std::size_t n_boundary_faces() const noexcept {
        std::size_t count = 0;
        for (const auto n : neighbour_) {
            if (n == BOUNDARY) ++count;
        }
        return count;
    }

private:
    void check_index(std::size_t i) const {
        if (i >= size()) {
            throw std::out_of_range("FaceOwnership: index out of range");
        }
    }

    std::vector<Owner> owner_;
    std::vector<Neighbour> neighbour_;
};

}  // namespace core
}  // namespace cfdx