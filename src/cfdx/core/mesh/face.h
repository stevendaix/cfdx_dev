// M0.1-T02 — Face CSR connectivity
//
// Spécification CFDX v0.7 §11 :
//   /mesh/faces/vertices
//   /mesh/faces/vertices_offsets
//
// Représentation CSR (Compressed Sparse Row) :
//   vertices[]       : flat list of vertex indices, one face after another
//   vertices_offsets[i] = index in vertices[] where face i starts
//   vertices_offsets[n_faces] = total length of vertices[]
//
// Permet de représenter des faces de nombre de sommets arbitraire
// sans allocation par face (§11.1).

#pragma once

#include "index_types.h"
#include <vector>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace cfdx {
namespace core {

class FaceConnectivity {
public:
    using Index = FaceIndex;
    using Offset = ::cfdx::core::Offset;

    FaceConnectivity() = default;

    // --- Construction progressive ---

    void reserve(std::size_t n_faces, std::size_t n_vertices) {
        vertices_.reserve(n_vertices);
        offsets_.reserve(n_faces + 1);
    }

    // Ajoute une face à partir d'une liste d'indices de sommets.
    void push_face(const std::vector<Index>& face_vertices) {
        if (offsets_.empty()) {
            offsets_.push_back(0);
        }
        for (const auto v : face_vertices) {
            vertices_.push_back(v);
        }
        offsets_.push_back(static_cast<Offset>(vertices_.size()));
    }

    // --- Accès ---

    // Nombre de faces.
    std::size_t n_faces() const noexcept {
        return offsets_.empty() ? 0 : offsets_.size() - 1;
    }

    // Nombre total de sommets (arêtes) references.
    std::size_t n_vertices() const noexcept { return vertices_.size(); }

    // Nombre de sommets de la face i.
    std::size_t face_size(std::size_t i) const {
        check_face(i);
        return offsets_[i + 1] - offsets_[i];
    }

    // Début de la face i dans le tableau vertices.
    Offset face_offset(std::size_t i) const {
        check_face(i);
        return offsets_[i];
    }

    // Accès direct aux données brutes.
    const Index* vertices_data() const noexcept { return vertices_.data(); }
    const Offset* offsets_data() const noexcept { return offsets_.data(); }

    Index* vertices_data() noexcept { return vertices_.data(); }
    Offset* offsets_data() noexcept { return offsets_.data(); }

    const std::vector<Index>& vertices() const noexcept { return vertices_; }
    const std::vector<Offset>& offsets() const noexcept { return offsets_; }

    // --- Validation ---

    // Vérifie la cohérence interne des offsets.
    // offsets doit être strictement croissant et offsets.back() == vertices.size().
    bool is_consistent() const noexcept {
        if (offsets_.empty() && vertices_.empty()) return true;
        if (offsets_.empty() || vertices_.empty()) return false;
        if (offsets_.size() < 2) return false;
        if (offsets_[0] != 0) return false;
        if (offsets_.back() != static_cast<Offset>(vertices_.size())) return false;
        for (std::size_t i = 1; i < offsets_.size(); ++i) {
            if (offsets_[i] <= offsets_[i - 1]) return false;
        }
        return true;
    }

    // Vérifie que chaque index de sommet est dans [0, n_points).
    bool indices_valid(std::size_t n_points) const noexcept {
        for (const auto v : vertices_) {
            if (v >= static_cast<Index>(n_points)) return false;
        }
        return true;
    }

    void clear() {
        vertices_.clear();
        offsets_.clear();
    }

    void build_from_scratch(const std::vector<std::vector<Index>>& all_faces_vertices) {
        offsets_.clear();
        vertices_.clear();
        offsets_.reserve(all_faces_vertices.size() + 1);
        offsets_.push_back(0);
        for (const auto& face : all_faces_vertices) {
            for (auto v : face) {
                vertices_.push_back(v);
            }
            offsets_.push_back(static_cast<Offset>(vertices_.size()));
        }
    }

private:
    void check_face(std::size_t i) const {
        if (i >= n_faces()) {
            throw std::out_of_range("FaceConnectivity: face index out of range");
        }
    }

    std::vector<FaceIndex> vertices_;
    std::vector<Offset> offsets_;
};

}  // namespace core
}  // namespace cfdx
