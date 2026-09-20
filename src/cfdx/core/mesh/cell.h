// M0.1-T04 — Cell connectivity (CSR)
//
// Spécification CFDX v0.7 §13 :
//   /mesh/cells/faces
//   /mesh/cells/faces_offsets
//
// Format CSR. Une cellule peut posséder un nombre quelconque de faces.

#pragma once

#include "index_types.h"
#include <vector>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace cfdx {
namespace core {

class CellConnectivity {
public:
    using FaceId = FaceIndex;
    using Offset = ::cfdx::core::Offset;

    CellConnectivity() = default;

    void reserve(std::size_t n_cells, std::size_t n_faces) {
        faces_.reserve(n_faces);
        offsets_.reserve(n_cells + 1);
    }

    void push_cell(const std::vector<FaceIndex>& cell_faces) {
        if (offsets_.empty()) {
            offsets_.push_back(0);
        }
        for (const auto f : cell_faces) {
            faces_.push_back(f);
        }
        offsets_.push_back(static_cast<Offset>(faces_.size()));
    }

    std::size_t n_cells() const noexcept {
        return offsets_.empty() ? 0 : offsets_.size() - 1;
    }

    std::size_t n_face_refs() const noexcept { return faces_.size(); }

    std::size_t cell_size(std::size_t i) const {
        check_cell(i);
        return offsets_[i + 1] - offsets_[i];
    }

    Offset cell_offset(std::size_t i) const {
        check_cell(i);
        return offsets_[i];
    }

    const FaceIndex* faces_data() const noexcept { return faces_.data(); }
    const Offset* offsets_data() const noexcept { return offsets_.data(); }

    FaceIndex* faces_data() noexcept { return faces_.data(); }
    Offset* offsets_data() noexcept { return offsets_.data(); }

    const std::vector<FaceIndex>& faces() const noexcept { return faces_; }
    const std::vector<Offset>& offsets() const noexcept { return offsets_; }

    // Vérifie la cohérence interne.
    bool is_consistent() const noexcept {
        if (offsets_.empty() && faces_.empty()) return true;
        if (offsets_.empty() || faces_.empty()) return false;
        if (offsets_.size() < 2) return false;
        if (offsets_[0] != 0) return false;
        if (offsets_.back() != static_cast<Offset>(faces_.size())) return false;
        for (std::size_t i = 1; i < offsets_.size(); ++i) {
            if (offsets_[i] <= offsets_[i - 1]) return false;
        }
        return true;
    }

    // Vérifie que chaque face reference est dans [0, n_faces).
    bool face_ids_valid(std::size_t n_faces) const noexcept {
        for (const auto f : faces_) {
            if (f >= static_cast<FaceIndex>(n_faces)) return false;
        }
        return true;
    }

    void clear() {
        faces_.clear();
        offsets_.clear();
    }

private:
    void check_cell(std::size_t i) const {
        if (i >= n_cells()) {
            throw std::out_of_range("CellConnectivity: cell index out of range");
        }
    }

    std::vector<FaceIndex> faces_;
    std::vector<Offset> offsets_;
};

}  // namespace core
}  // namespace cfdx
