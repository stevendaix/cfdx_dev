// M0.5-T01 — BoundaryField / PatchField
//
// Spécification CFDX v0.7 §26 :
//   Le champ volumique et sa représentation frontière sont séparés.
//   Field
//   BoundaryField
//   PatchField
//
// Exemple :
//   p.internal
//   p.boundary["inlet"]
//   p.boundary["wall"]
//
// Les conditions limites sont donc attachées aux patches et non au maillage.

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/boundary.h"
#include <vector>
#include <string>
#include <map>
#include <stdexcept>
#include <cstddef>

namespace cfdx {
namespace core {

// Un PatchField est une portion de champ sur un patch spécifique.
class PatchField {
public:
    PatchField() : patch_name_(), type_("zero"), n_faces_(0), data_() {}

    PatchField(const std::string& patch_name,
               std::size_t n_faces,
               const std::string& type = "zero")
        : patch_name_(patch_name),
          type_(type),
          n_faces_(n_faces),
          data_(n_faces, 0.0) {}

    // --- Accès ---

    const std::string& patch_name() const noexcept { return patch_name_; }
    const std::string& type() const noexcept { return type_; }
    std::size_t size() const noexcept { return n_faces_; }
    bool empty() const noexcept { return n_faces_ == 0; }

    double operator()(std::size_t i) const {
        check_index(i);
        return data_[i];
    }
    double& operator()(std::size_t i) {
        check_index(i);
        return data_[i];
    }

    void set_type(const std::string& type) { type_ = type; }
    void fill(double value) { std::fill(data_.begin(), data_.end(), value); }
    void resize(std::size_t n) { n_faces_ = n; data_.resize(n, 0.0); }
    void clear() { n_faces_ = 0; data_.clear(); }

    const double* data() const noexcept { return data_.data(); }
    double* data() noexcept { return data_.data(); }

private:
    void check_index(std::size_t i) const {
        if (i >= n_faces_) {
            throw std::out_of_range("PatchField: index out of range");
        }
    }

    std::string patch_name_;
    std::string type_;
    std::size_t n_faces_;
    std::vector<double> data_;
};

// BoundaryField gère les PatchField pour un champ donné.
// Il est associé à un BoundaryPatches (par son nom de patch).
class BoundaryField {
public:
    BoundaryField() = default;

    // --- Accès ---

    std::size_t n_patches() const noexcept { return patches_.size(); }

    PatchField& patch(std::size_t i) {
        check_index(i);
        return patches_[i];
    }
    const PatchField& patch(std::size_t i) const {
        check_index(i);
        return patches_[i];
    }

    // Trouve un patch par son nom. Retourne n_patches() si absent.
    std::size_t find(const std::string& name) const {
        for (std::size_t i = 0; i < patches_.size(); ++i) {
            if (patches_[i].patch_name() == name) return i;
        }
        return patches_.size();
    }

    bool has_patch(const std::string& name) const {
        return find(name) < patches_.size();
    }

    PatchField& operator[](const std::string& name) {
        const std::size_t idx = find(name);
        if (idx >= patches_.size()) {
            throw std::runtime_error("BoundaryField: patch '" + name + "' not found");
        }
        return patches_[idx];
    }
    const PatchField& operator[](const std::string& name) const {
        const std::size_t idx = find(name);
        if (idx >= patches_.size()) {
            throw std::runtime_error("BoundaryField: patch '" + name + "' not found");
        }
        return patches_[idx];
    }

    // --- Modification ---

    std::size_t add_patch(const std::string& name, std::size_t n_faces,
                          const std::string& type = "zero") {
        if (has_patch(name)) {
            throw std::runtime_error("BoundaryField: duplicate patch '" + name + "'");
        }
        patches_.push_back(PatchField(name, n_faces, type));
        return patches_.size() - 1;
    }

    void remove_patch(std::size_t i) {
        check_index(i);
        patches_.erase(patches_.begin() + static_cast<std::ptrdiff_t>(i));
    }

    void clear() { patches_.clear(); }

    // --- Initialisation depuis un BoundaryPatches ---

void init_from_patches(const BoundaryPatches& bp, const std::string& type = "zero") {
        patches_.clear();
        patches_.reserve(bp.n_patches());
        for (std::size_t i = 0; i < bp.n_patches(); ++i) {
            const auto& p = bp.patch(i);
            patches_.push_back(PatchField(p.name, p.size(), type));
        }
    }

    // --- Validation ---

    // Vérifie que chaque patch a la taille attendue (cohérente avec le mesh).
    bool is_consistent(const BoundaryPatches& bp) const noexcept {
        if (patches_.size() != bp.n_patches()) return false;
        for (std::size_t i = 0; i < patches_.size(); ++i) {
            if (patches_[i].patch_name() != bp.patch(i).name) return false;
            if (patches_[i].size() != bp.patch(i).size()) return false;
        }
        return true;
    }

private:
    void check_index(std::size_t i) const {
        if (i >= n_patches()) {
            throw std::out_of_range("BoundaryField: patch index out of range");
        }
    }

    std::vector<PatchField> patches_;
};

}  // namespace core
}  // namespace cfdx