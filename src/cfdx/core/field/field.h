// M0.4-T01 — Field<T, Location> template
//
// Spécification CFDX v0.7 §24 :
//   Les champs sont génériques.
//   Concept : Field<T, Location>
//   Locations : CELL, FACE, POINT, BOUNDARY
//
// Exemples :
//   Field<double, Cell> p;
//   Field<Vec3, Cell> U;
//   Field<double, Face> phi;

#pragma once

#include <vector>
#include <cstddef>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <cmath>

namespace cfdx {
namespace core {

enum class Location : std::uint8_t {
    CELL = 0,
    FACE,
    POINT,
    BOUNDARY,
    UNKNOWN
};

inline const char* to_string(Location loc) {
    switch (loc) {
        case Location::CELL:     return "cell";
        case Location::FACE:     return "face";
        case Location::POINT:    return "point";
        case Location::BOUNDARY: return "boundary";
        default:                 return "unknown";
    }
}

inline Location location_from_string(const std::string& s) {
    if (s == "cell")     return Location::CELL;
    if (s == "face")     return Location::FACE;
    if (s == "point")    return Location::POINT;
    if (s == "boundary") return Location::BOUNDARY;
    return Location::UNKNOWN;
}

// --- Field metadata (§25) ---
struct FieldMetadata {
    std::string name;
    Location location = Location::UNKNOWN;
    std::size_t dimension = 1;   // scalaire = 1, vecteur = 3
    std::string unit;
    std::string precision;       // "float32" ou "float64"
    std::string storage;         // "host", "device", "working_set"
};

// --- Field générique (§24) ---
// Le stockage est SoA pour les vecteurs (§38).
// Pour un champ scalaire, T = double, storage est un std::vector<double>.
// Pour un champ vecteur (Vec3), storage est 3 std::vector<double> (x, y, z).

template <typename T, Location L>
class Field {
public:
    using Value = T;
    static constexpr Location location = L;

    Field() = default;

    explicit Field(std::size_t n, const std::string& name = "",
                   const std::string& unit = "", std::size_t dim = 1)
        : meta_{name, L, dim, unit, "float64", "host"},
          n_(n),
          dim_(dim)
    {
        if (dim_ == 0) dim_ = 1;
        data_.resize(n_ * dim_, 0.0);
    }

    // --- Accès ---

    std::size_t size() const noexcept { return n_; }
    std::size_t dimension() const noexcept { return dim_; }
    bool empty() const noexcept { return n_ == 0; }

    // Accès scalaire (dim == 1).
    double operator()(std::size_t i) const {
        check_index(i);
        return data_[i];
    }
    double& operator()(std::size_t i) {
        check_index(i);
        return data_[i];
    }

    // Accès vecteur (dim == 3).
    void get(std::size_t i, double& x, double& y, double& z) const {
        check_index(i);
        const std::size_t base = i * dim_;
        x = data_[base];
        y = data_[base + 1];
        z = data_[base + 2];
    }
    void set(std::size_t i, double x, double y, double z) {
        check_index(i);
        const std::size_t base = i * dim_;
        data_[base] = x;
        data_[base + 1] = y;
        data_[base + 2] = z;
    }

    // --- Métadonnées ---

    const FieldMetadata& metadata() const noexcept { return meta_; }
    FieldMetadata& metadata() { return meta_; }
    const std::string& name() const noexcept { return meta_.name; }
    Location loc() const noexcept { return meta_.location; }

    // --- Opérations ---

    void fill(double value) {
        std::fill(data_.begin(), data_.end(), value);
    }

    void resize(std::size_t n) {
        n_ = n;
        data_.resize(n_ * dim_, 0.0);
    }

    void clear() {
        n_ = 0;
        data_.clear();
    }

    // --- Accès bulk ---

    const double* data() const noexcept { return data_.data(); }
    double* data() noexcept { return data_.data(); }

    // --- Validation ---

    bool is_valid() const {
        for (std::size_t i = 0; i < data_.size(); ++i) {
            if (!(std::isfinite(data_[i]))) return false;
        }
        return true;
    }

private:
    void check_index(std::size_t i) const {
        if (i >= n_) {
            throw std::out_of_range("Field: index out of range");
        }
    }

    FieldMetadata meta_;
    std::size_t n_ = 0;
    std::size_t dim_ = 1;
    std::vector<double> data_;
};

// Types courants (§24).
using ScalarCellField = Field<double, Location::CELL>;
using ScalarFaceField = Field<double, Location::FACE>;
using ScalarPointField = Field<double, Location::POINT>;
using Vec3CellField = Field<double, Location::CELL>;  // dim=3
using Vec3FaceField = Field<double, Location::FACE>;   // dim=3

}  // namespace core
}  // namespace cfdx