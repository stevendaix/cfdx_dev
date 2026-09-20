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
//   Field<float, Face> phi;
//
// Stockage SoA (Structure of Arrays) pour compatibilité SIMD/GPU :
//   - Scalaire : 1 array contigu de T
//   - Vecteur  : 3 arrays contigus de T (x[], y[], z[])

#pragma once

#include <vector>
#include <cstddef>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include <type_traits>

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

// Types de précision supportés
enum class Precision : std::uint8_t {
    FLOAT32 = 0,
    FLOAT64
};

inline const char* to_string(Precision p) {
    return (p == Precision::FLOAT32) ? "float32" : "float64";
}

inline Precision precision_from_string(const std::string& s) {
    if (s == "float32" || s == "float") return Precision::FLOAT32;
    if (s == "float64" || s == "double") return Precision::FLOAT64;
    return Precision::FLOAT64;
}

// --- Field metadata (§25) ---
struct FieldMetadata {
    std::string name;
    Location location = Location::UNKNOWN;
    std::size_t dimension = 1;   // scalaire = 1, vecteur = 3
    std::string unit;
    Precision precision = Precision::FLOAT64;
    std::string storage;         // "host", "device", "working_set"
};

// Vector type pour stockage SoA
struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
    Vec3() = default;
    Vec3(double v) : x(v), y(v), z(v) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }

    double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    double mag() const { return std::sqrt(x * x + y * y + z * z); }
    double mag2() const { return x * x + y * y + z * z; }

    Vec3 normalized() const {
        double m = mag();
        if (m == 0.0) return {0.0, 0.0, 0.0};
        return {x / m, y / m, z / m};
    }
};

// Helper pour détecter si T est un type vectoriel (Vec3, etc.)
template <typename T>
struct is_vector_type : std::false_type {};

template <>
struct is_vector_type<Vec3> : std::true_type {};

// --- Field générique (§24) ---
// Stockage SoA (Structure of Arrays) :
//   - Pour scalaire (dim=1) : data_[0] = array de n éléments
//   - Pour vecteur (dim=3) : data_[0] = x[], data_[1] = y[], data_[2] = z[]
// Le template T contrôle le type de précision (float, double, etc.)

template <typename T = double, Location L = Location::CELL>
class Field {
public:
    using ValueType = T;
    static constexpr Location location = L;

    static_assert(std::is_arithmetic_v<T> || std::is_same_v<T, Vec3>,
                  "Field value type must be arithmetic or Vec3");

    Field() = default;

    explicit Field(std::size_t n, const std::string& name = "",
                   const std::string& unit = "", std::size_t dim = 1,
                   Precision precision = Precision::FLOAT64)
        : meta_{name, L, dim, unit, precision, "host"},
          n_(n),
          dim_(dim)
    {
        if (dim_ == 0) dim_ = 1;
        allocate_storage();
    }

    // --- Accès ---

    std::size_t size() const noexcept { return n_; }
    std::size_t dimension() const noexcept { return dim_; }
    void set_dimension(std::size_t dim) {
        if (dim != dim_) {
            dim_ = dim;
            allocate_storage();
        }
    }
    bool empty() const noexcept { return n_ == 0; }

    Precision precision() const noexcept { return meta_.precision; }

    // Accès scalaire (dim == 1).
    T operator()(std::size_t i) const {
        check_index(i);
        return data_[0][i];
    }
    T& operator()(std::size_t i) {
        check_index(i);
        return data_[0][i];
    }

    // Accès vecteur (dim == 3) - SoA layout
    void get(std::size_t i, T& x, T& y, T& z) const {
        check_index(i);
        if (dim_ < 3) {
            throw std::runtime_error("Field::get: dimension < 3, cannot get vector components");
        }
        x = data_[0][i];
        y = data_[1][i];
        z = data_[2][i];
    }

    void set(std::size_t i, T x, T y, T z) {
        check_index(i);
        if (dim_ < 3) {
            throw std::runtime_error("Field::set: dimension < 3, cannot set vector components");
        }
        data_[0][i] = x;
        data_[1][i] = y;
        data_[2][i] = z;
    }

    // Accès par composante pour vecteur
    T& operator()(std::size_t i, std::size_t comp) {
        check_index(i);
        if (comp >= dim_) {
            throw std::out_of_range("Field: component index out of range");
        }
        return data_[comp][i];
    }

    T operator()(std::size_t i, std::size_t comp) const {
        check_index(i);
        if (comp >= dim_) {
            throw std::out_of_range("Field: component index out of range");
        }
        return data_[comp][i];
    }

    // --- Métadonnées ---

    const FieldMetadata& metadata() const noexcept { return meta_; }
    FieldMetadata& metadata() { return meta_; }
    const std::string& name() const noexcept { return meta_.name; }
    Location loc() const noexcept { return meta_.location; }

    // --- Opérations ---

    void fill(T value) {
        for (std::size_t c = 0; c < dim_; ++c) {
            std::fill(data_[c].begin(), data_[c].end(), value);
        }
    }

    void resize(std::size_t n) {
        n_ = n;
        allocate_storage();
    }

    void clear() {
        n_ = 0;
        for (auto& arr : data_) arr.clear();
    }

    // --- Accès bulk (pour kernels) ---

    // Retourne pointeurs vers les arrays SoA (data_[0], data_[1], data_[2] pour vecteur)
    const T* const* data() const noexcept { return data_ptrs_.data(); }
    T* const* data() noexcept { return data_ptrs_.data(); }

    // Accès direct à une composante
    const T* component_data(std::size_t comp) const noexcept {
        return (comp < dim_) ? data_[comp].data() : nullptr;
    }
    T* component_data(std::size_t comp) noexcept {
        return (comp < dim_) ? data_[comp].data() : nullptr;
    }

    // --- Validation ---

    bool is_valid() const {
        for (std::size_t c = 0; c < dim_; ++c) {
            for (std::size_t i = 0; i < n_; ++i) {
                if (!(std::isfinite(data_[c][i]))) return false;
            }
        }
        return true;
    }

private:
    void check_index(std::size_t i) const {
        if (i >= n_) {
            throw std::out_of_range("Field: index out of range");
        }
    }

    void allocate_storage() {
        data_.clear();
        data_.resize(dim_);
        for (std::size_t c = 0; c < dim_; ++c) {
            data_[c].resize(n_, T{0});
        }
        // Update pointer array for kernel access
        data_ptrs_.resize(dim_);
        for (std::size_t c = 0; c < dim_; ++c) {
            data_ptrs_[c] = data_[c].data();
        }
    }

    FieldMetadata meta_;
    std::size_t n_ = 0;
    std::size_t dim_ = 1;
    std::vector<std::vector<T>> data_;      // SoA: data_[comp][index]
    std::vector<T*> data_ptrs_;              // Pointers for kernel access
};

// Type aliases courants (§24).
using ScalarCellField = Field<double, Location::CELL>;
using ScalarFaceField = Field<double, Location::FACE>;
using ScalarPointField = Field<double, Location::POINT>;
using Vec3CellField = Field<Vec3, Location::CELL>;   // T=Vec3, dim=3 (stored as 3 arrays)
using Vec3FaceField = Field<Vec3, Location::FACE>;
using Float32CellField = Field<float, Location::CELL>;
using Float64CellField = Field<double, Location::CELL>;

}  // namespace core
}  // namespace cfdx
