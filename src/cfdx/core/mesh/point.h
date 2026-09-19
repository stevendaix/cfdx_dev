#pragma once

#include <vector>
#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace cfdx {
namespace core {

// ---------------------------------------------------------------------------
// M0.1-T01 — Point storage
//
// Spécification CFDX v0.7 §10 :
//   /mesh/points
//   float64[nPoints, 3]
//   Chaque point possède une position x y z.
//
// Le stockage est SoA (Structure of Arrays) pour une excellente localité
// mémoire dans les kernels FVM (§38).
// ---------------------------------------------------------------------------

class PointCloud {
public:
    PointCloud() = default;

    explicit PointCloud(std::size_t n_points)
        : x_(n_points, 0.0),
          y_(n_points, 0.0),
          z_(n_points, 0.0) {}

    // --- Accès par index ---

    double x(std::size_t i) const {
        check_index(i);
        return x_[i];
    }

    double y(std::size_t i) const {
        check_index(i);
        return y_[i];
    }

    double z(std::size_t i) const {
        check_index(i);
        return z_[i];
    }

    void set(std::size_t i, double x, double y, double z) {
        check_index(i);
        x_[i] = x;
        y_[i] = y;
        z_[i] = z;
    }

    // --- Accès bulk (pour les kernels) ---

    const double* x_data() const noexcept { return x_.data(); }
    const double* y_data() const noexcept { return y_.data(); }
    const double* z_data() const noexcept { return z_.data(); }

    double* x_data() noexcept { return x_.data(); }
    double* y_data() noexcept { return y_.data(); }
    double* z_data() noexcept { return z_.data(); }

    // --- Dimensions ---

    std::size_t size() const noexcept { return x_.size(); }
    bool empty() const noexcept { return x_.empty(); }

    void resize(std::size_t n) {
        x_.resize(n, 0.0);
        y_.resize(n, 0.0);
        z_.resize(n, 0.0);
    }

    void clear() {
        x_.clear();
        y_.clear();
        z_.clear();
    }

    // --- Validation ---

    // Vérifie qu'aucun point n'est NaN ou Inf (§18, §19).
    bool is_valid() const {
        for (std::size_t i = 0; i < size(); ++i) {
            if (!is_finite(x_[i]) || !is_finite(y_[i]) || !is_finite(z_[i])) {
                return false;
            }
        }
        return true;
    }

private:
    static bool is_finite(double v) noexcept {
        return !(std::isnan(v) || std::isinf(v));
    }

    void check_index(std::size_t i) const {
        if (i >= size()) {
            throw std::out_of_range("PointCloud: index out of range");
        }
    }

    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> z_;
};

}  // namespace core
}  // namespace cfdx