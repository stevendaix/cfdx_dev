#pragma once

#include "cfdx/physics/radiation.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

// -----------------------------------------------------------------------------
// Radiation material properties
// -----------------------------------------------------------------------------

struct RadiationOpticalProperties {
    double absorption = 0.0;   // 1/m
    double scattering = 0.0;   // 1/m
    double refractive_index = 1.0;
    double emissivity = 1.0;
    double anisotropy = 0.0;   // Henyey-Greenstein g, [-1,1]

    void validate() const {
        if (!std::isfinite(absorption) || !std::isfinite(scattering) ||
            !std::isfinite(refractive_index) || !std::isfinite(emissivity) ||
            !std::isfinite(anisotropy) ||
            absorption < 0.0 || scattering < 0.0 ||
            refractive_index <= 0.0 || emissivity < 0.0 || emissivity > 1.0 ||
            anisotropy < -1.0 || anisotropy > 1.0) {
            throw std::invalid_argument("invalid radiation optical properties");
        }
    }

    double extinction() const { return absorption + scattering; }

    double optical_thickness(double length) const {
        validate();
        if (!std::isfinite(length) || length < 0.0)
            throw std::invalid_argument("radiation path length must be finite and non-negative");
        return extinction() * length;
    }
};

struct RadiationOpticalPropertyField {
    std::vector<RadiationOpticalProperties> cells;

    void validate() const {
        for (const auto& p : cells) p.validate();
    }

    RadiationOpticalProperties& operator[](std::size_t i) { return cells.at(i); }
    const RadiationOpticalProperties& operator[](std::size_t i) const { return cells.at(i); }
};

// -----------------------------------------------------------------------------
// Spectral / non-gray infrastructure
// -----------------------------------------------------------------------------

struct RadiationBand {
    double wavelength_min = 0.0; // m
    double wavelength_max = 0.0; // m
    double weight = 1.0;
    RadiationOpticalProperties properties;

    void validate() const {
        if (!std::isfinite(wavelength_min) || !std::isfinite(wavelength_max) ||
            !std::isfinite(weight) || wavelength_min < 0.0 ||
            wavelength_max <= wavelength_min || weight < 0.0) {
            throw std::invalid_argument("invalid radiation spectral band");
        }
        properties.validate();
    }
};

inline void validate_radiation_bands(const std::vector<RadiationBand>& bands)
{
    if (bands.empty())
        throw std::invalid_argument("at least one radiation spectral band is required");
    double weight_sum = 0.0;
    double previous = 0.0;
    for (const auto& b : bands) {
        b.validate();
        if (b.wavelength_min < previous)
            throw std::invalid_argument("radiation spectral bands must be ordered");
        previous = b.wavelength_max;
        weight_sum += b.weight;
    }
    if (!(weight_sum > 0.0) || !std::isfinite(weight_sum))
        throw std::invalid_argument("radiation spectral band weights must have positive sum");
}

inline double weighted_band_absorption(
    const std::vector<RadiationBand>& bands,
    double temperature)
{
    validate_radiation_bands(bands);
    if (!std::isfinite(temperature) || temperature <= 0.0)
        throw std::invalid_argument("invalid temperature for spectral radiation");
    double numerator = 0.0, denominator = 0.0;
    for (const auto& b : bands) {
        numerator += b.weight * b.properties.absorption;
        denominator += b.weight;
    }
    return numerator / denominator;
}

// -----------------------------------------------------------------------------
// Diffuse-gray DOM wall operator
// -----------------------------------------------------------------------------

inline double diffuse_gray_wall_intensity(
    double emissivity,
    double temperature,
    double incoming_irradiation)
{
    return gray_diffuse_wall_intensity(emissivity, temperature,
                                       incoming_irradiation);
}

inline double hemispherical_irradiation(
    const std::vector<DiscreteDirection>& directions,
    const std::vector<double>& intensities,
    const std::array<double,3>& normal)
{
    if (directions.size() != intensities.size())
        throw std::invalid_argument("DOM direction/intensity size mismatch");
    const double nn = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] +
                                normal[2]*normal[2]);
    if (!std::isfinite(nn) || nn <= 0.0)
        throw std::invalid_argument("wall normal must be finite and non-zero");
    double G = 0.0;
    for (std::size_t m = 0; m < directions.size(); ++m) {
        const auto& d = directions[m];
        const double mu = (d.dx*normal[0] + d.dy*normal[1] + d.dz*normal[2]) / nn;
        if (mu < 0.0) {
            if (!std::isfinite(intensities[m]) || intensities[m] < 0.0)
                throw std::invalid_argument("incoming DOM intensity must be finite and non-negative");
            G += directions[m].weight * (-mu) * intensities[m];
        }
    }
    return G;
}

inline void apply_diffuse_gray_wall(
    const std::vector<DiscreteDirection>& directions,
    const std::vector<double>& incoming_intensities,
    const std::array<double,3>& outward_normal,
    double emissivity,
    double temperature,
    std::vector<double>& boundary_intensities)
{
    validate_discrete_directions(directions);
    if (incoming_intensities.size() != directions.size() ||
        boundary_intensities.size() != directions.size())
        throw std::invalid_argument("DOM wall direction/intensity size mismatch");
    const double G = hemispherical_irradiation(
        directions, incoming_intensities, outward_normal);
    const double Iw = gray_diffuse_wall_intensity(emissivity, temperature, G);
    const double nn = std::sqrt(outward_normal[0]*outward_normal[0] +
                                outward_normal[1]*outward_normal[1] +
                                outward_normal[2]*outward_normal[2]);
    for (std::size_t m = 0; m < directions.size(); ++m) {
        const double mu = (directions[m].dx*outward_normal[0] +
                           directions[m].dy*outward_normal[1] +
                           directions[m].dz*outward_normal[2]) / nn;
        if (mu > 0.0)
            boundary_intensities[m] = Iw;
    }
}

// -----------------------------------------------------------------------------
// View-factor matrix construction for a supplied geometric kernel.
//
// The geometry backend deliberately operates on a generic differential-patch
// representation so it can be fed by any mesh/face accelerator. The kernel
// uses the centroid differential-area approximation and is conservative only
// after closure. It is therefore a deterministic geometry estimator, not a
// replacement for a visibility/ray-tracing backend.
// -----------------------------------------------------------------------------

struct ViewFactorPatch {
    std::array<double,3> center{0.0,0.0,0.0};
    std::array<double,3> normal{0.0,0.0,1.0};
    double area = 0.0;
};

inline double patch_pair_view_factor(
    const ViewFactorPatch& source,
    const ViewFactorPatch& target)
{
    const double rx = target.center[0]-source.center[0];
    const double ry = target.center[1]-source.center[1];
    const double rz = target.center[2]-source.center[2];
    const double r2 = rx*rx+ry*ry+rz*rz;
    if (!std::isfinite(r2) || r2 <= 0.0)
        return 0.0;
    const double r = std::sqrt(r2);
    const double ns = std::sqrt(source.normal[0]*source.normal[0] +
                                 source.normal[1]*source.normal[1] +
                                 source.normal[2]*source.normal[2]);
    const double nt = std::sqrt(target.normal[0]*target.normal[0] +
                                 target.normal[1]*target.normal[1] +
                                 target.normal[2]*target.normal[2]);
    if (!(source.area > 0.0) || !(target.area > 0.0) ||
        !(ns > 0.0) || !(nt > 0.0))
        throw std::invalid_argument("invalid view-factor patch");
    const double cs = (source.normal[0]*rx + source.normal[1]*ry +
                       source.normal[2]*rz)/(ns*r);
    const double ct = -(target.normal[0]*rx + target.normal[1]*ry +
                        target.normal[2]*rz)/(nt*r);
    if (cs <= 0.0 || ct <= 0.0)
        return 0.0;
    return std::max(0.0, std::min(1.0, target.area * cs * ct /
                                        (M_PI * r2)));
}

inline std::vector<double> estimate_view_factor_matrix(
    const std::vector<ViewFactorPatch>& patches)
{
    if (patches.empty())
        throw std::invalid_argument("view-factor patch set is empty");
    const std::size_t n = patches.size();
    std::vector<double> F(n*n, 0.0);
    for (std::size_t i=0; i<n; ++i)
        for (std::size_t j=0; j<n; ++j)
            if (i != j)
                F[i*n+j] = patch_pair_view_factor(patches[i], patches[j]);

    // Enforce enclosure closure while preserving the physically required
    // reciprocity as far as the centroid approximation permits.
    for (std::size_t i=0; i<n; ++i) {
        double sum = 0.0;
        for (std::size_t j=0; j<n; ++j) sum += F[i*n+j];
        if (sum > 1.0) {
            for (std::size_t j=0; j<n; ++j) F[i*n+j] /= sum;
        }
    }
    return F;
}

// -----------------------------------------------------------------------------
// Radiation diagnostics / conservation
// -----------------------------------------------------------------------------

struct RadiationBalance {
    double emitted = 0.0;
    double absorbed = 0.0;
    double net = 0.0;
    double relative_error = 0.0;
};

inline RadiationBalance radiation_balance(double emitted, double absorbed)
{
    if (!std::isfinite(emitted) || !std::isfinite(absorbed) ||
        emitted < 0.0 || absorbed < 0.0)
        throw std::invalid_argument("invalid radiation balance");
    RadiationBalance b;
    b.emitted = emitted;
    b.absorbed = absorbed;
    b.net = emitted - absorbed;
    b.relative_error = std::abs(b.net) /
        std::max(1.0, std::max(emitted, absorbed));
    return b;
}

} // namespace cfdx::physics
