// Issue #461 — N2 gradient/reconstruction V&V campaign.
//
// Closes the "complete skew/non-orthogonal/high-aspect-ratio/polyhedral V&V"
// gap for the gradient family by running the same gradient operators over an
// orthogonal, an affine-sheared (non-orthogonal) and a stretched
// (high-aspect-ratio) mesh family:
//
//   - constant field  -> exact zero gradient (Green-Gauss and least-squares);
//   - linear field    -> exact interior gradient on every affine family;
//   - quadratic/smooth fields on a refinement sequence -> measured observed
//     order on interior cells with L1/L2/Linf reporting.
//
// Polyhedral remains out of scope here (see the #461 audit); this test makes
// the skew/non-orthogonal/high-aspect-ratio part explicit rather than relying
// on the Laplacian-only skew campaign that was previously cited as N2 evidence.

#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/field/field.h"
#include "verification_metrics.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::verification;

namespace {

struct Grid {
    Mesh mesh;
    GeometryCache geometry;
    std::vector<std::size_t> interior;
};

// Unit-cube Cartesian topology, then an affine map applied to every point:
//   x' = x + shear * y,  y' = stretch * y,  z' = z
// Affine transforms preserve the closed-cell geometric identity that makes
// Green-Gauss linear-exact, while shear introduces non-orthogonality and
// stretch introduces high aspect ratio.
Grid make_affine_cube(std::size_t n, double shear, double stretch)
{
    if (n < 2) throw std::invalid_argument("make_affine_cube: n must be >= 2");
    if (!(stretch > 0.0)) throw std::invalid_argument("make_affine_cube: stretch must be positive");

    const double h = 1.0 / static_cast<double>(n);
    const std::size_t nc = n * n * n;
    Mesh m;
    m.points().resize(8 * nc);

    const auto cell_id = [n](std::size_t i, std::size_t j, std::size_t k) {
        return (k * n + j) * n + i;
    };
    const auto transform = [shear, stretch](double x, double y, double z) {
        return std::array<double, 3>{x + shear * y, stretch * y, z};
    };

    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t c = cell_id(i, j, k);
                const std::size_t b = 8 * c;
                const double x0 = i * h, x1 = (i + 1) * h;
                const double y0 = j * h, y1 = (j + 1) * h;
                const double z0 = k * h, z1 = (k + 1) * h;
                const double raw[8][3] = {
                    {x0,y0,z0}, {x1,y0,z0}, {x1,y1,z0}, {x0,y1,z0},
                    {x0,y0,z1}, {x1,y0,z1}, {x1,y1,z1}, {x0,y1,z1}
                };
                for (std::size_t q = 0; q < 8; ++q) {
                    const auto p = transform(raw[q][0], raw[q][1], raw[q][2]);
                    m.points().set(b + q, p[0], p[1], p[2]);
                }
            }
        }
    }

    std::vector<std::vector<std::size_t>> faces_per_cell(nc);

    struct FaceKey {
        std::size_t i, j, k, axis;
        bool operator==(const FaceKey& other) const {
            return i == other.i && j == other.j && k == other.k && axis == other.axis;
        }
    };
    struct FaceKeyHash {
        std::size_t operator()(const FaceKey& f) const {
            std::size_t h = f.i;
            h = h * 1315423911u + f.j;
            h = h * 1315423911u + f.k;
            h = h * 1315423911u + f.axis;
            return h;
        }
    };
    std::unordered_map<FaceKey, std::size_t, FaceKeyHash> face_map;

    auto add_or_get_face = [&](FaceKey key, std::initializer_list<std::size_t> vertices) {
        auto it = face_map.find(key);
        if (it != face_map.end()) return it->second;
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face(std::vector<FaceIndex>(vertices.begin(), vertices.end()));
        face_map.emplace(key, id);
        return id;
    };

    const auto point = [](std::size_t c, std::size_t q) { return 8 * c + q; };

    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t c = cell_id(i, j, k);
                auto face = [&](std::size_t axis, std::size_t fi, std::size_t fj, std::size_t fk,
                                std::initializer_list<std::size_t> v) {
                    return add_or_get_face({fi, fj, fk, axis}, v);
                };
                const std::size_t xm = face(0, i, j, k, {point(c,0),point(c,4),point(c,7),point(c,3)});
                const std::size_t xp = face(0, i + 1, j, k, {point(c,1),point(c,2),point(c,6),point(c,5)});
                const std::size_t ym = face(1, i, j, k, {point(c,0),point(c,1),point(c,5),point(c,4)});
                const std::size_t yp = face(1, i, j + 1, k, {point(c,3),point(c,7),point(c,6),point(c,2)});
                const std::size_t zm = face(2, i, j, k, {point(c,0),point(c,3),point(c,2),point(c,1)});
                const std::size_t zp = face(2, i, j, k + 1, {point(c,4),point(c,5),point(c,6),point(c,7)});
                faces_per_cell[c] = {xm, xp, ym, yp, zm, zp};
            }
        }
    }

    m.ownership().resize(m.n_faces());
    for (const auto& kv : face_map) {
        const FaceKey key = kv.first;
        const std::size_t f = kv.second;
        const bool lower_boundary = (key.axis == 0 ? key.i == 0 :
                                     key.axis == 1 ? key.j == 0 : key.k == 0);
        const bool upper_boundary = (key.axis == 0 ? key.i == n :
                                     key.axis == 1 ? key.j == n : key.k == n);
        const std::size_t owner_i = key.axis == 0 ? (lower_boundary ? 0 : key.i - 1) : key.i;
        const std::size_t owner_j = key.axis == 1 ? (lower_boundary ? 0 : key.j - 1) : key.j;
        const std::size_t owner_k = key.axis == 2 ? (lower_boundary ? 0 : key.k - 1) : key.k;
        const std::size_t owner = cell_id(owner_i, owner_j, owner_k);
        m.ownership().set_owner(f, owner);
        if (lower_boundary || upper_boundary) {
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        } else {
            std::size_t ni = owner_i, nj = owner_j, nk = owner_k;
            if (key.axis == 0) ++ni;
            else if (key.axis == 1) ++nj;
            else ++nk;
            m.ownership().set_neighbour(f, static_cast<std::int64_t>(cell_id(ni, nj, nk)));
        }
    }

    for (std::size_t c = 0; c < nc; ++c) m.cells().push_cell(faces_per_cell[c]);

    Grid grid;
    grid.mesh = std::move(m);
    grid.geometry = make_geometry_cache(grid.mesh);
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i)
                if (i >= 1 && j >= 1 && k >= 1 && i + 1 < n && j + 1 < n && k + 1 < n)
                    grid.interior.push_back(cell_id(i, j, k));
    return grid;
}


} // namespace
