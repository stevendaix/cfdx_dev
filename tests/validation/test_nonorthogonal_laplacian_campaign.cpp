// Issue #461 — N3 quantitative non-orthogonal Laplacian campaign.
//
// The N3 unit tests pin the over-relaxed invariants on a two-cell case. This
// campaign quantifies the diffusion operators on a controlled
// non-orthogonality ladder, which is what is required before over-relaxed can
// be promoted from "implemented" to "verified":
//
//   linear manufactured field (exact Laplacian = 0)
//     -> uncorrected/orthogonal error grows with skew;
//     -> corrected/limited/over-relaxed stay at machine zero.
//
//   smooth manufactured field (exact Laplacian = -3 pi^2 phi)
//     -> measured observed order on interior cells for each scheme.

#include "cfdx/core/numerics/laplacian.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/field/field.h"
#include "verification_metrics.h"

#include <algorithm>
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

Grid make_affine_cube(std::size_t n, double shear)
{
    if (n < 2) throw std::invalid_argument("make_affine_cube: n must be >= 2");

    const double h = 1.0 / static_cast<double>(n);
    const std::size_t nc = n * n * n;
    Mesh m;
    m.points().resize(8 * nc);

    const auto cell_id = [n](std::size_t i, std::size_t j, std::size_t k) {
        return (k * n + j) * n + i;
    };
    const auto transform = [shear](double x, double y, double z) {
        return std::array<double, 3>{x + shear * y, y, z};
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
    // Deep interior: cells whose neighbours are all interior too. The
    // non-orthogonal correction uses the neighbour cell gradients, and
    // boundary-cell Green-Gauss gradients carry the zero-gradient boundary
    // approximation, so including boundary-adjacent cells would measure that
    // approximation rather than the diffusion operator under test.
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i)
                if (i >= 2 && j >= 2 && k >= 2 && i + 2 < n && j + 2 < n && k + 2 < n)
                    grid.interior.push_back(cell_id(i, j, k));
    return grid;
}

double smooth_value(const Vec3& p) {
    const double pi = std::acos(-1.0);
    return std::sin(pi * p.x) * std::cos(pi * p.y) * std::sin(pi * p.z);
}
double smooth_laplacian(const Vec3& p) {
    const double pi = std::acos(-1.0);
    return -3.0 * pi * pi * smooth_value(p);
}

std::vector<double> sample(const Grid& grid, double (*fn)(const Vec3&)) {
    std::vector<double> out(grid.mesh.n_cells(), 0.0);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c)
        out[c] = fn(grid.geometry.cell_centres[c]);
    return out;
}

ErrorMetrics interior_error(const Grid& grid, const Field<double, Location::CELL>& lap,
                            const std::vector<double>& exact)
{
    std::vector<double> got, want, weights;
    for (const std::size_t c : grid.interior) {
        got.push_back(lap(c));
        want.push_back(exact[c]);
        weights.push_back(grid.geometry.cell_volumes[c]);
    }
    return error_norms(got, want, weights);
}

Field<double, Location::CELL> laplacian_of(const Grid& grid, const std::vector<double>& phi,
                                           LaplacianScheme scheme)
{
    ScalarCellField f(grid.mesh.n_cells(), "phi", "1", 1);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c) f(c) = phi[c];
    return compute_laplacian(f, grid.mesh, grid.geometry, scheme, 0.5);
}

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

// Linear field: exact Laplacian is zero. Corrected schemes must be exact at the
// machine level on every affine mesh; the uncorrected two-point operator must
// degrade as the skew grows.
void linear_campaign()
{
    const std::vector<double> skews = {0.0, 0.25, 0.5, 1.0, 2.0};
    for (const double shear : skews) {
        const Grid grid = make_affine_cube(8, shear);
        std::vector<double> phi(grid.mesh.n_cells(), 0.0);
        for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c)
            phi[c] = grid.geometry.cell_centres[c].x;  // Laplacian = 0

        const std::vector<double> zero(grid.mesh.n_cells(), 0.0);
        const auto orth = interior_error(grid, laplacian_of(grid, phi, LaplacianScheme::ORTHOGONAL), zero);
        const auto corr = interior_error(grid, laplacian_of(grid, phi, LaplacianScheme::CORRECTED), zero);
        const auto lim = interior_error(grid, laplacian_of(grid, phi, LaplacianScheme::LIMITED), zero);
        const auto over = interior_error(grid, laplacian_of(grid, phi, LaplacianScheme::OVER_RELAXED), zero);

        std::cout << "N3_LINEAR shear=" << shear
                  << " uncorrected_linf=" << orth.linf
                  << " corrected_linf=" << corr.linf
                  << " limited_linf=" << lim.linf
                  << " over_relaxed_linf=" << over.linf << "\n";

        const double tol = 1e-9;
        require(corr.linf <= tol, "corrected must be linear-exact under skew");
        require(lim.linf <= tol, "limited must be linear-exact under skew");
        require(over.linf <= tol, "over-relaxed must be linear-exact under skew");
    }
}

void smooth_campaign()
{
    const std::vector<double> skews = {0.0, 0.5, 1.0};
    const std::vector<std::string> labels = {"orthogonal", "moderate", "strong"};
    for (std::size_t s = 0; s < skews.size(); ++s) {
        std::vector<double> err_uncorrected, err_corrected, err_over;
        for (const std::size_t n : {8u, 16u, 32u}) {
            const Grid grid = make_affine_cube(n, skews[s]);
            const auto exact = sample(grid, smooth_laplacian);
            const auto phi = sample(grid, smooth_value);
            const auto e_u = interior_error(grid, laplacian_of(grid, phi, LaplacianScheme::UNCORRECTED), exact);
            const auto e_c = interior_error(grid, laplacian_of(grid, phi, LaplacianScheme::CORRECTED), exact);
            const auto e_o = interior_error(grid, laplacian_of(grid, phi, LaplacianScheme::OVER_RELAXED), exact);
            err_uncorrected.push_back(e_u.l2);
            err_corrected.push_back(e_c.l2);
            err_over.push_back(e_o.l2);
            std::cout << "N3_SMOOTH level=" << labels[s] << " shear=" << skews[s] << " n=" << n
                      << " uncorrected_L2=" << e_u.l2 << " corrected_L2=" << e_c.l2
                      << " over_relaxed_L2=" << e_o.l2;
            if (err_corrected.size() > 1)
                std::cout << " corrected_order="
                          << observed_order(err_corrected[err_corrected.size() - 2], err_corrected.back())
                          << " over_relaxed_order="
                          << observed_order(err_over[err_over.size() - 2], err_over.back());
            std::cout << "\n";
        }
        require(err_corrected.back() < err_corrected.front(),
                "corrected smooth refinement must reduce the error");
        require(err_over.back() < err_over.front(),
                "over-relaxed smooth refinement must reduce the error");
        require_order(err_corrected, 1.0, 1.0, labels[s] + "/corrected");
        require_order(err_over, 1.0, 1.0, labels[s] + "/over_relaxed");
        // Under skew the uncorrected two-point operator must lose accuracy
        // relative to the corrected/over-relaxed operators on a curved field.
        if (skews[s] > 0.0) {
            require(err_uncorrected.back() > err_corrected.back(),
                    labels[s] + ": uncorrected must be worse than corrected under skew");
            require(err_uncorrected.back() > err_over.back(),
                    labels[s] + ": uncorrected must be worse than over-relaxed under skew");
        }
    }
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(12);
        linear_campaign();
        smooth_campaign();
        std::cout << "NONORTHOGONAL_LAPLACIAN_CAMPAIGN: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "NONORTHOGONAL_LAPLACIAN_CAMPAIGN: FAIL: " << e.what() << "\n";
        return 1;
    }
}
