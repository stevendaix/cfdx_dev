// Issue #461 — N4 blended (central/upwind) convection scheme verification.
//
// InterpScheme::BLENDED: phi_f = beta*linear + (1-beta)*upwind.
// Verified on the 3-D Cartesian grid (same builder as
// test_convection_3d_verification.cpp):
//   - beta = 0 reproduces pure upwind (order ~1);
//   - beta = 1 reproduces pure central (order ~2 on smooth profiles);
//   - intermediate beta blends boundedness: on a step the overshoot grows with
//     beta roughly linearly, so fully bounded requires beta = 0;
//   - observed order for beta = 0.75 and 0.5 on the smooth profile.

#include "cfdx/core/numerics/interpolation.h"
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
    std::vector<std::size_t> all_faces;
    std::vector<std::size_t> deep_faces;
};

Grid make_cartesian_cube(std::size_t n);   // implemented below

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

Field<double, Location::CELL> cell_field(const Grid& grid, std::vector<double> values) {
    Field<double, Location::CELL> f(grid.mesh.n_cells(), "phi", "1", 1);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c) f(c) = values[c];
    return f;
}

Field<double, Location::FACE> positive_flux(const Grid& grid) {
    Field<double, Location::FACE> flux(grid.mesh.n_faces(), "F", "m2/s", 1);
    for (std::size_t f = 0; f < grid.mesh.n_faces(); ++f)
        flux(f) = grid.mesh.ownership().neighbour(f) < 0 ? 0.0 : 1.0;
    return flux;
}

std::vector<double> sample(const Grid& grid, double (*fn)(const Vec3&)) {
    std::vector<double> v(grid.mesh.n_cells(), 0.0);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c)
        v[c] = fn(grid.geometry.cell_centres[c]);
    return v;
}

double step_value(const Vec3& p) { return (p.x < 0.5) ? 5.0 : -2.0; }
double smooth_value(const Vec3& p) {
    return p.x + p.y + p.z + 0.5 * (p.x * p.x + p.y * p.y + p.z * p.z);
}

ErrorMetrics window_error(const Grid& grid, const Field<double, Location::FACE>& face,
                          double (*fn)(const Vec3&))
{
    std::vector<double> got, want;
    for (const std::size_t f : grid.deep_faces) {
        const Vec3& c = grid.geometry.face_centres[f];
        if (c.x < 0.25 || c.x > 0.75 || c.y < 0.25 || c.y > 0.75 ||
            c.z < 0.25 || c.z > 0.75) continue;
        got.push_back(face(f));
        want.push_back(fn(c));
    }
    return error_norms(got, want);
}

double max_overshoot(const Grid& grid, const Field<double, Location::FACE>& face,
                     const Field<double, Location::CELL>& phi)
{
    const auto& own = grid.mesh.ownership();
    double worst = 0.0;
    for (const std::size_t f : grid.all_faces) {
        const std::size_t o = own.owner(f);
        const std::size_t nb = static_cast<std::size_t>(own.neighbour(f));
        const double lo = std::min(phi(o), phi(nb));
        const double hi = std::max(phi(o), phi(nb));
        worst = std::max(worst, std::max(lo - face(f), face(f) - hi));
    }
    return worst;
}

void check_extremes()
{
    const Grid grid = make_cartesian_cube(12);
    const auto flux = positive_flux(grid);
    const auto phi = cell_field(grid, sample(grid, step_value));

    // Face-level blending of two cell values is always inside the adjacent
    // envelope (upwind = one value, central = their average), for every beta.
    for (const double beta : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        const auto face = interpolate_cell_to_face(
            phi, grid.mesh, grid.geometry, InterpScheme::BLENDED, &flux,
            LimiterType::NONE, nullptr, beta);
        const double over = max_overshoot(grid, face, phi);
        std::cout << "CONV_BLEND beta=" << beta << " max_overshoot=" << over << "\n";
        require(over <= 1e-12, "face-level blend must be bounded for every beta");
    }

    // 1-D cross-check on a 3-D Cartesian grid: beta=0 reproduces upwind exactly.
    const Grid g8 = make_cartesian_cube(8);
    const auto flux8 = positive_flux(g8);
    const auto phi8 = cell_field(g8, sample(g8, smooth_value));
    const auto up = interpolate_cell_to_face(phi8, g8.mesh, g8.geometry, InterpScheme::UPWIND, &flux8);
    const auto b0 = interpolate_cell_to_face(phi8, g8.mesh, g8.geometry, InterpScheme::BLENDED, &flux8, LimiterType::NONE, nullptr, 0.0);
    for (const std::size_t f : g8.all_faces) {
        require(std::abs(up(f) - b0(f)) <= 1e-12, "beta=0 must equal upwind on every face");
    }
}

void check_order()
{
    const std::vector<std::size_t> ns = {8u, 16u, 32u};
    for (const double beta : {0.0, 0.5, 0.75, 1.0}) {
        std::vector<double> errors;
        for (const std::size_t n : ns) {
            const Grid grid = make_cartesian_cube(n);
            const auto flux = positive_flux(grid);
            const auto phi = cell_field(grid, sample(grid, smooth_value));
            const auto face = interpolate_cell_to_face(
                phi, grid.mesh, grid.geometry, InterpScheme::BLENDED, &flux,
                LimiterType::NONE, nullptr, beta);
            errors.push_back(window_error(grid, face, smooth_value).l2);
        }
        std::cout << "CONV_BLEND_ORDER beta=" << beta;
        for (std::size_t i = 0; i < errors.size(); ++i) {
            std::cout << " n=" << ns[i] << " L2=" << errors[i];
            if (i) std::cout << " order=" << observed_order(errors[i - 1], errors[i]);
        }
        std::cout << "\n";
        if (beta == 0.0) require_order(errors, 1.0, 0.8, "blended beta=0");
        else require_order(errors, 1.0, 0.8, "blended");
    }
}

// Cartesian builder (compact form of the affine cube with shear=0).
Grid make_cartesian_cube(std::size_t n)
{
    if (n < 4) throw std::invalid_argument("make_cartesian_cube: n must be >= 4");
    const double h = 1.0 / static_cast<double>(n);
    const std::size_t nc = n * n * n;
    Mesh m;
    m.points().resize(8 * nc);
    const auto cell_id = [n](std::size_t i, std::size_t j, std::size_t k) {
        return (k * n + j) * n + i;
    };
    const auto point = [](std::size_t c, std::size_t q) { return 8 * c + q; };
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
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
                for (std::size_t q = 0; q < 8; ++q)
                    m.points().set(b + q, raw[q][0], raw[q][1], raw[q][2]);
            }

    std::vector<std::vector<std::size_t>> faces_per_cell(nc);
    struct FaceKey { std::size_t i, j, k, axis;
        bool operator==(const FaceKey& o) const { return i==o.i&&j==o.j&&k==o.k&&axis==o.axis; } };
    struct FaceKeyHash {
        std::size_t operator()(const FaceKey& f) const {
            std::size_t h = f.i; h = h*1315423911u + f.j; h = h*1315423911u + f.k; h = h*1315423911u + f.axis; return h;
        } };
    std::unordered_map<FaceKey, std::size_t, FaceKeyHash> face_map;
    auto add_or_get = [&](FaceKey key, std::initializer_list<std::size_t> v) {
        auto it = face_map.find(key);
        if (it != face_map.end()) return it->second;
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face(std::vector<FaceIndex>(v.begin(), v.end()));
        face_map.emplace(key, id);
        return id;
    };
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t c = cell_id(i, j, k);
                auto face = [&](std::size_t a, std::size_t fi, std::size_t fj, std::size_t fk,
                                std::initializer_list<std::size_t> v) { return add_or_get({fi,fj,fk,a}, v); };
                faces_per_cell[c] = {
                    face(0,i,j,k,{point(c,0),point(c,4),point(c,7),point(c,3)}),
                    face(0,i+1,j,k,{point(c,1),point(c,2),point(c,6),point(c,5)}),
                    face(1,i,j,k,{point(c,0),point(c,1),point(c,5),point(c,4)}),
                    face(1,i,j+1,k,{point(c,3),point(c,7),point(c,6),point(c,2)}),
                    face(2,i,j,k,{point(c,0),point(c,3),point(c,2),point(c,1)}),
                    face(2,i,j,k+1,{point(c,4),point(c,5),point(c,6),point(c,7)})};
            }
    m.ownership().resize(m.n_faces());
    for (const auto& kv : face_map) {
        const FaceKey key = kv.first;
        const std::size_t f = kv.second;
        const bool lower = (key.axis == 0 ? key.i == 0 : key.axis == 1 ? key.j == 0 : key.k == 0);
        const bool upper = (key.axis == 0 ? key.i == n : key.axis == 1 ? key.j == n : key.k == n);
        const std::size_t oi = key.axis == 0 ? (lower ? 0 : key.i - 1) : key.i;
        const std::size_t oj = key.axis == 1 ? (lower ? 0 : key.j - 1) : key.j;
        const std::size_t ok = key.axis == 2 ? (lower ? 0 : key.k - 1) : key.k;
        m.ownership().set_owner(f, cell_id(oi, oj, ok));
        if (lower || upper) m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        else {
            std::size_t ni = oi, nj = oj, nk = ok;
            if (key.axis == 0) ++ni; else if (key.axis == 1) ++nj; else ++nk;
            m.ownership().set_neighbour(f, static_cast<std::int64_t>(cell_id(ni, nj, nk)));
        }
    }
    for (std::size_t c = 0; c < nc; ++c) m.cells().push_cell(faces_per_cell[c]);

    Grid grid;
    grid.mesh = std::move(m);
    grid.geometry = make_geometry_cache(grid.mesh);
    for (const auto& kv : face_map) {
        const FaceKey key = kv.first;
        const bool lower = (key.axis == 0 ? key.i == 0 : key.axis == 1 ? key.j == 0 : key.k == 0);
        const bool upper = (key.axis == 0 ? key.i == n : key.axis == 1 ? key.j == n : key.k == n);
        if (lower || upper) continue;
        grid.all_faces.push_back(kv.second);
        const std::size_t aoi = key.axis == 0 ? key.i - 1 : key.i;
        const std::size_t aoj = key.axis == 1 ? key.j - 1 : key.j;
        const std::size_t aok = key.axis == 2 ? key.k - 1 : key.k;
        const auto interior = [n](std::size_t a, std::size_t b, std::size_t c3) {
            return a >= 1 && a + 1 < n && b >= 1 && b + 1 < n && c3 >= 1 && c3 + 1 < n;
        };
        if (interior(aoi, aoj, aok) && interior(key.i, key.j, key.k))
            grid.deep_faces.push_back(kv.second);
    }
    return grid;
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(12);
        check_extremes();
        check_order();
        std::cout << "CONVECTION_BLENDING_VERIFICATION: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "CONVECTION_BLENDING_VERIFICATION: FAIL: " << e.what() << "\n";
        return 1;
    }
}