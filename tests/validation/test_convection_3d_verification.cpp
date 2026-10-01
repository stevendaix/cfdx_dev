// Issue #461 — N4 multidimensional convection scheme verification.
//
// The 1-D slab campaign (test_convection_scheme_verification.cpp) qualified
// the TVD/MUSCL face reconstructions; this campaign qualifies them on a
// genuinely 3-D Cartesian mesh:
//   - boundedness: every TVD limiter keeps face values in the adjacent-cell
//     envelope on a 3-D step;
//   - linear exactness: every TVD limiter reproduces a linear field exactly on
//     deep interior faces (psi(1) = 1 with an exact Green-Gauss gradient);
//   - observed order: first-order upwind vs the TVD limiters on a smooth field;
//   - constant advection: the discrete divergence of a constant field is zero.

#include "cfdx/core/numerics/convection.h"
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
    std::vector<std::size_t> deep_faces;   // internal faces with interior owner/neighbour
    std::vector<std::size_t> all_faces;    // every internal face
};

// Cartesian n x n x n hexahedral grid on [0,1]^3 (affine cube builder).
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
    struct FaceKey {
        std::size_t i, j, k, axis;
        bool operator==(const FaceKey& o) const {
            return i == o.i && j == o.j && k == o.k && axis == o.axis;
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

    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t c = cell_id(i, j, k);
                auto face = [&](std::size_t a, std::size_t fi, std::size_t fj, std::size_t fk,
                                std::initializer_list<std::size_t> v) {
                    return add_or_get_face({fi, fj, fk, a}, v);
                };
                faces_per_cell[c] = {
                    face(0, i, j, k, {point(c,0),point(c,4),point(c,7),point(c,3)}),
                    face(0, i + 1, j, k, {point(c,1),point(c,2),point(c,6),point(c,5)}),
                    face(1, i, j, k, {point(c,0),point(c,1),point(c,5),point(c,4)}),
                    face(1, i, j + 1, k, {point(c,3),point(c,7),point(c,6),point(c,2)}),
                    face(2, i, j, k, {point(c,0),point(c,3),point(c,2),point(c,1)}),
                    face(2, i, j, k + 1, {point(c,4),point(c,5),point(c,6),point(c,7)})};
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
        const std::size_t owner = cell_id(oi, oj, ok);
        m.ownership().set_owner(f, owner);
        if (lower || upper) {
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        } else {
            std::size_t ni = oi, nj = oj, nk = ok;
            if (key.axis == 0) ++ni; else if (key.axis == 1) ++nj; else ++nk;
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
    for (const auto& kv : face_map) {
        const FaceKey key = kv.first;
        const bool lower = (key.axis == 0 ? key.i == 0 : key.axis == 1 ? key.j == 0 : key.k == 0);
        const bool upper = (key.axis == 0 ? key.i == n : key.axis == 1 ? key.j == n : key.k == n);
        if (lower || upper) continue;
        grid.all_faces.push_back(kv.second);
        // Internal face on plane (i,j,k) separates owner cell (i-1,j,k in its
        // axis) and neighbour cell (i,j,k in its axis); deep means both CELLS
        // are interior (indices in [1, n-2]).
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

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

Field<double, Location::CELL> cell_field(const Grid& grid, std::vector<double> values)
{
    Field<double, Location::CELL> f(grid.mesh.n_cells(), "phi", "1", 1);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c) f(c) = values[c];
    return f;
}

Field<double, Location::FACE> positive_flux(const Grid& grid, double value)
{
    Field<double, Location::FACE> flux(grid.mesh.n_faces(), "F", "m2/s", 1);
    for (std::size_t f = 0; f < grid.mesh.n_faces(); ++f) flux(f) = value;
    for (std::size_t f = 0; f < grid.mesh.n_faces(); ++f)
        if (grid.mesh.ownership().neighbour(f) < 0) flux(f) = 0.0;   // walls
    return flux;
}

std::vector<double> sample(const Grid& grid, double (*fn)(const Vec3&)) {
    std::vector<double> v(grid.mesh.n_cells(), 0.0);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c)
        v[c] = fn(grid.geometry.cell_centres[c]);
    return v;
}

// Faces whose centroid lies in a FIXED nested physical window, so the refined
// order measurement is converged on the same region.
std::vector<std::size_t> window_faces(const Grid& grid)
{
    std::vector<std::size_t> out;
    for (const std::size_t f : grid.all_faces) {
        const Vec3& c = grid.geometry.face_centres[f];
        if (c.x >= 0.25 && c.x <= 0.75 && c.y >= 0.25 && c.y <= 0.75 &&
            c.z >= 0.25 && c.z <= 0.75)
            out.push_back(f);
    }
    return out;
}

ErrorMetrics face_error(const Grid& grid, const Field<double, Location::FACE>& face,
                        double (*fn)(const Vec3&), bool deep)
{
    std::vector<double> got, want;
    const auto& faces = deep ? grid.deep_faces : grid.all_faces;
    for (const std::size_t f : faces) {
        got.push_back(face(f));
        want.push_back(fn(grid.geometry.face_centres[f]));
    }
    return error_norms(got, want);
}

double step_value(const Vec3& p) { return (p.x < 0.5) ? 5.0 : -2.0; }
double linear_value(const Vec3& p) { return 2.0 * p.x - 3.0 * p.y + 0.5 * p.z + 1.0; }
double smooth_value(const Vec3& p) {
    // Monotone smooth profile with a non-degenerate gradient in every axis, so
    // first-order (upwind) and second-order (TVD in smooth regions) face
    // reconstruction errors are cleanly resolved.
    return p.x + p.y + p.z
         + 0.5 * (p.x * p.x + p.y * p.y + p.z * p.z);
}

void check_boundedness_linear()
{
    const Grid grid = make_cartesian_cube(12);
    const auto flux = positive_flux(grid, 1.0);
    const auto phi = cell_field(grid, sample(grid, step_value));
    const auto grad = compute_gradient_gauss(phi, grid.mesh, grid.geometry);

    const std::vector<LimiterType> limiters = {
        LimiterType::MINMOD, LimiterType::VANLEER, LimiterType::SUPERBEE,
        LimiterType::VAN_ALBADA, LimiterType::MC};
    for (const LimiterType limiter : limiters) {
        const auto face = interpolate_cell_to_face(
            phi, grid.mesh, grid.geometry, InterpScheme::LIMITED, &flux, limiter, &grad);
        const auto& own = grid.mesh.ownership();
        for (const std::size_t f : grid.all_faces) {
            const std::size_t o = own.owner(f);
            const std::size_t nb = static_cast<std::size_t>(own.neighbour(f));
            const double fv = face(f);
            require(fv >= std::min(phi(o), phi(nb)) - 1e-12 &&
                    fv <= std::max(phi(o), phi(nb)) + 1e-12,
                    std::string(to_string(limiter)) + ": bounded in 3-D");
        }
        std::cout << "CONV3D_BOUNDED limiter=" << to_string(limiter) << " ok\n";
    }

    // Linear exactness on deep faces for all limiters.
    const auto lin = cell_field(grid, sample(grid, linear_value));
    const auto lgrad = compute_gradient_gauss(lin, grid.mesh, grid.geometry);
    for (const LimiterType limiter : limiters) {
        const auto face = interpolate_cell_to_face(
            lin, grid.mesh, grid.geometry, InterpScheme::LIMITED, &flux, limiter, &lgrad);
        const auto e = face_error(grid, face, linear_value, true);
        require(e.linf <= 1e-9, std::string(to_string(limiter)) + ": linear-exact in 3-D");
        std::cout << "CONV3D_LINEAR limiter=" << to_string(limiter) << " Linf=" << e.linf << "\n";
    }
}

void check_order()
{
    const std::vector<std::size_t> ns = {8u, 16u, 32u};
    auto measure = [](const Grid& grid, const Field<double, Location::FACE>& face) {
        const auto win = window_faces(grid);
        std::vector<double> got, want;
        for (const std::size_t f : win) {
            got.push_back(face(f));
            want.push_back(smooth_value(grid.geometry.face_centres[f]));
        }
        return error_norms(got, want).l2;
    };

    // First-order upwind.
    {
        std::vector<double> errors;
        for (const std::size_t n : ns) {
            const Grid grid = make_cartesian_cube(n);
            const auto flux = positive_flux(grid, 1.0);
            const auto phi = cell_field(grid, sample(grid, smooth_value));
            const auto face = interpolate_cell_to_face(phi, grid.mesh, grid.geometry, InterpScheme::UPWIND, &flux);
            errors.push_back(measure(grid, face));
        }
        std::cout << "CONV3D_ORDER scheme=upwind";
        for (std::size_t i = 0; i < errors.size(); ++i) {
            std::cout << " n=" << ns[i] << " L2=" << errors[i];
            if (i) std::cout << " order=" << observed_order(errors[i - 1], errors[i]);
        }
        std::cout << "\n";
        require_order(errors, 1.0, 0.8, "3-D upwind");
    }

    for (const LimiterType limiter : {LimiterType::MINMOD, LimiterType::VANLEER, LimiterType::MC}) {
        std::vector<double> errors;
        for (const std::size_t n : ns) {
            const Grid grid = make_cartesian_cube(n);
            const auto flux = positive_flux(grid, 1.0);
            const auto phi = cell_field(grid, sample(grid, smooth_value));
            const auto grad = compute_gradient_gauss(phi, grid.mesh, grid.geometry);
            const auto face = interpolate_cell_to_face(
                phi, grid.mesh, grid.geometry, InterpScheme::LIMITED, &flux, limiter, &grad);
            errors.push_back(measure(grid, face));
        }
        std::cout << "CONV3D_ORDER scheme=tvd_" << to_string(limiter);
        for (std::size_t i = 0; i < errors.size(); ++i) {
            std::cout << " n=" << ns[i] << " L2=" << errors[i];
            if (i) std::cout << " order=" << observed_order(errors[i - 1], errors[i]);
        }
        std::cout << "\n";
        require_order(errors, 1.0, 0.8, std::string("3-D tvd_") + to_string(limiter));
    }
}


void check_limiter_gradient_variants()
{
    // N2/N4 bridge: the limiter contract must remain bounded and convergent
    // when its reconstruction gradient is supplied by the production
    // least-squares family, not only by Green-Gauss.
    const std::vector<GradientScheme> gradients = {
        GradientScheme::LEAST_SQUARES,
        GradientScheme::WEIGHTED_LEAST_SQUARES,
        GradientScheme::LEAST_SQUARES_QUADRATIC};

    const std::vector<LimiterType> limiters = {
        LimiterType::MINMOD, LimiterType::VANLEER, LimiterType::SUPERBEE,
        LimiterType::VAN_ALBADA, LimiterType::MC};

    const Grid grid = make_cartesian_cube(12);
    const auto flux = positive_flux(grid, 1.0);
    const auto phi = cell_field(grid, sample(grid, step_value));

    for (const GradientScheme gs : gradients) {
        const auto grad = cell_gradient(phi, grid.mesh, grid.geometry, gs);
        for (const LimiterType limiter : limiters) {
            const auto face = interpolate_cell_to_face(
                phi, grid.mesh, grid.geometry, InterpScheme::LIMITED,
                &flux, limiter, &grad);
            const auto& own = grid.mesh.ownership();
            for (const std::size_t f : grid.all_faces) {
                const std::size_t o = own.owner(f);
                const std::size_t nb = static_cast<std::size_t>(own.neighbour(f));
                require(face(f) >= std::min(phi(o), phi(nb)) - 1e-12 &&
                        face(f) <= std::max(phi(o), phi(nb)) + 1e-12,
                        std::string(to_string(gs)) + "/" + to_string(limiter) +
                        ": WLS-family limiter reconstruction must remain bounded");
            }
        }
        std::cout << "CONV3D_WLS_BOUNDED gradient=" << to_string(gs) << " ok\n";
    }

    // Linear exactness is checked independently of the limiter choice. On
    // this affine Cartesian stencil the LS family should reproduce the exact
    // cell gradient, hence psi(1)=1 must reconstruct the exact face value.
    const auto linear = cell_field(grid, sample(grid, linear_value));
    for (const GradientScheme gs : gradients) {
        const auto grad = cell_gradient(linear, grid.mesh, grid.geometry, gs);
        const auto face = interpolate_cell_to_face(
            linear, grid.mesh, grid.geometry, InterpScheme::LIMITED,
            &flux, LimiterType::VANLEER, &grad);
        const auto e = face_error(grid, face, linear_value, true);
        require(e.linf <= 1e-9,
                std::string(to_string(gs)) + ": limiter + LS-family gradient must be linear-exact");
        std::cout << "CONV3D_WLS_LINEAR gradient=" << to_string(gs)
                  << " Linf=" << e.linf << " L2=" << e.l2 << "\n";
    }

    // Smooth-field refinement: require a genuinely second-order TVD path for
    // each LS-family gradient. No isolated one-level pass is sufficient.
    for (const GradientScheme gs : gradients) {
        std::vector<double> errors;
        for (const std::size_t n : {8u, 16u, 32u}) {
            const Grid g = make_cartesian_cube(n);
            const auto f = positive_flux(g, 1.0);
            const auto field = cell_field(g, sample(g, smooth_value));
            const auto grad = cell_gradient(field, g.mesh, g.geometry, gs);
            const auto face = interpolate_cell_to_face(
                field, g.mesh, g.geometry, InterpScheme::LIMITED,
                &f, LimiterType::VANLEER, &grad);
            std::vector<double> got, want;
            for (const std::size_t face_id : window_faces(g)) {
                got.push_back(face(face_id));
                want.push_back(smooth_value(g.geometry.face_centres[face_id]));
            }
            errors.push_back(error_norms(got, want).l2);
        }
        std::cout << "CONV3D_WLS_ORDER gradient=" << to_string(gs);
        for (std::size_t k = 0; k < errors.size(); ++k) {
            std::cout << " n=" << {8u,16u,32u}[k] << " L2=" << errors[k];
            if (k) std::cout << " order=" << observed_order(errors[k-1], errors[k]);
        }
        std::cout << "\n";
        require_order(errors, 2.0, 1.5,
                      std::string("limiter + ") + to_string(gs) + " smooth-field refinement");
    }
}

void check_constant_advection()
{
    const Grid grid = make_cartesian_cube(10);
    const auto flux = positive_flux(grid, 1.0);
    const auto phi = cell_field(grid, std::vector<double>(grid.mesh.n_cells(), 5.0));
    const auto conv = compute_convection(phi, flux, grid.mesh, InterpScheme::UPWIND);
    double mx = 0.0;
    for (const std::size_t c : grid.interior) mx = std::max(mx, std::abs(conv(c)));
    std::cout << "CONV3D_CONSTANT max_interior_div=" << mx << "\n";
    require(mx <= 1e-12, "3-D constant advection must have zero divergence in interior cells");
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(12);
        check_boundedness_linear();
        check_order();
        check_limiter_gradient_variants();
        check_constant_advection();
        std::cout << "CONVECTION_3D_VERIFICATION: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "CONVECTION_3D_VERIFICATION: FAIL: " << e.what() << "\n";
        return 1;
    }
}