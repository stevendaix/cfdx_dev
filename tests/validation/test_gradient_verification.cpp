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

struct FieldCase {
    const char* name;
    double (*value)(const Vec3&);
    Vec3 (*exact_gradient)(const Vec3&);
};

double constant_value(const Vec3&) { return 3.0; }
Vec3 constant_gradient(const Vec3&) { return Vec3{0.0, 0.0, 0.0}; }

double linear_value(const Vec3& p) { return 2.0 * p.x - 3.0 * p.y + 0.5 * p.z; }
Vec3 linear_gradient(const Vec3&) { return Vec3{2.0, -3.0, 0.5}; }

double quadratic_value(const Vec3& p) { return p.x * p.x + p.y * p.y + p.z * p.z; }
Vec3 quadratic_gradient(const Vec3& p) { return Vec3{2.0 * p.x, 2.0 * p.y, 2.0 * p.z}; }

double smooth_value(const Vec3& p) {
    const double pi = std::acos(-1.0);
    return std::sin(pi * p.x) * std::cos(pi * p.y) * std::sin(pi * p.z);
}
Vec3 smooth_gradient(const Vec3& p) {
    const double pi = std::acos(-1.0);
    return Vec3{
        pi * std::cos(pi * p.x) * std::cos(pi * p.y) * std::sin(pi * p.z),
        -pi * std::sin(pi * p.x) * std::sin(pi * p.y) * std::sin(pi * p.z),
        pi * std::sin(pi * p.x) * std::cos(pi * p.y) * std::cos(pi * p.z)};
}

inline const FieldCase kConstant{"constant", constant_value, constant_gradient};
inline const FieldCase kLinear{"linear", linear_value, linear_gradient};
inline const FieldCase kQuadratic{"quadratic", quadratic_value, quadratic_gradient};
inline const FieldCase kSmooth{"smooth", smooth_value, smooth_gradient};

// Gradient error over the requested cells, volume-weighted. The error is the
// full vector difference |g_num - g_exact| (not a magnitude difference), so a
// first-order direction error cannot hide behind a magnitude comparison.
ErrorMetrics gradient_error(const Grid& grid, const Field<double, Location::CELL>& grad,
                            const FieldCase& field, bool interior_only)
{
    const std::size_t nc = grid.mesh.n_cells();
    std::vector<double> got, want, weights;
    const auto collect = [&](std::size_t c) {
        const Vec3 gc = grid.geometry.cell_centres[c];
        const Vec3 e = field.exact_gradient(gc);
        const Vec3 g{grad(c, 0), grad(c, 1), grad(c, 2)};
        const double err = (g - e).mag();
        const double magnitude = e.mag();
        got.push_back(magnitude + err);   // error_norms recovers err = got - want
        want.push_back(magnitude);
        weights.push_back(grid.geometry.cell_volumes[c]);
    };
    if (interior_only) {
        for (const std::size_t c : grid.interior) collect(c);
    } else {
        for (std::size_t c = 0; c < nc; ++c) collect(c);
    }
    return error_norms(got, want, weights);
}

Field<double, Location::CELL> sample_field(const Grid& grid, const FieldCase& field)
{
    Field<double, Location::CELL> f(grid.mesh.n_cells(), "phi", "1", 1);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c)
        f(c) = field.value(grid.geometry.cell_centres[c]);
    return f;
}

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

// The three gradient schemes exercised by the campaign: cell-based
// Green-Gauss, vertex-based (secondary) Green-Gauss and least-squares.
enum class GradScheme {
    GaussCell, GaussVertex, GaussPoint, LeastSquares,
    LeastSquaresWeighted, LeastSquaresWeightedInvR
};

std::string scheme_name(GradScheme s)
{
    switch (s) {
        case GradScheme::GaussCell:   return "green_gauss";
        case GradScheme::GaussVertex: return "green_gauss_vertex";
        case GradScheme::GaussPoint:  return "green_gauss_point";
        case GradScheme::LeastSquares: return "least_squares";
        case GradScheme::LeastSquaresWeighted: return "weighted_least_squares_1_over_r2";
        case GradScheme::LeastSquaresWeightedInvR: return "weighted_least_squares_1_over_r";
    }
    return "unknown";
}

Field<double, Location::CELL> compute_scheme(const Grid& grid,
                                             const Field<double, Location::CELL>& phi,
                                             GradScheme s)
{
    switch (s) {
        case GradScheme::GaussCell:
            return compute_gradient_gauss(phi, grid.mesh, grid.geometry);
        case GradScheme::GaussVertex:
            return compute_gradient_gauss_vertex(phi, grid.mesh, grid.geometry);
        case GradScheme::GaussPoint:
            return compute_gradient_gauss_point(phi, grid.mesh, grid.geometry);
        case GradScheme::LeastSquares:
            return compute_gradient_least_squares(phi, grid.mesh);
        case GradScheme::LeastSquaresWeighted:
            return compute_gradient_weighted_least_squares_extended(
                phi, grid.mesh, GradientWeighting::INVERSE_DISTANCE_SQUARED);
        case GradScheme::LeastSquaresWeightedInvR:
            return compute_gradient_weighted_least_squares_extended(
                phi, grid.mesh, GradientWeighting::INVERSE_DISTANCE);
    }
    return Field<double, Location::CELL>();
}

std::vector<double> refinement_errors(double shear, double stretch, const FieldCase& field,
                                      GradScheme scheme)
{
    std::vector<double> errors;
    for (std::size_t n : {4u, 8u, 16u}) {
        const Grid grid = make_affine_cube(n, shear, stretch);
        const auto phi = sample_field(grid, field);
        const auto grad = compute_scheme(grid, phi, scheme);
        const auto e = gradient_error(grid, grad, field, true);
        errors.push_back(e.l2);
    }
    return errors;
}

void report(const std::string& family, const std::string& scheme, const FieldCase& field,
            const ErrorMetrics& e)
{
    std::cout << "GRADIENT_VV family=" << family << " scheme=" << scheme
              << " field=" << field.name
              << " L1=" << e.l1 << " L2=" << e.l2 << " Linf=" << e.linf;
    if (!std::isfinite(e.l1_relative) || !std::isfinite(e.l2_relative))
        std::cout << " relative=nan";
    std::cout << "\n";
}

void check_exactness(const std::string& family, double shear, double stretch)
{
    const Grid grid = make_affine_cube(6, shear, stretch);

    for (const GradScheme scheme : {GradScheme::GaussCell, GradScheme::GaussVertex,
                                    GradScheme::GaussPoint, GradScheme::LeastSquares,
                                    GradScheme::LeastSquaresWeighted,
                                    GradScheme::LeastSquaresWeightedInvR}) {
        const std::string name = scheme_name(scheme);

        const auto constant_field = sample_field(grid, kConstant);
        const auto cg = compute_scheme(grid, constant_field, scheme);
        const auto ce = gradient_error(grid, cg, kConstant, false);
        report(family, name, kConstant, ce);
        require(ce.linf <= 1e-12,
                family + "/" + name + ": constant field must have zero gradient");

        const auto linear_field = sample_field(grid, kLinear);
        const auto lg = compute_scheme(grid, linear_field, scheme);
        const auto le = gradient_error(grid, lg, kLinear, true);
        report(family, name, kLinear, le);
        require(le.linf <= 1e-9,
                family + "/" + name + ": linear field must be exact on the interior");
    }
}

void check_order(const std::string& family, double shear, double stretch, const FieldCase& field)
{
    for (const GradScheme scheme : {GradScheme::GaussCell, GradScheme::GaussVertex,
                                    GradScheme::GaussPoint, GradScheme::LeastSquares,
                                    GradScheme::LeastSquaresWeighted,
                                    GradScheme::LeastSquaresWeightedInvR}) {
        const std::string name = scheme_name(scheme);
        const auto errors = refinement_errors(shear, stretch, field, scheme);
        for (std::size_t i = 0; i < errors.size(); ++i) {
            std::cout << "GRADIENT_VV_ORDER family=" << family << " scheme=" << name
                      << " field=" << field.name << " n=" << (4u << i)
                      << " L2=" << errors[i];
            if (i > 0)
                std::cout << " order=" << observed_order(errors[i - 1], errors[i]);
            std::cout << "\n";
        }
        require(errors.back() < errors.front(),
                family + "/" + name + "/" + field.name + ": refinement must reduce the error");
        require_order(errors, 1.0, 1.0, family + "/" + name + "/" + field.name);
    }
}

// Green-Gauss and least-squares are both exact for a quadratic field on an
// affine (uniform-lattice) mesh. For Green-Gauss the second-order
// face-interpolation offsets cancel between opposite faces. For least-squares
// the |Delta|^2 term in Delta(phi) is even in Delta and cancels over a
// centrally-symmetric +/- stencil, so the normal equations return the exact
// gradient. This is a stronger property than first-order exactness and is part
// of the reconstruction contract. The vertex-based Green-Gauss scheme is not
// quadratic-exact in general, so it is reported here and its reduction/order
// behaviour is covered by check_order.
void check_quadratic(const std::string& family, double shear, double stretch)
{
    const Grid grid = make_affine_cube(8, shear, stretch);

    const auto field = sample_field(grid, kQuadratic);

    const auto gg = compute_gradient_gauss(field, grid.mesh, grid.geometry);
    const auto gg_err = gradient_error(grid, gg, kQuadratic, true);
    report(family, scheme_name(GradScheme::GaussCell), kQuadratic, gg_err);
    require(gg_err.linf <= 1e-9,
            family + ": Green-Gauss must be quadratic-exact on an affine mesh");

    const auto vg = compute_gradient_gauss_vertex(field, grid.mesh, grid.geometry);
    const auto vg_err = gradient_error(grid, vg, kQuadratic, true);
    report(family, scheme_name(GradScheme::GaussVertex), kQuadratic, vg_err);
    require(std::isfinite(vg_err.l1) && std::isfinite(vg_err.l2) &&
            std::isfinite(vg_err.linf),
            family + ": vertex Green-Gauss quadratic error must be finite");

    for (const GradScheme scheme : {GradScheme::LeastSquares,
                                    GradScheme::LeastSquaresWeighted,
                                    GradScheme::LeastSquaresWeightedInvR}) {
        const auto ls = compute_scheme(grid, field, scheme);
        const auto ls_err = gradient_error(grid, ls, kQuadratic, true);
        report(family, scheme_name(scheme), kQuadratic, ls_err);
        require(ls_err.linf <= 1e-9,
                family + "/" + scheme_name(scheme) +
                ": least-squares family must be quadratic-exact on an affine mesh");
    }
}

// Boundary-neighbour policy: the cell-based Green-Gauss operator treats a
// boundary face by the zero-gradient rule (phi_f = phi_owner), represented by
// a zero geometric interpolation weight. This contract is what makes interior
// gradients exact while boundary-adjacent cells carry the boundary
// approximation; pin it so it cannot silently change.
void check_boundary_policy(const std::string& family, double shear, double stretch)
{
    const Grid grid = make_affine_cube(6, shear, stretch);

    const FaceOwnership& own = grid.mesh.ownership();
    std::size_t boundary_face = 0;
    for (std::size_t f = 0; f < grid.mesh.n_faces(); ++f) {
        if (own.neighbour(f) < 0) { boundary_face = f; break; }
    }
    const std::size_t owner = own.owner(boundary_face);
    const double w = geometric_interpolation_weight(
        grid.geometry.face_centres[boundary_face],
        grid.geometry.cell_centres[owner], nullptr);
    require(w == 0.0, family + ": boundary face interpolation weight must be 0 (zero-gradient)");

    // The boundary approximation pollutes boundary-adjacent cells: their
    // cell-based gradient must differ from the exact linear gradient.
    const auto phi = sample_field(grid, kLinear);
    const auto grad = compute_gradient_gauss(phi, grid.mesh, grid.geometry);
    const Vec3 gc = grid.geometry.cell_centres[owner];
    const Vec3 exact = kLinear.exact_gradient(gc);
    const Vec3 g{grad(owner, 0), grad(owner, 1), grad(owner, 2)};
    const double err = (g - exact).mag();
    require(err > 1e-9,
            family + ": boundary-adjacent gradient must carry the zero-gradient boundary approximation");
    std::cout << "GRADIENT_BOUNDARY family=" << family
              << " scheme=green_gauss boundary_cell_err=" << err << "\n";
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(12);

        check_exactness("orthogonal", 0.0, 1.0);
        check_exactness("sheared", 0.5, 1.0);
        check_exactness("stretched", 0.0, 4.0);

        check_boundary_policy("orthogonal", 0.0, 1.0);
        check_boundary_policy("sheared", 0.5, 1.0);
        check_boundary_policy("stretched", 0.0, 4.0);

        check_quadratic("orthogonal", 0.0, 1.0);
        check_quadratic("sheared", 0.5, 1.0);
        check_quadratic("stretched", 0.0, 4.0);

        check_order("orthogonal", 0.0, 1.0, kSmooth);
        check_order("sheared", 0.5, 1.0, kSmooth);
        check_order("stretched", 0.0, 4.0, kSmooth);

        std::cout << "GRADIENT_VERIFICATION: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "GRADIENT_VERIFICATION: FAIL: " << e.what() << "\n";
        return 1;
    }
}
