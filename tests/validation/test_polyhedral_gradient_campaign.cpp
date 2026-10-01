// Issue #461 — N2 polyhedral (tetrahedral) gradient campaign.
//
// The gradient V&V matrix so far covered hexahedral meshes (orthogonal,
// affine-sheared, stretched). This campaign verifies the gradient operators on
// an unstructured tetrahedral mesh: a `n x n x n` cube grid triangulated into
// 6 Kuhn tets per cube with a globally consistent diagonal pattern.
//
// Assertions:
//   - constant field -> exact zero gradient for Green-Gauss, vertex Green-Gauss
//     and least-squares;
//   - linear field    -> exact interior gradient for all three schemes;
//   - smooth field    -> measured observed order on interior tets.

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
#include <utility>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::verification;

namespace {

struct Grid {
    Mesh mesh;
    GeometryCache geometry;
    std::vector<std::size_t> interior;   // tets whose cube is not on the boundary
    std::vector<std::size_t> interior_vertices;
};

// Local cube corners, same orientation for every cube so the Kuhn diagonal
// pattern matches across shared faces.
//   L0=(0,0,0) L1=(1,0,0) L2=(1,1,0) L3=(0,1,0)
//   L4=(0,0,1) L5=(1,0,1) L6=(1,1,1) L7=(0,1,1)
static const std::array<std::array<unsigned, 3>, 8> kCorner = {{
    {{0,0,0}}, {{1,0,0}}, {{1,1,0}}, {{0,1,0}},
    {{0,0,1}}, {{1,0,1}}, {{1,1,1}}, {{0,1,1}}
}};
// The 6 Kuhn tets per cube, all sharing the main diagonal L0-L6.
static const std::array<std::array<unsigned, 4>, 6> kTets = {{
    {{0,1,2,6}}, {{0,1,5,6}}, {{0,3,2,6}},
    {{0,3,7,6}}, {{0,4,5,6}}, {{0,4,7,6}}
}};

std::string tri_key(std::size_t a, std::size_t b, std::size_t c)
{
    std::array<std::size_t, 3> v{a, b, c};
    std::sort(v.begin(), v.end());
    return std::to_string(v[0]) + "/" + std::to_string(v[1]) + "/" + std::to_string(v[2]);
}

Grid make_tet_grid(std::size_t n)
{
    if (n < 1) throw std::invalid_argument("make_tet_grid: n must be >= 1");

    const std::size_t nv = n + 1;
    const auto vid = [nv](std::size_t x, std::size_t y, std::size_t z) {
        return x + nv * (y + nv * z);
    };

    Mesh m;
    m.points().resize(nv * nv * nv);
    for (std::size_t z = 0; z <= n; ++z)
        for (std::size_t y = 0; y <= n; ++y)
            for (std::size_t x = 0; x <= n; ++x)
                m.points().set(vid(x, y, z),
                    static_cast<double>(x) / n, static_cast<double>(y) / n,
                    static_cast<double>(z) / n);

    std::unordered_map<std::string, std::size_t> face_ids;
    std::vector<std::vector<std::size_t>> face_cells;
    std::vector<std::vector<std::size_t>> cell_faces;   // per tet, the 4 face ids

    auto get_face = [&](std::size_t a, std::size_t b, std::size_t c) {
        const std::string key = tri_key(a, b, c);
        auto it = face_ids.find(key);
        if (it != face_ids.end()) return it->second;
        const std::size_t id = face_cells.size();
        m.faces().push_face({a, b, c});
        face_ids.emplace(key, id);
        face_cells.emplace_back();
        return id;
    };

    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                std::array<std::size_t, 8> corner;
                for (unsigned q = 0; q < 8; ++q)
                    corner[q] = vid(i + kCorner[q][0], j + kCorner[q][1], k + kCorner[q][2]);

                for (unsigned t = 0; t < kTets.size(); ++t) {
                    const auto& tet = kTets[t];
                    std::array<std::size_t, 4> v{
                        corner[tet[0]], corner[tet[1]], corner[tet[2]], corner[tet[3]]};
                    std::vector<std::size_t> faces(4);
                    faces[0] = get_face(v[0], v[1], v[2]);
                    faces[1] = get_face(v[0], v[1], v[3]);
                    faces[2] = get_face(v[0], v[2], v[3]);
                    faces[3] = get_face(v[1], v[2], v[3]);
                    const std::size_t cell = cell_faces.size();
                    for (const std::size_t f : faces) face_cells[f].push_back(cell);
                    cell_faces.push_back(faces);
                }
            }
        }
    }

    m.ownership().resize(m.n_faces());
    for (std::size_t f = 0; f < m.n_faces(); ++f) {
        const auto& cells = face_cells[f];
        if (cells.empty()) throw std::runtime_error("make_tet_grid: face with no owner");
        m.ownership().set_owner(f, cells[0]);
        if (cells.size() == 1) {
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        } else {
            if (cells.size() != 2)
                throw std::runtime_error("make_tet_grid: face shared by more than two tets");
            m.ownership().set_neighbour(f, static_cast<std::int64_t>(cells[1]));
        }
    }

    for (const auto& faces : cell_faces)
        m.cells().push_cell(faces);

    Grid grid;
    grid.mesh = std::move(m);
    grid.geometry = make_geometry_cache(grid.mesh);

    // Interior = tets whose cube is not adjacent to the domain boundary.
    for (std::size_t k = 1; k + 1 < n; ++k)
        for (std::size_t j = 1; j + 1 < n; ++j)
            for (std::size_t i = 1; i + 1 < n; ++i) {
                const std::size_t base = 6 * (i + n * (j + n * k));
                for (unsigned t = 0; t < kTets.size(); ++t)
                    grid.interior.push_back(base + t);
            }

    for (std::size_t z = 1; z < n; ++z)
        for (std::size_t y = 1; y < n; ++y)
            for (std::size_t x = 1; x < n; ++x)
                grid.interior_vertices.push_back(vid(x, y, z));
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

double quadratic_value(const Vec3& p) { return p.x*p.x + 2.0*p.x*p.y + 0.5*p.y*p.y + 0.75*p.z*p.z; }
Vec3 quadratic_gradient(const Vec3& p) { return Vec3{2.0*p.x + 2.0*p.y, 2.0*p.x + p.y, 1.5*p.z}; }

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

double cubic_x3(const Vec3& p) { return p.x*p.x*p.x; }
Vec3 cubic_x3_gradient(const Vec3& p) { return Vec3{3.0*p.x*p.x, 0.0, 0.0}; }
double cubic_y3(const Vec3& p) { return p.y*p.y*p.y; }
Vec3 cubic_y3_gradient(const Vec3& p) { return Vec3{0.0, 3.0*p.y*p.y, 0.0}; }
double cubic_z3(const Vec3& p) { return p.z*p.z*p.z; }
Vec3 cubic_z3_gradient(const Vec3& p) { return Vec3{0.0, 0.0, 3.0*p.z*p.z}; }
double cubic_x2y(const Vec3& p) { return p.x*p.x*p.y; }
Vec3 cubic_x2y_gradient(const Vec3& p) { return Vec3{2.0*p.x*p.y, p.x*p.x, 0.0}; }
double cubic_x2z(const Vec3& p) { return p.x*p.x*p.z; }
Vec3 cubic_x2z_gradient(const Vec3& p) { return Vec3{2.0*p.x*p.z, 0.0, p.x*p.x}; }
double cubic_y2x(const Vec3& p) { return p.y*p.y*p.x; }
Vec3 cubic_y2x_gradient(const Vec3& p) { return Vec3{p.y*p.y, 2.0*p.x*p.y, 0.0}; }
double cubic_y2z(const Vec3& p) { return p.y*p.y*p.z; }
Vec3 cubic_y2z_gradient(const Vec3& p) { return Vec3{0.0, 2.0*p.y*p.z, p.y*p.y}; }
double cubic_z2x(const Vec3& p) { return p.z*p.z*p.x; }
Vec3 cubic_z2x_gradient(const Vec3& p) { return Vec3{p.z*p.z, 0.0, 2.0*p.x*p.z}; }
double cubic_z2y(const Vec3& p) { return p.z*p.z*p.y; }
Vec3 cubic_z2y_gradient(const Vec3& p) { return Vec3{0.0, p.z*p.z, 2.0*p.y*p.z}; }
double cubic_xyz(const Vec3& p) { return p.x*p.y*p.z; }
Vec3 cubic_xyz_gradient(const Vec3& p) { return Vec3{p.y*p.z, p.x*p.z, p.x*p.y}; }

enum class GradScheme { GaussCell, GaussVertex, GaussPoint, LeastSquares, LeastSquares2, LeastSquaresQuad };
std::string scheme_name(GradScheme s)
{
    switch (s) {
        case GradScheme::GaussCell:   return "green_gauss";
        case GradScheme::GaussVertex: return "green_gauss_vertex";
        case GradScheme::GaussPoint:  return "green_gauss_point";
        case GradScheme::LeastSquares: return "least_squares";
        case GradScheme::LeastSquares2: return "least_squares_2ring";
        case GradScheme::LeastSquaresQuad: return "least_squares_quadratic";
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
        case GradScheme::LeastSquares2:
            return compute_gradient_least_squares_extended(phi, grid.mesh);
        case GradScheme::LeastSquaresQuad:
            return compute_gradient_least_squares_quadratic(phi, grid.mesh);
    }
    return Field<double, Location::CELL>();
}

ErrorMetrics gradient_error(const Grid& grid, const Field<double, Location::CELL>& grad,
                            const FieldCase& field, bool interior_only)
{
    std::vector<double> got, want, weights;
    std::vector<std::size_t> cells;
    if (interior_only) {
        cells = grid.interior;
    } else {
        cells.resize(grid.mesh.n_cells());
        for (std::size_t c = 0; c < cells.size(); ++c) cells[c] = c;
    }
    for (const std::size_t c : cells) {
        const Vec3 gc = grid.geometry.cell_centres[c];
        const Vec3 e = field.exact_gradient(gc);
        const Vec3 g{grad(c, 0), grad(c, 1), grad(c, 2)};
        const double err = (g - e).mag();
        const double mag = e.mag();
        got.push_back(mag + err);
        want.push_back(mag);
        weights.push_back(grid.geometry.cell_volumes[c]);
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

void check_exactness(std::size_t n)
{
    const Grid grid = make_tet_grid(n);
    for (const GradScheme scheme : {GradScheme::GaussCell, GradScheme::GaussVertex,
                                    GradScheme::GaussPoint, GradScheme::LeastSquares, GradScheme::LeastSquares2,
                                    GradScheme::LeastSquaresQuad}) {
        const std::string name = scheme_name(scheme);

        const auto cf = sample_field(grid, kConstant);
        const auto ce = gradient_error(grid, compute_scheme(grid, cf, scheme), kConstant, false);
        std::cout << "POLY_GRAD n=" << n << " scheme=" << name << " field=constant"
                  << " L2=" << ce.l2 << " Linf=" << ce.linf << "\n";
        require(ce.linf <= 1e-12, name + ": constant field must have zero gradient");

        const auto lf = sample_field(grid, kLinear);
        const auto le = gradient_error(grid, compute_scheme(grid, lf, scheme), kLinear, true);
        std::cout << "POLY_GRAD n=" << n << " scheme=" << name << " field=linear"
                  << " L2=" << le.l2 << " Linf=" << le.linf << "\n";
        // Least-squares, vertex GG and point-linear GG are linear-exact on
        // interior tetrahedra; the two-point Gauss scheme is not (face
        // centroid off the owner-neighbour line).
        if (scheme != GradScheme::GaussCell)
            require(le.linf <= 1e-9, name + ": must be linear-exact on interior tets");
        else
            require(std::isfinite(le.l1) && std::isfinite(le.l2),
                    name + ": linear gradient error must be finite");
    }
}

struct StencilDiagnostics {
    std::size_t size = 0;
    int rank = 0;
    double sigma_min = 0.0;
    double sigma_max = 0.0;
    double condition = 0.0;
    double radius_min = 0.0;
    double radius_max = 0.0;
};

StencilDiagnostics quadratic_stencil_diagnostics(const Grid& grid, std::size_t c)
{
    const std::size_t n_cells = grid.mesh.n_cells();
    const auto* cell_faces = grid.mesh.cells().faces_data();
    const auto* cell_offsets = grid.mesh.cells().offsets_data();
    const FaceOwnership& own = grid.mesh.ownership();
    std::vector<std::vector<std::size_t>> ring1(n_cells);
    for (std::size_t cell = 0; cell < n_cells; ++cell) {
        auto& ring = ring1[cell];
        for (Offset k = cell_offsets[cell]; k < cell_offsets[cell + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            const std::size_t owner = own.owner(f);
            const std::int64_t nraw = own.neighbour(f);
            std::size_t nb = n_cells;
            if (owner == cell) {
                if (nraw >= 0) nb = static_cast<std::size_t>(nraw);
            } else {
                nb = owner;
            }
            if (nb < n_cells && nb != cell) ring.push_back(nb);
        }
        std::sort(ring.begin(), ring.end());
        ring.erase(std::unique(ring.begin(), ring.end()), ring.end());
    }

    std::vector<unsigned char> seen(n_cells, 0);
    std::vector<std::size_t> touched;
    std::vector<std::size_t> frontier = ring1[c];
    std::vector<std::size_t> stencil;
    auto add = [&](std::size_t nb) {
        if (nb == c || seen[nb]) return;
        seen[nb] = 1;
        touched.push_back(nb);
        stencil.push_back(nb);
    };
    for (const std::size_t nb : frontier) add(nb);
    for (int depth = 0; depth < 2; ++depth) {
        std::vector<std::size_t> next;
        for (const std::size_t nb : frontier)
            for (const std::size_t nb2 : ring1[nb])
                if (nb2 != c && !seen[nb2]) { add(nb2); next.push_back(nb2); }
        frontier.swap(next);
        if (frontier.empty()) break;
    }

    double h = 0.0;
    double rmin = std::numeric_limits<double>::max();
    double rmax = 0.0;
    for (const std::size_t nb : stencil) {
        const double r = (grid.geometry.cell_centres[nb] - grid.geometry.cell_centres[c]).mag();
        h += r;
        rmin = std::min(rmin, r);
        rmax = std::max(rmax, r);
    }
    h /= static_cast<double>(stencil.size());

    double A[9][9] = {};
    for (const std::size_t nb : stencil) {
        const Vec3 q = (grid.geometry.cell_centres[nb] - grid.geometry.cell_centres[c]) * (1.0 / h);
        const double r2 = q.mag2();
        const double w = 1.0 / r2;
        const double b[9] = {q.x,q.y,q.z,q.x*q.x,q.y*q.y,q.z*q.z,q.x*q.y,q.x*q.z,q.y*q.z};
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j)
                A[i][j] += w*b[i]*b[j];
    }

    // Jacobi eigensolve of the symmetric positive-semidefinite normal matrix.
    double D[9][9];
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            D[i][j] = A[i][j];
    for (int sweep = 0; sweep < 100; ++sweep) {
        int p = 0, q = 1;
        double largest = 0.0;
        for (int i = 0; i < 9; ++i)
            for (int j = i + 1; j < 9; ++j)
                if (std::abs(D[i][j]) > largest) { largest = std::abs(D[i][j]); p = i; q = j; }
        if (largest <= 1e-13) break;
        const double tau = (D[q][q] - D[p][p]) / (2.0 * D[p][q]);
        const double t = (tau >= 0.0 ? 1.0 : -1.0) /
                         (std::abs(tau) + std::sqrt(1.0 + tau*tau));
        const double cs = 1.0 / std::sqrt(1.0 + t*t);
        const double sn = t * cs;
        const double app = D[p][p], aqq = D[q][q], apq = D[p][q];
        D[p][p] = app - t*apq;
        D[q][q] = aqq + t*apq;
        D[p][q] = D[q][p] = 0.0;
        for (int k = 0; k < 9; ++k) {
            if (k == p || k == q) continue;
            const double dkp = D[k][p], dkq = D[k][q];
            D[k][p] = D[p][k] = cs*dkp - sn*dkq;
            D[k][q] = D[q][k] = sn*dkp + cs*dkq;
        }
    }

    std::array<double,9> lambda{};
    for (int i = 0; i < 9; ++i) lambda[i] = std::max(0.0, D[i][i]);
    std::sort(lambda.begin(), lambda.end());
    const double scale = lambda.back();
    const double tol = 4096.0 * std::numeric_limits<double>::epsilon() * scale;
    int rank = 0;
    for (const double v : lambda) if (v > tol) ++rank;

    StencilDiagnostics out;
    out.size = stencil.size();
    out.rank = rank;
    out.sigma_min = std::sqrt(lambda.front());
    out.sigma_max = std::sqrt(lambda.back());
    out.condition = out.sigma_min > 0.0 ? out.sigma_max / out.sigma_min
                                        : std::numeric_limits<double>::infinity();
    out.radius_min = rmin;
    out.radius_max = rmax;
    for (const std::size_t nb : touched) seen[nb] = 0;
    return out;
}

void check_stencil_diagnostics()
{
    for (const std::size_t n : {4u, 6u, 8u}) {
        const Grid grid = make_tet_grid(n);
        std::size_t min_size = std::numeric_limits<std::size_t>::max(), max_size = 0;
        int min_rank = 9, max_rank = 0;
        double min_sigma = std::numeric_limits<double>::max(), max_sigma = 0.0;
        double max_condition = 0.0, min_radius = std::numeric_limits<double>::max(), max_radius = 0.0;
        for (const std::size_t c : grid.interior) {
            const auto d = quadratic_stencil_diagnostics(grid, c);
            min_size = std::min(min_size, d.size); max_size = std::max(max_size, d.size);
            min_rank = std::min(min_rank, d.rank); max_rank = std::max(max_rank, d.rank);
            min_sigma = std::min(min_sigma, d.sigma_min); max_sigma = std::max(max_sigma, d.sigma_max);
            max_condition = std::max(max_condition, d.condition);
            min_radius = std::min(min_radius, d.radius_min); max_radius = std::max(max_radius, d.radius_max);
        }
        std::cout << "POLY_GRAD_DIAG n=" << n
                  << " stencil_size=[" << min_size << "," << max_size << "]"
                  << " rank=[" << min_rank << "," << max_rank << "]"
                  << " sigma_min=" << min_sigma << " sigma_max=" << max_sigma
                  << " cond_max=" << max_condition
                  << " radius=[" << min_radius << "," << max_radius << "]\n";
        require(min_rank == 9 && max_rank == 9,
                "least_squares_quadratic: all interior stencils must have full rank 9");
        require(std::isfinite(max_condition),
                "least_squares_quadratic: stencil condition numbers must be finite");
    }
}

void check_quadratic_reconstruction()
{
    for (const std::size_t n : {4u, 6u, 8u}) {
        const Grid grid = make_tet_grid(n);
        const auto phi = sample_field(grid, kQuadratic);
        const auto g = compute_scheme(grid, phi, GradScheme::LeastSquaresQuad);
        const auto e = gradient_error(grid, g, kQuadratic, true);
        std::cout << "POLY_GRAD scheme=least_squares_quadratic field=quadratic n=" << n
                  << " L2=" << e.l2 << " Linf=" << e.linf << "\n";
        require(e.linf <= 1e-8,
                "least_squares_quadratic: quadratic field must be exact on interior tetrahedra");
    }
}

void check_cubic_basis()
{
    const std::array<FieldCase,10> fields{{
        {"x3", cubic_x3, cubic_x3_gradient},
        {"y3", cubic_y3, cubic_y3_gradient},
        {"z3", cubic_z3, cubic_z3_gradient},
        {"x2y", cubic_x2y, cubic_x2y_gradient},
        {"x2z", cubic_x2z, cubic_x2z_gradient},
        {"y2x", cubic_y2x, cubic_y2x_gradient},
        {"y2z", cubic_y2z, cubic_y2z_gradient},
        {"z2x", cubic_z2x, cubic_z2x_gradient},
        {"z2y", cubic_z2y, cubic_z2y_gradient},
        {"xyz", cubic_xyz, cubic_xyz_gradient}
    }};
    for (const auto& field : fields) {
        std::vector<double> errors;
        for (const std::size_t n : {4u, 6u, 8u}) {
            const Grid grid = make_tet_grid(n);
            const auto phi = sample_field(grid, field);
            const auto e = gradient_error(grid,
                compute_scheme(grid, phi, GradScheme::LeastSquaresQuad), field, true);
            errors.push_back(e.l2);
            std::cout << "POLY_GRAD_CUBIC field=" << field.name << " n=" << n
                      << " L2=" << e.l2 << " Linf=" << e.linf;
            if (errors.size() > 1)
                std::cout << " order=" << observed_order(
                    errors[errors.size()-2], errors.back(),
                    static_cast<double>(n) / static_cast<double>(
                        std::array<std::size_t,3>{4u,6u,8u}[errors.size()-2]));
            std::cout << "\n";
        }
        const double o1 = observed_order(errors[0], errors[1], 6.0/4.0);
        const double o2 = observed_order(errors[1], errors[2], 8.0/6.0);
        require(std::isfinite(errors[2]) && o1 >= 1.8 && o2 >= 1.8,
                std::string("least_squares_quadratic: cubic field ") + field.name +
                " must show >= 1.8 gradient order");
    }
}

void check_order()
{
    for (const GradScheme scheme : {GradScheme::GaussCell, GradScheme::GaussVertex,
                                    GradScheme::GaussPoint, GradScheme::LeastSquares, GradScheme::LeastSquares2,
                                    GradScheme::LeastSquaresQuad}) {
        std::vector<double> errors;
        for (const std::size_t n : {4u, 6u, 8u}) {
            const Grid grid = make_tet_grid(n);
            const auto phi = sample_field(grid, kSmooth);
            const auto e = gradient_error(grid, compute_scheme(grid, phi, scheme), kSmooth, true);
            errors.push_back(e.l2);
            std::cout << "POLY_GRAD_ORDER scheme=" << scheme_name(scheme) << " n=" << n
                      << " L2=" << e.l2 << " Linf=" << e.linf;
            if (errors.size() > 1)
                std::cout << " order=" << observed_order(errors[errors.size() - 2], errors.back());
            std::cout << "\n";
        }
        require(std::isfinite(errors.back()),
                scheme_name(scheme) + ": tet gradient errors must be finite");
        if (scheme == GradScheme::LeastSquaresQuad) {
            // The tetrahedral campaign refines n=4 -> 6 -> 8, so h is
            // proportional to 1/n and the refinement ratios are 6/4 and 8/6.
            // Do not use the verification helper default (ratio=2): that would
            // report a false sub-second-order result for this non-dyadic mesh sequence.
            const double o1 = observed_order(errors[0], errors[1], 6.0 / 4.0);
            const double o2 = observed_order(errors[1], errors[2], 8.0 / 6.0);
            require(o1 >= 1.8 && o2 >= 1.8,
                    "least_squares_quadratic: tetrahedral smooth-field order must be >= 1.8");
        } else if (scheme == GradScheme::GaussCell) {
            // Two-point Gauss is the documented skew-inconsistent scheme.
        } else {
            require(errors.back() < errors.front(),
                    scheme_name(scheme) + ": tet refinement must reduce the error");
        }
        // The measured smooth-field orders on this tetrahedral stencil
        // (cell Green-Gauss ~0.1-0.2 and non-consistent; vertex ~0.4-0.65;
        // least-squares ~0.4; skew-corrected point-linear ~0.5-0.9) are
        // documented as a finding: linear-consistency is achieved by LS/vertex/
        // point, but none reaches second order at these resolutions, which is
        // exactly why the polyhedral scheme-accuracy gap stays open. Only the
        // robust invariants above are asserted.
    }
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(12);
        check_exactness(3);
        check_quadratic_reconstruction();
        check_stencil_diagnostics();
        check_cubic_basis();
        check_order();
        std::cout << "POLYHEDRAL_GRADIENT_CAMPAIGN: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "POLYHEDRAL_GRADIENT_CAMPAIGN: FAIL: " << e.what() << "\n";
        return 1;
    }
}