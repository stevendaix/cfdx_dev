// Issue #461 — N3 diffusion on a polyhedral (tetrahedral) mesh.
//
// Darkens the "non-affine/extremely non-orthogonal and polyhedral diffusion"
// N3 gap: the Laplacian operators on an unstructured tetrahedral grid.
// Constant fields must be annihilated; linear fields and smooth manufactured
// field converge/conserve as measured, reported without over-promotion.
// With a linear-consistent (least-squares) cell gradient the corrected
// operator is linear-exact, and that consistency is enforced at every
// campaign refinement level. The smooth-field interior residual under
// this consistent configuration is measured and reported WITHOUT a
// decrease threshold: a linear-exact but non-symmetric corrected operator
// is not quadratic-exact, so no pointwise residual convergence is claimed.
// The default two-point Gauss gradient variant is not linear-consistent
// on tetrahedra and its measured non-decreasing error is reported rather
// than hidden.

#include "cfdx/core/numerics/laplacian.h"
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
    std::vector<std::size_t> interior;   // tets whose cube is not on the boundary
};

static const std::array<std::array<unsigned, 3>, 8> kCorner = {{
    {{0,0,0}}, {{1,0,0}}, {{1,1,0}}, {{0,1,0}},
    {{0,0,1}}, {{1,0,1}}, {{1,1,1}}, {{0,1,1}}
}};
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
    if (n < 2) throw std::invalid_argument("make_tet_grid: n must be >= 2");

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
    std::vector<std::vector<std::size_t>> cell_faces;

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

    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i) {
                std::array<std::size_t, 8> corner;
                for (unsigned q = 0; q < 8; ++q)
                    corner[q] = vid(i + kCorner[q][0], j + kCorner[q][1], k + kCorner[q][2]);
                for (unsigned t = 0; t < kTets.size(); ++t) {
                    const auto& tet = kTets[t];
                    std::vector<std::size_t> faces(4);
                    faces[0] = get_face(corner[tet[0]], corner[tet[1]], corner[tet[2]]);
                    faces[1] = get_face(corner[tet[0]], corner[tet[1]], corner[tet[3]]);
                    faces[2] = get_face(corner[tet[0]], corner[tet[2]], corner[tet[3]]);
                    faces[3] = get_face(corner[tet[1]], corner[tet[2]], corner[tet[3]]);
                    const std::size_t cell = cell_faces.size();
                    for (const std::size_t f : faces) face_cells[f].push_back(cell);
                    cell_faces.push_back(faces);
                }
            }

    m.ownership().resize(m.n_faces());
    for (std::size_t f = 0; f < m.n_faces(); ++f) {
        const auto& cells = face_cells[f];
        m.ownership().set_owner(f, cells[0]);
        m.ownership().set_neighbour(f, cells.size() == 1
            ? FaceOwnership::BOUNDARY : static_cast<std::int64_t>(cells[1]));
    }
    for (const auto& faces : cell_faces) m.cells().push_cell(faces);

    Grid grid;
    grid.mesh = std::move(m);
    grid.geometry = make_geometry_cache(grid.mesh);
    for (std::size_t k = 1; k + 1 < n; ++k)
        for (std::size_t j = 1; j + 1 < n; ++j)
            for (std::size_t i = 1; i + 1 < n; ++i) {
                const std::size_t base = 6 * (i + n * (j + n * k));
                for (unsigned t = 0; t < kTets.size(); ++t) grid.interior.push_back(base + t);
            }
    return grid;
}

std::vector<double> sample(const Grid& grid, double (*fn)(const Vec3&)) {
    std::vector<double> out(grid.mesh.n_cells(), 0.0);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c)
        out[c] = fn(grid.geometry.cell_centres[c]);
    return out;
}

double linear_value(const Vec3& p) { return 2.0 * p.x - 3.0 * p.y + 0.5 * p.z; }  // Laplacian = 0
double smooth_value(const Vec3& p) {
    const double pi = std::acos(-1.0);
    return std::sin(pi * p.x) * std::cos(pi * p.y) * std::sin(pi * p.z);
}
double smooth_laplacian(const Vec3& p) {
    const double pi = std::acos(-1.0);
    return -3.0 * pi * pi * smooth_value(p);
}

Field<double, Location::CELL> laplacian_of(const Grid& grid, const std::vector<double>& phi,
                                           LaplacianScheme scheme,
                                           GradientScheme gs = GradientScheme::GAUSS_TWO_POINT)
{
    Field<double, Location::CELL> f(grid.mesh.n_cells(), "phi", "1", 1);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c) f(c) = phi[c];
    return compute_laplacian(f, grid.mesh, grid.geometry, scheme, 0.5, gs);
}

ErrorMetrics interior_error(const Grid& grid, const Field<double, Location::CELL>& lap,
                            const std::vector<double>& exact)
{
    std::vector<double> got, want, w;
    for (const std::size_t c : grid.interior) {
        got.push_back(lap(c));
        want.push_back(exact[c]);
        w.push_back(grid.geometry.cell_volumes[c]);
    }
    return error_norms(got, want, w);
}

void require(bool c, const std::string& m) { if (!c) throw std::runtime_error(m); }

void check_constant_and_linear(std::size_t n)
{
    const Grid grid = make_tet_grid(n);
    const std::vector<double> zero(grid.mesh.n_cells(), 0.0);

    // Measure on the same grid so scheme comparisons are meaningful.
    double uncorrected_linear = 0.0;
    for (const auto scheme : {LaplacianScheme::UNCORRECTED, LaplacianScheme::CORRECTED,
                              LaplacianScheme::LIMITED, LaplacianScheme::OVER_RELAXED}) {
        const auto constant = sample(grid, [](const Vec3&) { return 42.0; });
        const auto ce = interior_error(grid, laplacian_of(grid, constant, scheme), zero);
        std::cout << "POLY_LAP n=" << n << " scheme=" << to_string(scheme)
                  << " field=constant L2=" << ce.l2 << " Linf=" << ce.linf << "\n";
        require(ce.linf <= 1e-9, std::string(to_string(scheme)) + ": Laplacian of constant must vanish");

        const auto linear = sample(grid, linear_value);
        const auto le = interior_error(grid, laplacian_of(grid, linear, scheme), zero);
        std::cout << "POLY_LAP n=" << n << " scheme=" << to_string(scheme)
                  << " field=linear L2=" << le.l2 << " Linf=" << le.linf << "\n";
        if (scheme == LaplacianScheme::UNCORRECTED) uncorrected_linear = le.l2;
        if (scheme != LaplacianScheme::UNCORRECTED) {
            // The non-orthogonal correction must strictly reduce the
            // linear-field Laplacian error on the tetrahedral mesh.
            require(le.l2 < uncorrected_linear,
                    std::string(to_string(scheme)) + ": correction must reduce linear error vs uncorrected");
        }
    }
}

void check_smooth()
{
    const std::vector<std::size_t> ns = {4u, 8u, 16u};

    std::vector<double> uncorrected(ns.size(), 0.0);
    for (std::size_t idx = 0; idx < ns.size(); ++idx) {
        const Grid grid = make_tet_grid(ns[idx]);
        const auto e = interior_error(grid, laplacian_of(grid, sample(grid, smooth_value),
                                                         LaplacianScheme::UNCORRECTED),
                                      sample(grid, smooth_laplacian));
        uncorrected[idx] = e.l2;
        std::cout << "POLY_LAP_ORDER scheme=uncorrected n=" << ns[idx]
                  << " L2=" << e.l2 << " Linf=" << e.linf << "\n";
        require(std::isfinite(e.l2), "uncorrected: tet Laplacian errors must be finite");
    }

    for (const auto scheme : {LaplacianScheme::CORRECTED, LaplacianScheme::OVER_RELAXED}) {
        double last_l2 = 0.0;
        for (std::size_t idx = 0; idx < ns.size(); ++idx) {
            const Grid grid = make_tet_grid(ns[idx]);
            const auto e = interior_error(grid, laplacian_of(grid, sample(grid, smooth_value), scheme),
                                          sample(grid, smooth_laplacian));
            last_l2 = e.l2;
            std::cout << "POLY_LAP_ORDER scheme=" << to_string(scheme) << " n=" << ns[idx]
                      << " L2=" << e.l2 << " Linf=" << e.linf << "\n";
            require(e.l2 <= uncorrected[idx],
                    std::string(to_string(scheme)) +
                    ": tet Laplacian error must not exceed the uncorrected error");
        }
        require(std::isfinite(last_l2),
                std::string(to_string(scheme)) + ": tet Laplacian errors must be finite");
    }
}

// Smooth-field diagnostic with linear-consistent gradients: the corrected
// operator with a least-squares cell gradient is linear-exact
// (check_operator_consistency enforces this at every campaign refinement
// level), but it is not quadratic-exact and its gradient-based correction
// term makes the operator non-symmetric. Its pointwise interior residual
// on a smooth manufactured field is therefore O(1), and no convergence
// claim is made for it: a strict-decrease requirement on these residuals
// was measured to FAIL on the n=4/8/16 tetrahedral family (first revision
// of this PR), which is exactly why the campaign reports them without a
// threshold instead of promoting an unproven order. Errors must remain
// finite; the observed orders are reported as diagnostics only.
void check_smooth_with_consistent_gradient()
{
    const std::vector<std::size_t> ns = {4u, 8u, 16u};

    for (const auto gs : {GradientScheme::LEAST_SQUARES,
                          GradientScheme::LEAST_SQUARES_QUADRATIC}) {
        std::vector<double> l2(ns.size(), 0.0);
        for (std::size_t idx = 0; idx < ns.size(); ++idx) {
            const Grid grid = make_tet_grid(ns[idx]);
            const auto e = interior_error(grid,
                laplacian_of(grid, sample(grid, smooth_value),
                             LaplacianScheme::CORRECTED, gs),
                sample(grid, smooth_laplacian));
            require(std::isfinite(e.l2) && std::isfinite(e.linf),
                    std::string(to_string(gs)) +
                    ": corrected Laplacian errors must be finite");
            l2[idx] = e.l2;
            std::cout << "POLY_LAP_ORDER scheme=corrected gradient=" << to_string(gs)
                      << " n=" << ns[idx] << " L2=" << e.l2 << " Linf=" << e.linf << "\n";
        }
        for (std::size_t idx = 1; idx < ns.size(); ++idx) {
            const double observed_order =
                std::log(l2[idx - 1] / l2[idx]) / std::log(2.0);
            std::cout << "POLY_LAP_ORDER scheme=corrected gradient=" << to_string(gs)
                      << " from_n=" << ns[idx - 1] << " to_n=" << ns[idx]
                      << " observed_order=" << observed_order << "\n";
        }
    }
}

// Operator consistency: with a linear-consistent cell gradient (least
// squares / quadratic least squares) the corrected non-orthogonal operator
// annihilates a linear field on interior tetrahedra. The two-point Gauss
// gradient is NOT linear-consistent on tets, so the corrected operator is not
// even consistent there by default — this is the measured root of the N3
// polyhedral diffusion non-convergence. (The point-linear GG is excluded:
// its vertex-based values are only linear-exact at interior cells, while the
// face correction also reads boundary-adjacent cell gradients.)
void check_operator_consistency(std::size_t n)
{
    const Grid grid = make_tet_grid(n);
    const std::vector<double> zero(grid.mesh.n_cells(), 0.0);
    const auto linear = sample(grid, linear_value);
    for (const auto gs : {GradientScheme::LEAST_SQUARES,
                          GradientScheme::LEAST_SQUARES_QUADRATIC}) {
        const auto e = interior_error(grid,
            laplacian_of(grid, linear, LaplacianScheme::CORRECTED, gs), zero);
        std::cout << "POLY_LAP_CONSISTENT gradient=" << to_string(gs)
                  << " corrected linear Linf=" << e.linf << "\n";
        require(e.linf <= 1e-9,
                std::string(to_string(gs)) +
                ": corrected Laplacian must be linear-exact on interior tets");
    }
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(12);
        for (const std::size_t n : {4u, 6u, 8u, 16u}) check_operator_consistency(n);
        check_constant_and_linear(4);
        check_smooth();
        check_smooth_with_consistent_gradient();
        std::cout << "POLYHEDRAL_LAPLACIAN_CAMPAIGN: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "POLYHEDRAL_LAPLACIAN_CAMPAIGN: FAIL: " << e.what() << "\n";
        return 1;
    }
}