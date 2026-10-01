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
#include <limits>
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

Grid make_tet_grid(std::size_t n, double sx = 1.0, double sy = 1.0, double sz = 1.0,
                    double shear_xy = 0.0, double shear_xz = 0.0)
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
            for (std::size_t x = 0; x <= n; ++x) {
                const double qx = static_cast<double>(x) / static_cast<double>(n);
                const double qy = static_cast<double>(y) / static_cast<double>(n);
                const double qz = static_cast<double>(z) / static_cast<double>(n);
                m.points().set(vid(x, y, z),
                    sx * qx + shear_xy * qy + shear_xz * qz,
                    sy * qy,
                    sz * qz);
            }

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
inline const FieldCase kSmooth{"smooth", smooth_value, smooth_gradient};

enum class GradScheme { GaussCell, GaussVertex, GaussPoint, LeastSquares, LeastSquares2, LeastSquaresWeighted, LeastSquaresWeightedInvR, LeastSquaresQuad };
std::string scheme_name(GradScheme s)
{
    switch (s) {
        case GradScheme::GaussCell:   return "green_gauss";
        case GradScheme::GaussVertex: return "green_gauss_vertex";
        case GradScheme::GaussPoint:  return "green_gauss_point";
        case GradScheme::LeastSquares: return "least_squares";
        case GradScheme::LeastSquares2: return "least_squares_2ring";
        case GradScheme::LeastSquaresWeighted: return "weighted_least_squares_1_over_r2";
        case GradScheme::LeastSquaresWeightedInvR: return "weighted_least_squares_1_over_r";
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
        case GradScheme::LeastSquaresWeighted:
            return compute_gradient_weighted_least_squares_extended(phi, grid.mesh, GradientWeighting::INVERSE_DISTANCE_SQUARED);
        case GradScheme::LeastSquaresWeightedInvR:
            return compute_gradient_weighted_least_squares_extended(phi, grid.mesh, GradientWeighting::INVERSE_DISTANCE);
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
                                    GradScheme::LeastSquaresWeighted, GradScheme::LeastSquaresWeightedInvR,
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




double symmetric_condition_reference(const Vec3& centre,
                                     const std::vector<Vec3>& neighbours,
                                     GradientWeighting weighting)
{
    double a[3][3] = {};
    for (const Vec3& p : neighbours) {
        const Vec3 d = p - centre;
        const double r2 = d.mag2();
        if (!(r2 > 0.0) || !std::isfinite(r2))
            return std::numeric_limits<double>::infinity();
        const double r = std::sqrt(r2);
        double w = 1.0;
        if (weighting == GradientWeighting::INVERSE_DISTANCE) w = 1.0 / r;
        if (weighting == GradientWeighting::INVERSE_DISTANCE_SQUARED) w = 1.0 / r2;
        a[0][0] += w*d.x*d.x; a[0][1] += w*d.x*d.y; a[0][2] += w*d.x*d.z;
        a[1][0] += w*d.y*d.x; a[1][1] += w*d.y*d.y; a[1][2] += w*d.y*d.z;
        a[2][0] += w*d.z*d.x; a[2][1] += w*d.z*d.y; a[2][2] += w*d.z*d.z;
    }

    // Jacobi diagonalisation of the symmetric normal matrix.  This is an
    // independent reference for the spectral 2-norm condition number.
    for (int sweep = 0; sweep < 50; ++sweep) {
        int p = 0, q = 1;
        double max_off = std::abs(a[0][1]);
        if (std::abs(a[0][2]) > max_off) { max_off = std::abs(a[0][2]); p = 0; q = 2; }
        if (std::abs(a[1][2]) > max_off) { max_off = std::abs(a[1][2]); p = 1; q = 2; }
        if (!(max_off > 0.0) || !std::isfinite(max_off)) break;
        const double tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
        const double t = (tau >= 0.0 ? 1.0 : -1.0) /
                         (std::abs(tau) + std::sqrt(1.0 + tau*tau));
        const double cs = 1.0 / std::sqrt(1.0 + t*t);
        const double sn = t * cs;
        const double app = a[p][p], aqq = a[q][q], apq = a[p][q];
        a[p][p] = app - t*apq;
        a[q][q] = aqq + t*apq;
        a[p][q] = a[q][p] = 0.0;
        for (int k = 0; k < 3; ++k) {
            if (k == p || k == q) continue;
            const double akp = a[k][p], akq = a[k][q];
            a[k][p] = a[p][k] = cs*akp - sn*akq;
            a[k][q] = a[q][k] = sn*akp + cs*akq;
        }
    }
    double lo = a[0][0], hi = a[0][0];
    for (int i = 1; i < 3; ++i) {
        lo = std::min(lo, a[i][i]);
        hi = std::max(hi, a[i][i]);
    }
    if (!(lo > 0.0) || !std::isfinite(lo) || !std::isfinite(hi))
        return std::numeric_limits<double>::infinity();
    return hi / lo;
}

void check_wls_conditioning_reference_campaign()
{
    struct Family { const char* name; std::vector<Vec3> points; };
    const double e = 1.0e-6;
    const std::vector<Family> families = {
        {"isotropic", {
            {1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}}},
        {"anisotropic", {
            {1,0,0},{-1,0,0},{0,0.1,0},{0,-0.1,0},{0,0,1.7},{0,0,-1.7}}},
        {"quasi_coplanar", {
            {1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0.3,0.7,e},{-0.3,-0.7,-e}}},
        {"quasi_collinear", {
            {1,0,0},{-1,0,0},{2,e,0},{-2,-e,0},{3,0,e},{-3,0,-e}}}
    };
    const std::array<double,3> scales = {1.0, 1.0e-4, 1.0e4};
    const std::array<GradientWeighting,3> weights = {
        GradientWeighting::UNIFORM,
        GradientWeighting::INVERSE_DISTANCE,
        GradientWeighting::INVERSE_DISTANCE_SQUARED
    };
    const char* weight_names[] = {"uniform", "1/r", "1/r2"};

    std::cout << "WLS_REFERENCE family weighting scale rank pivot_ratio spectral_condition rejected\n";
    for (const auto& family : families) {
        for (const auto weighting : weights) {
            for (const double scale : scales) {
                std::vector<Vec3> pts;
                pts.reserve(family.points.size());
                for (const Vec3& p : family.points)
                    pts.push_back(Vec3{scale*p.x, scale*p.y, scale*p.z});

                StencilQuality q;
                bool rejected = false;
                try {
                    (void)weighted_least_squares_gradient(
                        Vec3{0,0,0}, 0.0, pts,
                        std::vector<double>(pts.size(), 0.0),
                        weighting, &q);
                } catch (const std::runtime_error&) {
                    rejected = true;
                }

                const double reference = symmetric_condition_reference(
                    Vec3{0,0,0}, pts, weighting);
                std::cout << "WLS_REFERENCE family=" << family.name
                          << " weighting=" << weight_names[static_cast<int>(weighting)]
                          << " scale=" << scale
                          << " rank=" << q.rank
                          << " pivot_ratio=" << q.condition_estimate
                          << " spectral_condition=" << reference
                          << " rejected=" << (rejected ? "true" : "false") << "\n";

                require(std::isfinite(reference),
                        "conditioning reference: spectral condition must remain finite for characterized stencils");
                require(q.condition_estimate >= 1.0 || std::isinf(q.condition_estimate),
                        "conditioning estimator must be >= 1 or infinite");

                const bool expected_rank_rejection =
                    std::string(family.name) == "quasi_collinear";
                if (expected_rank_rejection) {
                    require(q.rank < 3,
                            "conditioning reference: quasi-collinear family must expose numerical rank loss");
                    require(rejected,
                            "conditioning reference: quasi-collinear family must be rejected by WLS");
                    require(std::isinf(q.condition_estimate),
                            "conditioning reference: numerical rank loss must report infinite pivot diagnostic");
                } else {
                    require(q.rank == 3,
                            "conditioning reference: non-collinear family must remain full rank");
                    require(!rejected,
                            "conditioning reference: full-rank characterized family must not be rejected");
                }
            }
        }
    }
}

void check_wls_conditioning()
{
    // Characterization campaign only. This exercises the monotonic conditioning
    // trend and the numerical-rank transition; it does not define a production
    // acceptance threshold.
    const std::array<double, 6> epsilons = {
        1.0, 1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10};

    double previous_pivot = 0.0;
    double previous_reference = 0.0;
    bool saw_rank_rejection = false;

    for (const double eps : epsilons) {
        const Vec3 centre{0.0, 0.0, 0.0};
        const std::vector<Vec3> neighbours = {
            {1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0},
            {0.0, eps, 0.0}, {0.0, -eps, 0.0},
            {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}};
        const std::vector<double> values(neighbours.size(), 0.0);

        StencilQuality q;
        bool rank_rejected = false;
        try {
            (void)weighted_least_squares_gradient(
                centre, 0.0, neighbours, values,
                GradientWeighting::UNIFORM, &q);
        } catch (const std::runtime_error&) {
            rank_rejected = true;
        }

        const double reference =
            symmetric_condition_reference(centre, neighbours, GradientWeighting::UNIFORM);
        require(std::isfinite(reference),
                "conditioning campaign: independent spectral reference must remain finite");

        if (rank_rejected) {
            std::cout << "WLS_CONDITIONING eps=" << eps
                      << " pivot=" << q.condition_estimate
                      << " spectral_condition=" << reference
                      << " class=REJECT_RANK rank=" << q.rank << "\n";
            require(q.rank < 3,
                    "conditioning campaign: rank rejection must report rank loss");
            require(std::isinf(q.condition_estimate),
                    "conditioning campaign: rank loss must report infinite pivot diagnostic");
            require(reference > previous_reference,
                    "conditioning campaign: independent reference must continue increasing at rank transition");
            saw_rank_rejection = true;
            break;
        }

        require(q.full_rank,
                "conditioning campaign: non-rejected stencil must be full rank");
        require(std::isfinite(q.condition_estimate),
                "conditioning campaign: full-rank stencil must have finite pivot diagnostic");
        if (previous_pivot > 0.0)
            require(q.condition_estimate > previous_pivot,
                    "conditioning campaign: pivot diagnostic must increase monotonically");
        if (previous_reference > 0.0)
            require(reference > previous_reference,
                    "conditioning campaign: independent spectral reference must increase monotonically");

        std::cout << "WLS_CONDITIONING eps=" << eps
                  << " pivot=" << q.condition_estimate
                  << " spectral_condition=" << reference
                  << " class=FULL_RANK rank=" << q.rank << "\n";
        previous_pivot = q.condition_estimate;
        previous_reference = reference;
    }

    require(saw_rank_rejection,
            "conditioning campaign: rank-based rejection was not exercised");
}

void check_distorted_order()
{
    // Fixed affine distortion: 10:1 aspect ratio in y plus x-y/x-z shear.
    // The physical field is evaluated at the transformed cell centres, so
    // this probes the gradient operator on a consistently refined,
    // high-aspect-ratio/skewed tetrahedral family.
    const std::array<GradScheme, 2> schemes = {
        GradScheme::LeastSquaresWeighted, GradScheme::LeastSquaresWeightedInvR};
    for (const GradScheme scheme : schemes) {
        std::vector<double> errors;
        for (const std::size_t n : {4u, 8u, 16u}) {
            const Grid grid = make_tet_grid(n, 1.0, 0.1, 1.8, 0.35, 0.20);
            const auto phi = sample_field(grid, kSmooth);
            const auto e = gradient_error(grid,
                compute_scheme(grid, phi, scheme), kSmooth, true);
            require(std::isfinite(e.l2), scheme_name(scheme) + ": distorted L2 must be finite");
            errors.push_back(e.l2);
            std::cout << "POLY_GRAD_DISTORTED scheme=" << scheme_name(scheme)
                      << " n=" << n << " L2=" << e.l2 << " Linf=" << e.linf;
            if (errors.size() > 1)
                std::cout << " order=" << observed_order(errors[errors.size()-2], errors.back());
            std::cout << "\n";
        }
        require(errors.back() < errors.front(),
                scheme_name(scheme) + ": distorted tetra refinement must reduce error");
    }
}

void check_boundary_reconstruction()
{
    const Grid grid = make_tet_grid(4, 1.0, 0.8, 1.4, 0.25, 0.15);
    const auto phi = sample_field(grid, kLinear);
    const auto& own = grid.mesh.ownership();

    std::vector<BoundaryGradientCondition> dirichlet(grid.mesh.n_faces());
    std::vector<BoundaryGradientCondition> neumann(grid.mesh.n_faces());
    std::vector<BoundaryGradientCondition> mixed(grid.mesh.n_faces());

    std::vector<std::size_t> boundary_cells;
    std::vector<std::size_t> corner_cells;
    std::vector<unsigned> boundary_count(grid.mesh.n_cells(), 0);

    for (std::size_t f = 0; f < grid.mesh.n_faces(); ++f) {
        if (own.neighbour(f) >= 0) continue;
        const std::size_t c = own.owner(f);
        ++boundary_count[c];

        const Vec3 n = grid.geometry.face_normals[f];
        const double face_value = linear_value(grid.geometry.face_centres[f]);
        const double normal_gradient =
            linear_gradient(grid.geometry.face_centres[f]).dot(n);

        dirichlet[f] = {
            BoundaryGradientConditionType::DIRICHLET, face_value};
        neumann[f] = {
            BoundaryGradientConditionType::NEUMANN, normal_gradient};

        // Alternate BC type by boundary face so cells at edges/corners
        // necessarily exercise mixed reconstruction.
        mixed[f] = (f % 2 == 0)
            ? BoundaryGradientCondition{
                  BoundaryGradientConditionType::DIRICHLET, face_value}
            : BoundaryGradientCondition{
                  BoundaryGradientConditionType::NEUMANN, normal_gradient};
    }

    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c) {
        if (boundary_count[c] > 0) boundary_cells.push_back(c);
        if (boundary_count[c] >= 2) corner_cells.push_back(c);
    }
    require(!boundary_cells.empty(),
            "tet boundary campaign requires boundary-adjacent cells");
    require(!corner_cells.empty(),
            "tet boundary campaign requires corner/edge cells");

    const auto boundary_error = [&](const Field<double, Location::CELL>& grad,
                                     const std::vector<std::size_t>& cells) {
        std::vector<double> got, want, weights;
        for (const std::size_t c : cells) {
            const Vec3 exact = kLinear.exact_gradient(grid.geometry.cell_centres[c]);
            const Vec3 numerical{grad(c, 0), grad(c, 1), grad(c, 2)};
            const double err = (numerical - exact).mag();
            got.push_back(exact.mag() + err);
            want.push_back(exact.mag());
            weights.push_back(grid.geometry.cell_volumes[c]);
        }
        return error_norms(got, want, weights);
    };

    const std::array<std::pair<const char*, const std::vector<BoundaryGradientCondition>*>, 3> cases = {{
        {"dirichlet", &dirichlet},
        {"neumann", &neumann},
        {"mixed", &mixed}
    }};

    for (const auto& [name, conditions] : cases) {
        const auto grad = compute_gradient_weighted_least_squares(
            phi, grid.mesh, GradientWeighting::INVERSE_DISTANCE_SQUARED,
            std::numeric_limits<double>::infinity(),
            BoundaryGradientPolicy::ZERO_GRADIENT_GHOST, conditions);

        const auto all = boundary_error(grad, boundary_cells);
        const auto corners = boundary_error(grad, corner_cells);
        std::cout << "POLY_GRAD_BOUNDARY case=" << name
                  << " cells=" << boundary_cells.size()
                  << " L1=" << all.l1 << " L2=" << all.l2
                  << " Linf=" << all.linf << "\n";
        std::cout << "POLY_GRAD_BOUNDARY_CORNER case=" << name
                  << " cells=" << corner_cells.size()
                  << " L1=" << corners.l1 << " L2=" << corners.l2
                  << " Linf=" << corners.linf << "\n";

        require(std::isfinite(all.l1) && std::isfinite(all.l2) &&
                std::isfinite(all.linf),
                std::string("tet boundary ") + name + ": errors must be finite");
        require(all.linf <= 1e-9,
                std::string("tet boundary ") + name +
                ": linear field must remain exact");
        require(corners.linf <= 1e-9,
                std::string("tet boundary ") + name +
                ": corner/edge cells must remain exact");
    }

    const auto constant = sample_field(grid, kConstant);
    const auto zg = compute_gradient_weighted_least_squares(
        constant, grid.mesh, GradientWeighting::INVERSE_DISTANCE_SQUARED,
        std::numeric_limits<double>::infinity(),
        BoundaryGradientPolicy::ZERO_GRADIENT_GHOST);
    // For the constant-field invariant the exact gradient is zero; do not
    // reuse boundary_error(), which is intentionally defined against the
    // linear-field reference gradient.
    double zero_linf = 0.0;
    for (const std::size_t c : boundary_cells) {
        const Vec3 numerical{zg(c, 0), zg(c, 1), zg(c, 2)};
        zero_linf = std::max(zero_linf, numerical.mag());
    }
    std::cout << "POLY_GRAD_BOUNDARY case=zero_gradient_default"
              << " Linf=" << zero_linf << "\n";
    require(zero_linf <= 1e-12,
            "tet zero-gradient default ghost must preserve constants");

    bool rejected = false;
    try {
        (void)compute_gradient_weighted_least_squares(
            phi, grid.mesh, GradientWeighting::INVERSE_DISTANCE_SQUARED,
            std::numeric_limits<double>::infinity(),
            BoundaryGradientPolicy::REJECT_BOUNDARY_STENCIL);
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    require(rejected,
            "tet REJECT_BOUNDARY_STENCIL must reject boundary stencils");
}

void check_order()
{
    for (const GradScheme scheme : {GradScheme::GaussCell, GradScheme::GaussVertex,
                                    GradScheme::GaussPoint, GradScheme::LeastSquares, GradScheme::LeastSquares2,
                                    GradScheme::LeastSquaresWeighted, GradScheme::LeastSquaresWeightedInvR,
                                    GradScheme::LeastSquaresQuad}) {
        std::vector<double> errors;
        for (const std::size_t n : {4u, 8u, 16u}) {
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
        if (scheme == GradScheme::GaussCell) {
            // Two-point Gauss is not even linear-consistent on tetrahedra.
        } else if (scheme == GradScheme::LeastSquaresQuad) {
            require(errors.back() < errors.front(),
                    "quadratic LS: tet refinement must reduce the error");
            require_order(errors, 2.0, 1.5,
                          "quadratic least-squares gradient on tetrahedra");
        } else {
            // Vertex GG, point-linear GG and the least-squares variants are
            // linear-consistent; refinement must reduce the error.
            require(errors.back() < errors.front(),
                    scheme_name(scheme) + ": tet refinement must reduce the error");
        }
// Smooth-field observed orders on the tetrahedral grid with a 2:1
    // refinement sweep (n = 4, 8, 16):
    //   two-point GG      ~0.25/0.08  (non-consistent)
    //   vertex GG         ~1.06/0.98  (converging, first order)
    //   point-linear GG   ~1.45/1.08
    //   least-squares     ~0.77/0.96
    //   2-ring LS         ~1.09/1.01
    //   quadratic LS      ~2.13/2.01  (**second order**)
    // So only the quadratic-basis least-squares reaches second order on
    // tetrahedra; the others are first-order-converging (except the two-point
    // Gauss scheme, which is not even linear-consistent there).
    }
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(12);
        check_exactness(3);
        check_order();
        check_distorted_order();
        check_boundary_reconstruction();
        check_wls_conditioning_reference_campaign();
        check_wls_conditioning();
        std::cout << "POLYHEDRAL_GRADIENT_CAMPAIGN: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "POLYHEDRAL_GRADIENT_CAMPAIGN: FAIL: " << e.what() << "\n";
        return 1;
    }
}