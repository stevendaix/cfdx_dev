// Issue #461 — N10 difficult-mesh numerical robustness campaign.
//
// N10 is deliberately a mesh-robustness package, not a new discretisation
// family.  This campaign provides one executable quality ladder covering:
//   orthogonal -> mildly skewed -> strongly skewed -> high-aspect-ratio
//   -> polyhedral -> near-degenerate-but-valid -> invalid/rejected.
//
// For every valid structured case the campaign records mesh quality and solves
// the corresponding two-point finite-volume diffusion matrix.  The linear
// system is assembled from the actual CFDX face geometry/topology; it is not a
// synthetic matrix.  A manufactured algebraic solution is used so the solver
// gate is independent of an external PDE oracle.
//
// The polyhedral member is used for geometry/topology robustness and is not
// promoted to a second-order accuracy claim here.  N10 must not silently
// convert mesh degradation into a different numerical method.

#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/geometry/mesh_validator.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/mesh/mesh.h"
#include "common/test_harness.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

struct Grid {
    Mesh mesh;
    GeometryCache geometry;
};

struct QualityCase {
    const char* name;
    double shear;
    double stretch;
    std::size_t n;
};

std::string face_key(std::vector<std::size_t> vertices)
{
    std::sort(vertices.begin(), vertices.end());
    std::ostringstream out;
    for (const auto v : vertices) out << v << '/';
    return out.str();
}

// validate_mesh() requires every boundary face to belong to a patch, so each
// generated mesh has to register its boundary explicitly.
void add_boundary_patches(Mesh& mesh, const std::string& name)
{
    Patch patch{name, PatchType::WALL, {}};
    for (std::size_t f = 0; f < mesh.n_faces(); ++f)
        if (mesh.ownership().neighbour(f) == FaceOwnership::BOUNDARY)
            patch.face_ids.push_back(f);
    if (patch.face_ids.empty())
        throw std::runtime_error("generated mesh has no boundary faces");
    mesh.boundary().add_patch(patch);
}

Grid make_hex_grid(std::size_t n, double shear, double stretch)
{
    if (n < 2) throw std::invalid_argument("make_hex_grid: n must be >= 2");
    const std::size_t nv = n + 1;
    const auto vid = [nv](std::size_t i, std::size_t j, std::size_t k) {
        return i + nv * (j + nv * k);
    };
    const auto transform = [shear, stretch](double x, double y, double z) {
        return std::array<double, 3>{
            x + shear * y,
            stretch * y,
            z
        };
    };

    Mesh mesh;
    mesh.points().resize(nv * nv * nv);
    for (std::size_t k = 0; k <= n; ++k)
        for (std::size_t j = 0; j <= n; ++j)
            for (std::size_t i = 0; i <= n; ++i) {
                const double x = static_cast<double>(i) / static_cast<double>(n);
                const double y = static_cast<double>(j) / static_cast<double>(n);
                const double z = static_cast<double>(k) / static_cast<double>(n);
                const auto p = transform(x, y, z);
                mesh.points().set(vid(i, j, k), p[0], p[1], p[2]);
            }

    std::unordered_map<std::string, std::size_t> face_ids;
    std::vector<std::vector<std::size_t>> cell_faces;
    const auto cell_id = [n](std::size_t i, std::size_t j, std::size_t k) {
        return (k * n + j) * n + i;
    };

    auto get_face = [&](std::vector<std::size_t> vertices) {
        const auto key = face_key(vertices);
        const auto it = face_ids.find(key);
        if (it != face_ids.end()) return it->second;
        const std::size_t id = mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        face_ids.emplace(key, id);
        return id;
    };

    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t c = cell_id(i, j, k);
                const std::size_t p000 = vid(i, j, k);
                const std::size_t p100 = vid(i + 1, j, k);
                const std::size_t p110 = vid(i + 1, j + 1, k);
                const std::size_t p010 = vid(i, j + 1, k);
                const std::size_t p001 = vid(i, j, k + 1);
                const std::size_t p101 = vid(i + 1, j, k + 1);
                const std::size_t p111 = vid(i + 1, j + 1, k + 1);
                const std::size_t p011 = vid(i, j + 1, k + 1);

                const std::array<std::vector<std::size_t>, 6> faces = {{
                    {p000, p001, p011, p010},
                    {p100, p110, p111, p101},
                    {p000, p100, p101, p001},
                    {p010, p011, p111, p110},
                    {p000, p010, p110, p100},
                    {p001, p101, p111, p011}
                }};
                std::vector<std::size_t> ids;
                ids.reserve(6);
                for (const auto& f : faces) ids.push_back(get_face(f));
                if (cell_faces.size() != c) throw std::runtime_error("cell numbering mismatch");
                cell_faces.push_back(std::move(ids));
            }

    mesh.ownership().resize(mesh.n_faces());
    for (const auto& entry : face_ids) {
        const std::size_t f = entry.second;
        std::vector<std::size_t> cells;
        for (std::size_t c = 0; c < cell_faces.size(); ++c)
            if (std::find(cell_faces[c].begin(), cell_faces[c].end(), f) != cell_faces[c].end())
                cells.push_back(c);

        if (cells.empty() || cells.size() > 2)
            throw std::runtime_error("invalid face adjacency in generated mesh");
        mesh.ownership().set_owner(f, cells[0]);
        mesh.ownership().set_neighbour(
            f, cells.size() == 1
                ? FaceOwnership::BOUNDARY
                : static_cast<std::int64_t>(cells[1]));
    }

    for (const auto& faces : cell_faces) mesh.cells().push_cell(faces);

    add_boundary_patches(mesh, "boundary");

    Grid grid{std::move(mesh), {}};
    grid.geometry = make_geometry_cache(grid.mesh);
    return grid;
}

struct TetraGrid {
    Mesh mesh;
    GeometryCache geometry;
};

TetraGrid make_tetra_grid(std::size_t n)
{
    if (n < 2) throw std::invalid_argument("make_tetra_grid: n must be >= 2");
    const std::size_t nv = n + 1;
    const auto vid = [nv](std::size_t i, std::size_t j, std::size_t k) {
        return i + nv * (j + nv * k);
    };
    const std::array<std::array<unsigned, 4>, 6> tets = {{
        {{0,1,2,6}}, {{0,1,5,6}}, {{0,3,2,6}},
        {{0,3,7,6}}, {{0,4,5,6}}, {{0,4,7,6}}
    }};

    Mesh mesh;
    mesh.points().resize(nv * nv * nv);
    for (std::size_t k = 0; k <= n; ++k)
        for (std::size_t j = 0; j <= n; ++j)
            for (std::size_t i = 0; i <= n; ++i)
                mesh.points().set(
                    vid(i,j,k),
                    static_cast<double>(i) / static_cast<double>(n),
                    static_cast<double>(j) / static_cast<double>(n),
                    static_cast<double>(k) / static_cast<double>(n));

    std::unordered_map<std::string, std::size_t> face_ids;
    std::vector<std::vector<std::size_t>> cell_faces;
    std::vector<std::vector<std::size_t>> face_cells;

    auto get_face = [&](std::size_t a, std::size_t b, std::size_t c,
                            const Vec3& owner_centre) {
        const auto key = face_key({a,b,c});
        const auto it = face_ids.find(key);
        if (it != face_ids.end()) return it->second;
        const Vec3 A{mesh.points().x(a), mesh.points().y(a), mesh.points().z(a)};
        const Vec3 B{mesh.points().x(b), mesh.points().y(b), mesh.points().z(b)};
        const Vec3 C{mesh.points().x(c), mesh.points().y(c), mesh.points().z(c)};
        // Half of the Kuhn subdivision tets are left-handed, so the vertex
        // order of a face is not by itself an outward orientation. Store the
        // winding that points away from the owning cell: compute_cell_geometry_
        // oriented() negates it for the neighbour, so one global winding per
        // face is enough to give every cell a positive signed volume.
        const Vec3 n = (B - A).cross(C - A);
        const Vec3 face_centre = (A + B + C) * (1.0 / 3.0);
        if (n.dot(face_centre - owner_centre) < 0.0)
            std::swap(b, c);
        const std::size_t id = mesh.faces().n_faces();
        mesh.faces().push_face({a,b,c});
        face_ids.emplace(key, id);
        face_cells.emplace_back();
        return id;
    };

    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i) {
                std::array<std::size_t,8> corner{};
                const std::array<std::array<unsigned,3>,8> corners = {{
                    {{0,0,0}}, {{1,0,0}}, {{1,1,0}}, {{0,1,0}},
                    {{0,0,1}}, {{1,0,1}}, {{1,1,1}}, {{0,1,1}}
                }};
                for (unsigned q = 0; q < 8; ++q)
                    corner[q] = vid(i+corners[q][0], j+corners[q][1], k+corners[q][2]);

                for (const auto& tet : tets) {
                    const std::array<std::size_t,4> v = {{
                        corner[tet[0]], corner[tet[1]],
                        corner[tet[2]], corner[tet[3]]
                    }};
                    const Vec3 tet_centre =
                        (Vec3{mesh.points().x(v[0]), mesh.points().y(v[0]), mesh.points().z(v[0])} +
                         Vec3{mesh.points().x(v[1]), mesh.points().y(v[1]), mesh.points().z(v[1])} +
                         Vec3{mesh.points().x(v[2]), mesh.points().y(v[2]), mesh.points().z(v[2])} +
                         Vec3{mesh.points().x(v[3]), mesh.points().y(v[3]), mesh.points().z(v[3])}) *
                        0.25;
                    const std::array<std::size_t,4> faces = {{
                        get_face(v[0],v[1],v[2],tet_centre),
                        get_face(v[0],v[1],v[3],tet_centre),
                        get_face(v[0],v[2],v[3],tet_centre),
                        get_face(v[1],v[2],v[3],tet_centre)
                    }};
                    const std::size_t c = cell_faces.size();
                    cell_faces.push_back({faces.begin(), faces.end()});
                    for (const auto f : faces) face_cells[f].push_back(c);
                }
            }

    mesh.ownership().resize(mesh.n_faces());
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (face_cells[f].empty() || face_cells[f].size() > 2)
            throw std::runtime_error("invalid tetra face adjacency");
        mesh.ownership().set_owner(f, face_cells[f][0]);
        mesh.ownership().set_neighbour(
            f, face_cells[f].size() == 1
                ? FaceOwnership::BOUNDARY
                : static_cast<std::int64_t>(face_cells[f][1]));
    }
    for (const auto& faces : cell_faces) mesh.cells().push_cell(faces);

    add_boundary_patches(mesh, "boundary");

    TetraGrid grid{std::move(mesh), {}};
    grid.geometry = make_geometry_cache(grid.mesh);
    return grid;
}

double cell_bbox_aspect_ratio(const Mesh& mesh, std::size_t cell)
{
    const auto offset = mesh.cells().cell_offset(cell);
    const auto size = mesh.cells().cell_size(cell);
    std::set<std::size_t> vertices;
    for (std::size_t k = 0; k < size; ++k) {
        const auto f = mesh.cells().faces_data()[offset + k];
        const auto fo = mesh.faces().face_offset(f);
        const auto fn = mesh.faces().face_size(f);
        for (std::size_t q = 0; q < fn; ++q)
            vertices.insert(mesh.faces().vertices_data()[fo + q]);
    }

    Vec3 lo{std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()};
    Vec3 hi{-std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity()};
    for (const auto v : vertices) {
        lo.x = std::min(lo.x, mesh.points().x(v));
        lo.y = std::min(lo.y, mesh.points().y(v));
        lo.z = std::min(lo.z, mesh.points().z(v));
        hi.x = std::max(hi.x, mesh.points().x(v));
        hi.y = std::max(hi.y, mesh.points().y(v));
        hi.z = std::max(hi.z, mesh.points().z(v));
    }

    const double ex = hi.x - lo.x;
    const double ey = hi.y - lo.y;
    const double ez = hi.z - lo.z;
    const double emin = std::min({ex, ey, ez});
    const double emax = std::max({ex, ey, ez});
    if (!(emin > 0.0) || !std::isfinite(emin))
        return std::numeric_limits<double>::infinity();
    return emax / emin;
}

SparseMatrix assemble_two_point_diffusion(const Grid& grid)
{
    const std::size_t n = grid.mesh.n_cells();
    SparseMatrix matrix(n, n);
    std::vector<double> diagonal(n, 0.0);

    for (std::size_t f = 0; f < grid.mesh.n_faces(); ++f) {
        const std::size_t owner = grid.mesh.ownership().owner(f);
        const auto neighbour_raw = grid.mesh.ownership().neighbour(f);
        const double area = grid.geometry.face_areas[f];

        if (!(area > 0.0)) throw std::runtime_error("non-positive face area");
        if (neighbour_raw >= 0) {
            const std::size_t neighbour = static_cast<std::size_t>(neighbour_raw);
            const double d = (grid.geometry.cell_centres[neighbour] -
                              grid.geometry.cell_centres[owner]).mag();
            if (!(d > 0.0)) throw std::runtime_error("zero owner-neighbour distance");
            const double w = area / d;
            diagonal[owner] += w;
            diagonal[neighbour] += w;
            matrix.push_back(owner, neighbour, -w);
            matrix.push_back(neighbour, owner, -w);
        } else {
            const double d = (grid.geometry.face_centres[f] -
                              grid.geometry.cell_centres[owner]).mag();
            if (!(d > 0.0)) throw std::runtime_error("zero boundary distance");
            diagonal[owner] += area / d;
        }
    }

    for (std::size_t c = 0; c < n; ++c) matrix.push_back(c, c, diagonal[c]);
    matrix.finalize();
    return matrix;
}

double relative_true_residual(const SparseMatrix& matrix,
                              const Vector& x,
                              const Vector& rhs)
{
    const auto ax = matrix.matvec(x);
    double r2 = 0.0;
    double b2 = 0.0;
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        const double r = rhs(i) - ax[i];
        r2 += r * r;
        b2 += rhs(i) * rhs(i);
    }
    return std::sqrt(r2 / std::max(b2, 1e-300));
}

void check_solver_robustness(const QualityCase& qc)
{
    const Grid grid = make_hex_grid(qc.n, qc.shear, qc.stretch);
    const auto matrix = assemble_two_point_diffusion(grid);

    Vector exact(matrix.n_rows());
    for (std::size_t i = 0; i < exact.size(); ++i)
        exact(i) = std::sin(0.37 * static_cast<double>(i + 1)) +
                   0.1 * std::cos(0.11 * static_cast<double>(i + 3));

    const auto rhs_values = matrix.matvec(exact);
    Vector rhs(rhs_values.size());
    for (std::size_t i = 0; i < rhs.size(); ++i) rhs(i) = rhs_values[i];

    Vector solution(matrix.n_rows(), 0.0);
    const auto result = solve_cg(matrix, rhs, solution, 5000, 1e-10);
    const double true_residual = relative_true_residual(matrix, solution, rhs);

    EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
    EXPECT_TRUE(std::isfinite(true_residual));
    EXPECT_TRUE(true_residual < 1e-8);

    std::cout << "N10_SOLVER case=" << qc.name
              << " cells=" << matrix.n_rows()
              << " iterations=" << result.iterations
              << " reported_relative_residual=" << result.residual_relative
              << " true_relative_residual=" << true_residual << "\n";
}

void check_hex_quality_ladder()
{
    const std::array<QualityCase,4> cases = {{
        {"orthogonal", 0.0, 1.0, 6},
        {"mild_skew", 0.20, 1.0, 6},
        {"strong_skew", 0.80, 1.0, 6},
        {"high_aspect", 0.0, 20.0, 6}
    }};

    double previous_nonorth = -1.0;
    for (const auto& qc : cases) {
        const Grid grid = make_hex_grid(qc.n, qc.shear, qc.stretch);
        const auto report = validate_mesh(grid.mesh);

        EXPECT_TRUE(report.ok);
        EXPECT_TRUE(std::isfinite(report.max_skewness));
        EXPECT_TRUE(std::isfinite(report.max_non_orthogonality_deg));
        EXPECT_TRUE(report.min_cell_volume > 0.0);
        EXPECT_TRUE(report.max_cell_volume >= report.min_cell_volume);

        double max_aspect = 0.0;
        for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c)
            max_aspect = std::max(max_aspect, cell_bbox_aspect_ratio(grid.mesh, c));
        EXPECT_TRUE(std::isfinite(max_aspect));

        if (qc.name == std::string("orthogonal")) {
            EXPECT_NEAR(report.max_non_orthogonality_deg, 0.0, 1e-12);
            EXPECT_NEAR(report.max_skewness, 0.0, 1e-12);
            EXPECT_NEAR(max_aspect, 1.0, 1e-12);
        }
        if (qc.name == std::string("mild_skew")) {
            EXPECT_TRUE(report.max_skewness > 0.0);
            EXPECT_TRUE(report.max_non_orthogonality_deg > previous_nonorth);
        }
        if (qc.name == std::string("strong_skew")) {
            EXPECT_TRUE(report.max_skewness > 0.0);
            EXPECT_TRUE(report.max_non_orthogonality_deg >
                        1.0);
        }
        if (qc.name == std::string("high_aspect"))
            EXPECT_TRUE(max_aspect >= 20.0 - 1e-12);

        previous_nonorth = report.max_non_orthogonality_deg;

        std::cout << "N10_QUALITY case=" << qc.name
                  << " cells=" << grid.mesh.n_cells()
                  << " max_skewness=" << report.max_skewness
                  << " max_nonorth_deg=" << report.max_non_orthogonality_deg
                  << " max_aspect=" << max_aspect
                  << " min_volume=" << report.min_cell_volume
                  << " max_volume=" << report.max_cell_volume << "\n";

        check_solver_robustness(qc);
    }
}

void check_polyhedral_quality()
{
    const auto grid = make_tetra_grid(4);
    const auto report = validate_mesh(grid.mesh);

    EXPECT_TRUE(report.ok);
    EXPECT_TRUE(std::isfinite(report.max_skewness));
    EXPECT_TRUE(std::isfinite(report.max_non_orthogonality_deg));
    EXPECT_TRUE(report.min_cell_volume > 0.0);

    double max_aspect = 0.0;
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c)
        max_aspect = std::max(max_aspect, cell_bbox_aspect_ratio(grid.mesh, c));

    EXPECT_TRUE(std::isfinite(max_aspect));
    EXPECT_TRUE(max_aspect >= 1.0);

    ScalarCellField phi(grid.mesh.n_cells(), "phi", "1", 1);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c)
        phi(c) = 2.0 * grid.geometry.cell_centres[c].x
               - 3.0 * grid.geometry.cell_centres[c].y
               + 0.5 * grid.geometry.cell_centres[c].z;

const auto grad = compute_gradient_least_squares(phi, grid.mesh);
    // Linear-exactness is only meaningful on interior cells. Boundary tets keep
    // a one-sided stencil (boundary faces are excluded), so their least-squares
    // system is rank-deficient and cannot reproduce the gradient.
    const std::size_t n = 4;
    double max_error = 0.0;
    std::size_t checked = 0;
    for (std::size_t k = 1; k + 1 < n; ++k)
        for (std::size_t j = 1; j + 1 < n; ++j)
            for (std::size_t i = 1; i + 1 < n; ++i) {
                const std::size_t base = 6 * (i + n * (j + n * k));
                for (std::size_t t = 0; t < 6; ++t) {
                    const std::size_t c = base + t;
                    const Vec3 got{grad(c,0), grad(c,1), grad(c,2)};
                    max_error = std::max(max_error,
                                         (got - Vec3{2.0,-3.0,0.5}).mag());
                    ++checked;
                }
            }
    EXPECT_TRUE(checked > 0);
    EXPECT_TRUE(std::isfinite(max_error));
    EXPECT_TRUE(max_error < 1e-8);

    std::cout << "N10_POLYHEDRAL cells=" << grid.mesh.n_cells()
              << " max_skewness=" << report.max_skewness
              << " max_nonorth_deg=" << report.max_non_orthogonality_deg
              << " max_aspect=" << max_aspect
              << " linear_gradient_Linf=" << max_error << "\n";
}

void check_near_degenerate_valid()
{
    const Grid grid = make_hex_grid(3, 0.0, 1.0e-2);
    const auto report = validate_mesh(grid.mesh);
    EXPECT_TRUE(report.ok);
    EXPECT_TRUE(report.min_cell_volume > 0.0);

    double max_aspect = 0.0;
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c)
        max_aspect = std::max(max_aspect, cell_bbox_aspect_ratio(grid.mesh, c));
    EXPECT_TRUE(max_aspect >= 100.0 - 1e-10);

    std::cout << "N10_NEAR_DEGENERATE_VALID max_aspect=" << max_aspect
              << " min_volume=" << report.min_cell_volume << "\n";
}

void check_invalid_rejection()
{
    Grid grid = make_hex_grid(2, 0.0, 1.0);
    const auto face = grid.mesh.cells().faces_data()[0];
    const auto offset = grid.mesh.faces().face_offset(face);
    const auto size = grid.mesh.faces().face_size(face);
    const auto first = grid.mesh.faces().vertices_data()[offset];
    const double x = grid.mesh.points().x(first);
    const double y = grid.mesh.points().y(first);
    const double z = grid.mesh.points().z(first);

    // Collapse every vertex of the face onto a single point. Merging only two
    // of them would leave a non-degenerate triangle, which validate_mesh()
    // accepts, so the whole face has to be collapsed to get a zero-area face.
    for (std::size_t q = 0; q < size; ++q) {
        const auto v = grid.mesh.faces().vertices_data()[offset + q];
        grid.mesh.points().set(v, x, y, z);
    }

    const auto report = validate_mesh(grid.mesh);
    EXPECT_FALSE(report.ok);
    EXPECT_TRUE(!report.errors.empty());
    std::cout << "N10_INVALID rejected_errors=" << report.errors.size() << "\n";
}


} // namespace

int main()
{
    std::cout << std::setprecision(12);

    run_case("n10_quality_ladder", check_hex_quality_ladder);
    run_case("n10_polyhedral_quality", check_polyhedral_quality);
    run_case("n10_near_degenerate_valid", check_near_degenerate_valid);
    run_case("n10_invalid_rejection", check_invalid_rejection);

    return run_all();
}
