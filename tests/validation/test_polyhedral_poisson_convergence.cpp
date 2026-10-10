// Issue #461 — N10 solution-level convergence of the corrected polyhedral
// Laplacian.
//
// The polyhedral Laplacian campaign enforces constant/linear consistency of
// the corrected operator and measures its smooth-field interior residual
// WITHOUT a threshold: the linear-exact but non-symmetric corrected operator
// is not quadratic-exact, and a strict residual-decrease gate was measured to
// FAIL on the Kuhn tetrahedral family (negative result, PR #845). That result
// closes the pointwise-residual claim only; it does not answer whether the
// SOLUTION of the discrete Poisson problem converges under refinement.
//
// This test closes the remaining solution-level gap. It solves
//     Delta phi = -3 pi^2 phi,   phi = sin(pi x) sin(pi y) sin(pi z)
// on Kuhn-tetrahedralised unit cubes refined at n = 4/8/16 interior cubes,
// with the outer tet shell carrying exact Dirichlet ghost values.
//
// Fail-closed design:
//   * The interior-interior block of the discrete operator is assembled from
//     black-box applications of the production compute_laplacian on
//     indicator fields; the operator is never re-implemented here. A
//     distance-5 cube colouring guarantees that two interior cells of the
//     same colour are more than the two-hop dependency support of the
//     corrected operator apart (the least-squares cell gradient reads one
//     face-hop ring around each one-hop cell), so per-colour responses
//     cannot contaminate each other.
//   * The assembly is GATED: for a deterministic pseudo-random field that
//     vanishes on the ghost shell, A*psi must reproduce compute_laplacian(psi)
//     on every interior cell up to floating-point round-off. A wrong
//     colouring or an incomplete dependency support cannot pass this gate.
//   * The right-hand side is built through the production operator from the
//     exact field, so the solved system is exactly the discrete problem the
//     production operator defines with exact Dirichlet ghosts.
//   * The system is solved with the production BiCGStab (the corrected
//     operator is non-symmetric); CONVERGED status and an independently
//     recomputed true relative residual <= 1e-8 are required (same solver
//     gate as the N10 difficult-mesh campaign).
//   * ACCEPTANCE: the volume-weighted L2 solution error must strictly
//     decrease at every refinement level. Observed orders are reported as
//     diagnostics only; no accuracy order is claimed. If the strict decrease
//     fails, that is an honest negative result which blocks the N10 campaign;
//     no tolerance of any existing test is changed.

#include "cfdx/core/field/field.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/geometry/mesh_validator.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/laplacian.h"
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

constexpr std::size_t kNoRow = static_cast<std::size_t>(-1);

double manufactured(const Vec3& p)
{
    const double pi = std::acos(-1.0);
    return std::sin(pi * p.x) * std::sin(pi * p.y) * std::sin(pi * p.z);
}

double manufactured_source(const Vec3& p)
{
    const double pi = std::acos(-1.0);
    return -3.0 * pi * pi * manufactured(p);
}

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

std::string tri_key(std::size_t a, std::size_t b, std::size_t c)
{
    std::array<std::size_t, 3> v{a, b, c};
    std::sort(v.begin(), v.end());
    return std::to_string(v[0]) + "/" + std::to_string(v[1]) + "/" + std::to_string(v[2]);
}

// validate_mesh() requires every boundary face to belong to a patch.
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

struct PoissonGrid {
    Mesh mesh;
    GeometryCache geometry;
    std::size_t n = 0;       // interior cubes per direction
    std::size_t ncubes = 0;  // cubes per direction (n + 2)
    std::vector<std::size_t> cell_cube;
    std::vector<unsigned> cell_tet;
    std::vector<std::size_t> row_of_cell;
    std::vector<std::size_t> cell_of_row;
};

// Kuhn-tetrahedralised unit cube with n interior cubes per direction and a
// one-cube ghost shell. Cell ids follow the generation order
// 6 * (i + ncubes * (j + ncubes * k)) + tet, and the winding of every face is
// stored away from its first (owning) cell so all signed volumes are
// positive (same correction as the N10 difficult-mesh generator).
PoissonGrid make_poisson_grid(std::size_t n)
{
    if (n < 2) throw std::invalid_argument("make_poisson_grid: n must be >= 2");
    const std::size_t N = n + 2;
    const std::size_t nv = N + 1;
    const auto vid = [nv](std::size_t i, std::size_t j, std::size_t k) {
        return i + nv * (j + nv * k);
    };
    const std::array<std::array<unsigned, 4>, 6> tets = {{
        {{0,1,2,6}}, {{0,1,5,6}}, {{0,3,2,6}},
        {{0,3,7,6}}, {{0,4,5,6}}, {{0,4,7,6}}
    }};

    Mesh mesh;
    mesh.points().resize(nv * nv * nv);
    for (std::size_t k = 0; k <= N; ++k)
        for (std::size_t j = 0; j <= N; ++j)
            for (std::size_t i = 0; i <= N; ++i)
                mesh.points().set(vid(i, j, k),
                    static_cast<double>(i) / static_cast<double>(N),
                    static_cast<double>(j) / static_cast<double>(N),
                    static_cast<double>(k) / static_cast<double>(N));

    std::unordered_map<std::string, std::size_t> face_ids;
    std::vector<std::vector<std::size_t>> cell_faces;
    std::vector<std::vector<std::size_t>> face_cells;
    std::vector<std::size_t> cell_cube;
    std::vector<unsigned> cell_tet;

    auto get_face = [&](std::size_t a, std::size_t b, std::size_t c,
                        const Vec3& owner_centre) {
        const std::string key = tri_key(a, b, c);
        auto it = face_ids.find(key);
        if (it != face_ids.end()) return it->second;
        const Vec3 A{mesh.points().x(a), mesh.points().y(a), mesh.points().z(a)};
        const Vec3 B{mesh.points().x(b), mesh.points().y(b), mesh.points().z(b)};
        const Vec3 C{mesh.points().x(c), mesh.points().y(c), mesh.points().z(c)};
        // Half of the Kuhn subdivision tets are left-handed, so the vertex
        // order of a face is not by itself an outward orientation: store the
        // winding pointing away from the owning cell.
        const Vec3 fn = (B - A).cross(C - A);
        const Vec3 face_centre = (A + B + C) * (1.0 / 3.0);
        if (fn.dot(face_centre - owner_centre) < 0.0)
            std::swap(b, c);
        const std::size_t id = face_cells.size();
        mesh.faces().push_face({a, b, c});
        face_ids.emplace(key, id);
        face_cells.emplace_back();
        return id;
    };

    for (std::size_t k = 0; k < N; ++k)
        for (std::size_t j = 0; j < N; ++j)
            for (std::size_t i = 0; i < N; ++i) {
                std::array<std::size_t, 8> corner{};
                const std::array<std::array<unsigned, 3>, 8> corners = {{
                    {{0,0,0}}, {{1,0,0}}, {{1,1,0}}, {{0,1,0}},
                    {{0,0,1}}, {{1,0,1}}, {{1,1,1}}, {{0,1,1}}
                }};
                for (unsigned q = 0; q < 8; ++q)
                    corner[q] = vid(i + corners[q][0], j + corners[q][1], k + corners[q][2]);
                for (unsigned t = 0; t < tets.size(); ++t) {
                    const auto& tet = tets[t];
                    const std::array<std::size_t, 4> v = {{
                        corner[tet[0]], corner[tet[1]],
                        corner[tet[2]], corner[tet[3]]
                    }};
                    const Vec3 tet_centre =
                        (Vec3{mesh.points().x(v[0]), mesh.points().y(v[0]), mesh.points().z(v[0])} +
                         Vec3{mesh.points().x(v[1]), mesh.points().y(v[1]), mesh.points().z(v[1])} +
                         Vec3{mesh.points().x(v[2]), mesh.points().y(v[2]), mesh.points().z(v[2])} +
                         Vec3{mesh.points().x(v[3]), mesh.points().y(v[3]), mesh.points().z(v[3])}) * 0.25;
                    const std::array<std::size_t, 4> faces = {{
                        get_face(v[0], v[1], v[2], tet_centre),
                        get_face(v[0], v[1], v[3], tet_centre),
                        get_face(v[0], v[2], v[3], tet_centre),
                        get_face(v[1], v[2], v[3], tet_centre)
                    }};
                    const std::size_t c = cell_faces.size();
                    cell_faces.push_back({faces.begin(), faces.end()});
                    cell_cube.push_back(i + N * (j + N * k));
                    cell_tet.push_back(t);
                    for (const std::size_t f : faces) face_cells[f].push_back(c);
                }
            }

    mesh.ownership().resize(mesh.n_faces());
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (face_cells[f].empty() || face_cells[f].size() > 2)
            throw std::runtime_error("invalid tetra face adjacency");
        mesh.ownership().set_owner(f, face_cells[f][0]);
        mesh.ownership().set_neighbour(f, face_cells[f].size() == 1
            ? FaceOwnership::BOUNDARY
            : static_cast<std::int64_t>(face_cells[f][1]));
    }
    for (const auto& faces : cell_faces) mesh.cells().push_cell(faces);

    add_boundary_patches(mesh, "boundary");

    PoissonGrid grid;
    grid.mesh = std::move(mesh);
    grid.geometry = make_geometry_cache(grid.mesh);
    grid.n = n;
    grid.ncubes = N;
    grid.cell_cube = std::move(cell_cube);
    grid.cell_tet = std::move(cell_tet);

    const std::size_t n_cells = grid.mesh.n_cells();
    grid.row_of_cell.assign(n_cells, kNoRow);
    std::size_t row = 0;
    for (std::size_t k = 1; k <= n; ++k)
        for (std::size_t j = 1; j <= n; ++j)
            for (std::size_t i = 1; i <= n; ++i)
                for (unsigned t = 0; t < 6; ++t) {
                    const std::size_t c = 6 * (i + N * (j + N * k)) + t;
                    grid.row_of_cell[c] = row++;
                }
    grid.cell_of_row.resize(row);
    for (std::size_t c = 0; c < n_cells; ++c)
        if (grid.row_of_cell[c] != kNoRow)
            grid.cell_of_row[grid.row_of_cell[c]] = c;
    return grid;
}

std::vector<std::vector<std::size_t>> face_neighbour_lists(const Mesh& mesh)
{
    const std::size_t n_cells = mesh.n_cells();
    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    const FaceOwnership& own = mesh.ownership();
    std::vector<std::vector<std::size_t>> ring(n_cells);
    for (std::size_t c = 0; c < n_cells; ++c) {
        auto& list = ring[c];
        for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            const std::size_t owner = own.owner(f);
            const std::int64_t nraw = own.neighbour(f);
            std::size_t nb = n_cells;
            if (owner == c) {
                if (nraw >= 0) nb = static_cast<std::size_t>(nraw);
            } else {
                nb = owner;
            }
            if (nb >= n_cells || nb == c) continue;
            if (std::find(list.begin(), list.end(), nb) == list.end())
                list.push_back(nb);
        }
    }
    return ring;
}

// Dependency support of one interior unknown: every cell within two face
// hops. The corrected Laplacian at cell c reads the cell gradients of c and
// of its one-hop face neighbours, and each least-squares gradient reads the
// one-hop ring of its cell, so phi influences the operator output at c
// within two hops exactly. The assembly gate below fails closed if this
// support model ever disagrees with the production operator.
std::vector<std::vector<std::size_t>> dependency_supports(
    const PoissonGrid& grid,
    const std::vector<std::vector<std::size_t>>& ring)
{
    const std::size_t n_cells = grid.mesh.n_cells();
    const std::size_t n_rows = grid.cell_of_row.size();
    std::vector<std::vector<std::size_t>> supports(n_rows);
    std::vector<std::uint8_t> marked(n_cells, 0);
    std::vector<std::size_t> frontier;
    std::vector<std::size_t> next;
    for (std::size_t r = 0; r < n_rows; ++r) {
        const std::size_t d = grid.cell_of_row[r];
        marked[d] = 1;
        supports[r].push_back(d);
        frontier.clear();
        frontier.push_back(d);
        for (int depth = 0; depth < 2; ++depth) {
            next.clear();
            for (const std::size_t c : frontier)
                for (const std::size_t nb : ring[c]) {
                    if (marked[nb]) continue;
                    marked[nb] = 1;
                    supports[r].push_back(nb);
                    next.push_back(nb);
                }
            frontier.swap(next);
        }
        for (const std::size_t c : supports[r]) marked[c] = 0;
    }
    return supports;
}

// Interior-interior block of the discrete corrected Laplacian, assembled
// from black-box operator applications on indicator fields. Distance-5
// colouring: two interior cells share a colour only when their cubes differ
// by at least five in one direction, hence strictly more than the two-hop
// support apart; the tet index is part of the colour so the six tets of one
// cube never share a colour.
SparseMatrix assemble_operator_block(
    const PoissonGrid& grid,
    const std::vector<std::vector<std::size_t>>& supports)
{
    const std::size_t n_cells = grid.mesh.n_cells();
    const std::size_t n_rows = grid.cell_of_row.size();
    SparseMatrix matrix(n_rows, n_rows);

    std::vector<std::vector<std::size_t>> colour(750);
    for (std::size_t r = 0; r < n_rows; ++r) {
        const std::size_t c = grid.cell_of_row[r];
        const std::size_t cube = grid.cell_cube[c];
        const std::size_t ci = cube % grid.ncubes;
        const std::size_t cj = (cube / grid.ncubes) % grid.ncubes;
        const std::size_t ck = cube / (grid.ncubes * grid.ncubes);
        const std::size_t id =
            (((ci % 5) * 5 + (cj % 5)) * 5 + (ck % 5)) * 6 + grid.cell_tet[c];
        colour[id].push_back(r);
    }

    Field<double, Location::CELL> indicator(n_cells, "indicator", "1", 1);
    for (const auto& rows : colour) {
        if (rows.empty()) continue;
        for (std::size_t c = 0; c < n_cells; ++c) indicator(c) = 0.0;
        for (const std::size_t r : rows) indicator(grid.cell_of_row[r]) = 1.0;
        const auto response = compute_laplacian(
            indicator, grid.mesh, grid.geometry,
            LaplacianScheme::CORRECTED, 0.5, GradientScheme::LEAST_SQUARES);
        for (const std::size_t r : rows)
            for (const std::size_t c : supports[r]) {
                const std::size_t row = grid.row_of_cell[c];
                if (row == kNoRow) continue;
                matrix.push_back(row, r, response(c));
            }
    }
    matrix.finalize();
    return matrix;
}

struct LevelResult {
    std::size_t n = 0;
    double l2 = 0.0;
};

LevelResult run_level(std::size_t n)
{
    const PoissonGrid grid = make_poisson_grid(n);
    const auto report = validate_mesh(grid.mesh);
    require(report.ok, "the generated Poisson grid must be a valid mesh");
    const auto ring = face_neighbour_lists(grid.mesh);
    const auto supports = dependency_supports(grid, ring);
    const SparseMatrix matrix = assemble_operator_block(grid, supports);
    const std::size_t n_rows = grid.cell_of_row.size();
    require(matrix.n_rows() == n_rows && matrix.n_cols() == n_rows,
            "the assembled block must be square over the interior unknowns");

    // Assembly gate: A*psi must reproduce the production operator on a
    // deterministic field that vanishes on the ghost shell.
    Field<double, Location::CELL> psi(grid.mesh.n_cells(), "psi", "1", 1);
    std::uint32_t lcg = 0x9e3779b9u;
    for (std::size_t r = 0; r < n_rows; ++r) {
        lcg = lcg * 1664525u + 1013904223u;
        psi(grid.cell_of_row[r]) = static_cast<double>(lcg >> 8) / 8388608.0 - 1.0;
    }
    const auto direct = compute_laplacian(
        psi, grid.mesh, grid.geometry,
        LaplacianScheme::CORRECTED, 0.5, GradientScheme::LEAST_SQUARES);
    Vector psi_int(n_rows, 0.0);
    for (std::size_t r = 0; r < n_rows; ++r)
        psi_int(r) = psi(grid.cell_of_row[r]);
    const auto assembled = matrix.matvec(psi_int);
    double max_diff = 0.0;
    double max_direct = 0.0;
    for (std::size_t r = 0; r < n_rows; ++r) {
        const std::size_t c = grid.cell_of_row[r];
        max_diff = std::max(max_diff, std::abs(direct(c) - assembled[r]));
        max_direct = std::max(max_direct, std::abs(direct(c)));
    }
    const double assembly_rel = max_diff / std::max(1.0, max_direct);
    require(std::isfinite(assembly_rel) && assembly_rel <= 1e-9,
            "assembly self-check: A*psi must reproduce compute_laplacian(psi) "
            "on the interior cells");

    // Right-hand side built through the production operator from the exact
    // field: b = f - L(phi_exact) + A*phi_exact_interior.
    Field<double, Location::CELL> exact(grid.mesh.n_cells(), "phi_exact", "1", 1);
    for (std::size_t c = 0; c < grid.mesh.n_cells(); ++c)
        exact(c) = manufactured(grid.geometry.cell_centres[c]);
    const auto lap_exact = compute_laplacian(
        exact, grid.mesh, grid.geometry,
        LaplacianScheme::CORRECTED, 0.5, GradientScheme::LEAST_SQUARES);
    Vector exact_int(n_rows, 0.0);
    for (std::size_t r = 0; r < n_rows; ++r)
        exact_int(r) = exact(grid.cell_of_row[r]);
    const auto interior_part = matrix.matvec(exact_int);
    Vector b(n_rows, 0.0);
    for (std::size_t r = 0; r < n_rows; ++r) {
        const std::size_t c = grid.cell_of_row[r];
        b(r) = manufactured_source(grid.geometry.cell_centres[c])
             - lap_exact(c) + interior_part[r];
    }

    // Production non-symmetric Krylov solve.
    Vector x(n_rows, 0.0);
    const auto result = solve_bicgstab(matrix, b, x, 20000, 1e-10);
    require(result.status == SolverStatus::CONVERGED,
            "BiCGStab must converge on the assembled Poisson system");
    const auto ax = matrix.matvec(x);
    double r2 = 0.0;
    double b2 = 0.0;
    for (std::size_t r = 0; r < n_rows; ++r) {
        const double residual = b(r) - ax[r];
        r2 += residual * residual;
        b2 += b(r) * b(r);
    }
    const double true_relative_residual = std::sqrt(r2 / std::max(b2, 1e-300));
    require(std::isfinite(true_relative_residual) && true_relative_residual <= 1e-8,
            "the true relative residual must be finite and below 1e-8");

    // Solution error against the manufactured field at interior centroids.
    std::vector<double> got;
    std::vector<double> want;
    std::vector<double> weights;
    got.reserve(n_rows);
    want.reserve(n_rows);
    weights.reserve(n_rows);
    for (std::size_t r = 0; r < n_rows; ++r) {
        const std::size_t c = grid.cell_of_row[r];
        got.push_back(x(r));
        want.push_back(exact(c));
        weights.push_back(grid.geometry.cell_volumes[c]);
    }
    const auto metrics = error_norms(got, want, weights);
    require(std::isfinite(metrics.l2) && std::isfinite(metrics.linf),
            "the solution error norms must be finite");

    std::cout << "N10_POISSON n=" << n
              << " cells=" << grid.mesh.n_cells()
              << " unknowns=" << n_rows
              << " iterations=" << result.iterations
              << " l2=" << metrics.l2
              << " linf=" << metrics.linf
              << " assembly_max_rel_diff=" << assembly_rel
              << " true_relative_residual=" << true_relative_residual << "\n";

    LevelResult out;
    out.n = n;
    out.l2 = metrics.l2;
    return out;
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(12);
        const std::array<std::size_t, 3> levels = {{4u, 8u, 16u}};
        std::vector<LevelResult> results;
        for (const std::size_t n : levels) results.push_back(run_level(n));

        for (std::size_t idx = 1; idx < results.size(); ++idx) {
            const LevelResult& coarse = results[idx - 1];
            const LevelResult& fine = results[idx];
            const double order = observed_order(coarse.l2, fine.l2);
            std::cout << "N10_POISSON_ORDER from_n=" << coarse.n
                      << " to_n=" << fine.n
                      << " observed_order=" << order << "\n";
            if (coarse.l2 > 0.0 && fine.l2 > 0.0)
                require(std::isfinite(order),
                        "the observed solution order must be finite");
            require(fine.l2 < coarse.l2,
                    "the solution L2 error must strictly decrease under "
                    "refinement (n=" + std::to_string(coarse.n) + " -> n=" +
                    std::to_string(fine.n) + ")");
        }
        std::cout << "POLYHEDRAL_POISSON_CONVERGENCE: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "POLYHEDRAL_POISSON_CONVERGENCE: FAIL: " << e.what() << "\n";
        return 1;
    }
}
