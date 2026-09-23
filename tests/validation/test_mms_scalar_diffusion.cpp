// Phase 3.7 — Manufactured-solution verification for scalar diffusion.

#include "cfdx/core/solvers/scalar_diffusion.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include <unordered_map>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

Mesh make_channel(std::size_t n)
{
    if (n == 0) throw std::invalid_argument("make_channel: n must be positive");

    Mesh m;
    m.points().resize(8 * n);
    const double dy = 1.0 / static_cast<double>(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double y0 = dy * static_cast<double>(i);
        const double y1 = dy * static_cast<double>(i + 1);
        const std::size_t b = 8 * i;
        const double p[8][3] = {
            {0.0,y0,0.0}, {1.0,y0,0.0}, {1.0,y1,0.0}, {0.0,y1,0.0},
            {0.0,y0,1.0}, {1.0,y0,1.0}, {1.0,y1,1.0}, {0.0,y1,1.0}
        };
        for (std::size_t j = 0; j < 8; ++j)
            m.points().set(b + j, p[j][0], p[j][1], p[j][2]);
    }

    std::vector<std::size_t> bottom, top, x0, x1, z0, z1, internal;
    auto add_face = [&](std::initializer_list<std::size_t> vertices) {
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face(std::vector<FaceIndex>(vertices.begin(), vertices.end()));
        return id;
    };

    bottom.push_back(add_face({0,1,5,4}));
    top.push_back(add_face({8*(n-1)+3,8*(n-1)+7,8*(n-1)+6,8*(n-1)+2}));
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t b = 8 * i;
        x0.push_back(add_face({b,b+4,b+7,b+3}));
        x1.push_back(add_face({b+1,b+2,b+6,b+5}));
        z0.push_back(add_face({b,b+3,b+2,b+1}));
        z1.push_back(add_face({b+4,b+5,b+6,b+7}));
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        const std::size_t b = 8 * i;
        internal.push_back(add_face({b+3,b+7,b+6,b+2}));
    }

    m.ownership().resize(m.n_faces());

    // Ownership is a property of each face, not of each cell incidence.
    // Internal face i is shared by cells i and i+1, so its owner/neighbour
    // must be assigned once. Reassigning it while visiting the neighbour
    // would overwrite the canonical owner and corrupt the topology.
    m.ownership().set_owner(bottom[0], 0);
    m.ownership().set_neighbour(bottom[0], FaceOwnership::BOUNDARY);
    m.ownership().set_owner(top[0], n - 1);
    m.ownership().set_neighbour(top[0], FaceOwnership::BOUNDARY);

    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t boundary_faces[] = {
            x0[i], x1[i], z0[i], z1[i]
        };
        for (const auto f : boundary_faces) {
            m.ownership().set_owner(f, i);
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        }
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        m.ownership().set_owner(internal[i], i);
        m.ownership().set_neighbour(internal[i], static_cast<std::int64_t>(i + 1));
    }

    std::vector<std::vector<std::size_t>> cell_faces(n);
    for (std::size_t i = 0; i < n; ++i) {
        cell_faces[i] = {
            i == 0 ? bottom[0] : internal[i-1],
            i + 1 == n ? top[0] : internal[i],
            x0[i], x1[i], z0[i], z1[i]
        };
    }
    for (const auto& faces : cell_faces) m.cells().push_cell(faces);

    auto add_patch = [&](const char* name, const std::vector<std::size_t>& faces) {
        Patch p;
        p.name = name;
        p.type = PatchType::WALL;
        p.face_ids = faces;
        m.boundary().add_patch(p);
    };
    add_patch("bottom", bottom);
    add_patch("top", top);
    add_patch("x0", x0);
    add_patch("x1", x1);
    add_patch("z0", z0);
    add_patch("z1", z1);
    return m;
}

struct Error {
    double l2;
    double linf;
};

Error solve_mms(std::size_t n)
{
    const Mesh mesh = make_channel(n);
    const auto geometry = make_geometry_cache(mesh);
    constexpr double gamma = 1.0;

    // Manufactured solution: phi(y) = y^2.
    // Therefore -div(gamma grad(phi)) = -2 gamma.
    // The solver accepts a volumetric source density and integrates it
    // internally as S_c * V_c.
    std::vector<double> source(n, 0.0);
    for (std::size_t c = 0; c < n; ++c)
        source[c] = -2.0 * gamma;
    DirichletBoundary bc;
    bc.face_values.assign(
        mesh.n_faces(), std::numeric_limits<double>::quiet_NaN());
    bc.face_values[0] = 0.0;
    bc.face_values[1] = 1.0;

    ScalarDiffusionConfig cfg;
    cfg.diffusivity = gamma;
    cfg.max_iterations = 5000;
    cfg.tolerance = 1e-13;
    const auto result = solve_poisson_dirichlet(mesh, bc, source, cfg);
    if (result.linear_result.status != SolverStatus::CONVERGED)
        throw std::runtime_error("MMS linear solve did not converge");

    double sum_sq = 0.0;
    double sum_volume = 0.0;
    double linf = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        const double y = geometry.cell_centres[c].y;
        const double exact = y * y;
        const double error = result.solution(c) - exact;
        sum_sq += geometry.cell_volumes[c] * error * error;
        sum_volume += geometry.cell_volumes[c];
        linf = std::max(linf, std::abs(error));
    }
    return {std::sqrt(sum_sq / sum_volume), linf};
}


Mesh make_cartesian_cube(std::size_t n)
{
    if (n == 0) throw std::invalid_argument("make_cartesian_cube: n must be positive");

    const double h = 1.0 / static_cast<double>(n);
    const std::size_t nc = n * n * n;
    Mesh m;
    m.points().resize(8 * nc);

    const auto cell_id = [n](std::size_t i, std::size_t j, std::size_t k) {
        return (k * n + j) * n + i;
    };

    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t c = cell_id(i, j, k);
                const std::size_t b = 8 * c;
                const double x0 = i * h, x1 = (i + 1) * h;
                const double y0 = j * h, y1 = (j + 1) * h;
                const double z0 = k * h, z1 = (k + 1) * h;
                const double p[8][3] = {
                    {x0,y0,z0}, {x1,y0,z0}, {x1,y1,z0}, {x0,y1,z0},
                    {x0,y0,z1}, {x1,y0,z1}, {x1,y1,z1}, {x0,y1,z1}
                };
                for (std::size_t q = 0; q < 8; ++q)
                    m.points().set(b + q, p[q][0], p[q][1], p[q][2]);
            }
        }
    }

    std::vector<std::size_t> cell_faces(nc);
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
                const std::size_t b = 8 * c;
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
        const std::size_t owner_i = key.axis == 0 ? key.i - 1 : key.i;
        const std::size_t owner_j = key.axis == 1 ? key.j - 1 : key.j;
        const std::size_t owner_k = key.axis == 2 ? key.k - 1 : key.k;
        const bool lower_boundary = (key.axis == 0 ? key.i == 0 :
                                     key.axis == 1 ? key.j == 0 : key.k == 0);
        const bool upper_boundary = (key.axis == 0 ? key.i == n :
                                     key.axis == 1 ? key.j == n : key.k == n);
        const std::size_t owner = lower_boundary
            ? cell_id(0, key.axis == 0 ? key.j : owner_j, key.axis == 2 ? key.k : owner_k)
            : cell_id(owner_i, owner_j, owner_k);
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

    for (std::size_t c = 0; c < nc; ++c)
        m.cells().push_cell(faces_per_cell[c]);

    return m;
}

struct ThreeDError {
    double l2;
    double linf;
};

ThreeDError solve_3d_mms(std::size_t n)
{
    const Mesh mesh = make_cartesian_cube(n);
    const auto geometry = make_geometry_cache(mesh);
    constexpr double gamma = 1.0;
    constexpr double pi = 3.14159265358979323846;

    // phi = sin(pi*x) sin(pi*y) sin(pi*z)
    // -div(gamma grad(phi)) = 3*pi^2*gamma*phi.
    std::vector<double> source(mesh.n_cells(), 0.0);
    PoissonBoundaryCondition bc =
        PoissonBoundaryCondition::dirichlet(mesh.n_faces());

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const std::size_t owner = mesh.ownership().owner(f);
        const auto& fc = geometry.face_centres[f];
        const double phi_b = std::sin(pi * fc.x) *
                             std::sin(pi * fc.y) *
                             std::sin(pi * fc.z);
        bc.face_values[f] = phi_b;
        (void)owner;
    }

    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        const auto& cc = geometry.cell_centres[c];
        const double phi = std::sin(pi * cc.x) *
                           std::sin(pi * cc.y) *
                           std::sin(pi * cc.z);
        source[c] = 3.0 * pi * pi * gamma * phi;
    }

    ScalarDiffusionConfig cfg;
    cfg.diffusivity = gamma;
    cfg.max_iterations = 10000;
    cfg.tolerance = 1e-12;
    const auto result = solve_poisson_dirichlet(mesh, bc, source, cfg);
    if (result.linear_result.status != SolverStatus::CONVERGED)
        throw std::runtime_error("3D Poisson MMS linear solve did not converge");

    double sum_sq = 0.0;
    double sum_volume = 0.0;
    double linf = 0.0;
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        const auto& cc = geometry.cell_centres[c];
        const double exact = std::sin(pi * cc.x) *
                             std::sin(pi * cc.y) *
                             std::sin(pi * cc.z);
        const double error = result.solution(c) - exact;
        sum_sq += geometry.cell_volumes[c] * error * error;
        sum_volume += geometry.cell_volumes[c];
        linf = std::max(linf, std::abs(error));
    }
    return {std::sqrt(sum_sq / sum_volume), linf};
}

} // namespace

int main()
{
    std::vector<double> l2;
    const std::size_t levels[] = {8, 16, 32, 64};

    for (const auto n : levels) {
        const auto e = solve_mms(n);
        l2.push_back(e.l2);
        std::cout << "MMS phi=y^2 N=" << n
                  << " L2=" << e.l2
                  << " Linf=" << e.linf << "\n";
        EXPECT_TRUE(std::isfinite(e.l2));
        EXPECT_TRUE(std::isfinite(e.linf));
        EXPECT_TRUE(e.l2 > 0.0);
    }

    for (std::size_t i = 1; i < l2.size(); ++i) {
        const double order = std::log(l2[i-1] / l2[i]) / std::log(2.0);
        std::cout << "MMS observed order " << levels[i-1] << "->" << levels[i]
                  << ": " << order << "\n";
        EXPECT_TRUE(order >= 1.90);
    }


    std::vector<double> l2_3d;
    const std::size_t levels_3d[] = {4, 8, 16, 32};

    for (const auto n : levels_3d) {
        const auto e = solve_3d_mms(n);
        l2_3d.push_back(e.l2);
        std::cout << "Poisson 3D phi=sin(pi*x)sin(pi*y)sin(pi*z) N=" << n
                  << " L2=" << e.l2
                  << " Linf=" << e.linf << "\n";
        EXPECT_TRUE(std::isfinite(e.l2));
        EXPECT_TRUE(std::isfinite(e.linf));
        EXPECT_TRUE(e.l2 > 0.0);
    }

    for (std::size_t i = 1; i < l2_3d.size(); ++i) {
        const double order =
            std::log(l2_3d[i - 1] / l2_3d[i]) / std::log(2.0);
        std::cout << "Poisson 3D observed order "
                  << levels_3d[i - 1] << "->" << levels_3d[i]
                  << ": " << order << "\n";
        EXPECT_TRUE(order >= 1.80);
    }

    return run_all();
}
