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
    std::vector<std::vector<std::size_t>> cell_faces(n);
    for (std::size_t i = 0; i < n; ++i) {
        cell_faces[i] = {
            i == 0 ? bottom[0] : internal[i-1],
            i + 1 == n ? top[0] : internal[i],
            x0[i], x1[i], z0[i], z1[i]
        };
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (const auto f : cell_faces[i]) {
            m.ownership().set_owner(f, i);
            const bool is_internal =
                std::find(internal.begin(), internal.end(), f) != internal.end();
            m.ownership().set_neighbour(
                f, is_internal ? static_cast<int>(i + 1) : FaceOwnership::BOUNDARY);
        }
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
    std::vector<double> source(n, -2.0 * gamma);
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

    return run_all();
}
