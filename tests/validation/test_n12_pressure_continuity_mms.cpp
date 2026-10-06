#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
constexpr double PI = 3.1415926535897932384626433832795;

Mesh make_channel(std::size_t n)
{
    Mesh m;
    m.points().resize(8 * n);
    const double dy = 1.0 / static_cast<double>(n);

    for (std::size_t j = 0; j < n; ++j) {
        const double y0 = j * dy;
        const double y1 = (j + 1) * dy;
        const std::size_t b = 8 * j;
        const double p[8][3] = {
            {0.0,y0,0.0}, {1.0,y0,0.0}, {1.0,y1,0.0}, {0.0,y1,0.0},
            {0.0,y0,1.0}, {1.0,y0,1.0}, {1.0,y1,1.0}, {0.0,y1,1.0}
        };
        for (std::size_t q = 0; q < 8; ++q)
            m.points().set(b + q, p[q][0], p[q][1], p[q][2]);
    }

    std::vector<std::size_t> bottom, top, x0, x1, z0, z1, internal;
    auto add_face = [&](std::initializer_list<std::size_t> v) {
        const auto id = m.faces().n_faces();
        m.faces().push_face(std::vector<FaceIndex>(v.begin(), v.end()));
        return id;
    };

    bottom.push_back(add_face({0,1,5,4}));
    top.push_back(add_face({8*(n-1)+3,8*(n-1)+7,8*(n-1)+6,8*(n-1)+2}));
    for (std::size_t j = 0; j < n; ++j) {
        const auto b = 8 * j;
        x0.push_back(add_face({b,b+4,b+7,b+3}));
        x1.push_back(add_face({b+1,b+2,b+6,b+5}));
        z0.push_back(add_face({b,b+3,b+2,b+1}));
        z1.push_back(add_face({b+4,b+5,b+6,b+7}));
    }
    for (std::size_t j = 0; j + 1 < n; ++j) {
        const auto b = 8 * j;
        internal.push_back(add_face({b+3,b+7,b+6,b+2}));
    }

    m.ownership().resize(m.n_faces());
    m.ownership().set_owner(bottom[0], 0);
    m.ownership().set_neighbour(bottom[0], FaceOwnership::BOUNDARY);
    m.ownership().set_owner(top[0], n - 1);
    m.ownership().set_neighbour(top[0], FaceOwnership::BOUNDARY);
    for (std::size_t j = 0; j < n; ++j) {
        for (const auto f : {x0[j], x1[j], z0[j], z1[j]}) {
            m.ownership().set_owner(f, j);
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        }
    }
    for (std::size_t j = 0; j + 1 < n; ++j) {
        m.ownership().set_owner(internal[j], j);
        m.ownership().set_neighbour(internal[j], static_cast<std::int64_t>(j + 1));
    }
    for (std::size_t j = 0; j < n; ++j) {
        m.cells().push_cell({
            j == 0 ? bottom[0] : internal[j-1],
            j + 1 == n ? top[0] : internal[j],
            x0[j], x1[j], z0[j], z1[j]
        });
    }

    auto patch = [&](const char* name, const std::vector<std::size_t>& faces) {
        Patch p;
        p.name = name;
        p.type = PatchType::WALL;
        p.face_ids = faces;
        m.boundary().add_patch(p);
    };
    patch("bottom", bottom);
    patch("top", top);
    patch("x0", x0);
    patch("x1", x1);
    patch("z0", z0);
    patch("z1", z1);
    return m;
}

double ux(double y)
{
    return std::sin(PI * y) * std::sin(PI * y);
}

struct Norms {
    double l1;
    double l2;
    double linf;
};

Norms continuity_error(const Mesh& mesh, const cfdx::physics::FvGeometry& geometry)
{
    Field<double, Location::FACE> phi(mesh.n_faces(), "phi_mms", "kg/s", 1);

    // Exact face mass flux for U=(sin²(pi y),0,0).  The only non-zero
    // boundary fluxes are x-normal; the y/z faces carry exactly zero flux.
    for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
        const auto& fc = geometry.face_centres[face];
        phi(face) = ux(fc.y) * geometry.face_area_vectors[face].x;
    }

    double l1 = 0.0;
    double l2sq = 0.0;
    double linf = 0.0;
    double volume = 0.0;

    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        double net = 0.0;
        const auto off = mesh.cells().offsets_data()[c];
        const auto end = mesh.cells().offsets_data()[c + 1];
        for (auto k = off; k < end; ++k) {
            const auto face = mesh.cells().faces_data()[k];
            const bool owner = mesh.ownership().owner(face) == c;
            net += owner ? phi(face) : -phi(face);
        }
        const double r = net / geometry.cell_volumes[c];
        l1 += std::abs(r) * geometry.cell_volumes[c];
        l2sq += r * r * geometry.cell_volumes[c];
        linf = std::max(linf, std::abs(r));
        volume += geometry.cell_volumes[c];
    }

    return {l1 / volume, std::sqrt(l2sq / volume), linf};
}
} // namespace

int main()
{
    try {
        for (const std::size_t n : {16u, 32u, 64u}) {
            const Mesh mesh = make_channel(n);
            const auto geometry = cfdx::physics::build_fv_geometry(mesh);
            const auto e = continuity_error(mesh, geometry);

            std::cout << std::setprecision(12)
                      << "N12_CONTINUITY_MMS n=" << n
                      << " L1=" << e.l1
                      << " L2=" << e.l2
                      << " Linf=" << e.linf << "\n";

            if (!(e.l1 < 1e-13 && e.l2 < 1e-13 && e.linf < 1e-13))
                throw std::runtime_error(
                    "discrete continuity is not exactly satisfied by the manufactured divergence-free field");
        }

        std::cout << "N12_CONTINUITY_MMS: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "N12_CONTINUITY_MMS: FAIL: " << e.what() << "\n";
        return 1;
    }
}
