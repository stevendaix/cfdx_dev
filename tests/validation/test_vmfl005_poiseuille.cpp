#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {

Mesh make_pipe(std::size_t ns, std::size_t nz, double R, double L)
{
    Mesh m;
    const std::size_t p_layer = 1 + ns;
    m.points().resize(p_layer * (nz + 1));
    const double dz = L / static_cast<double>(nz);
    constexpr double pi = 3.1415926535897932384626433832795;

    for (std::size_t k = 0; k <= nz; ++k) {
        const std::size_t b = k * p_layer;
        const double z = dz * static_cast<double>(k);
        m.points().set(b, 0.0, 0.0, z);
        for (std::size_t s = 0; s < ns; ++s) {
            const double a = 2.0 * pi * static_cast<double>(s) / static_cast<double>(ns);
            m.points().set(b + 1 + s, R * std::cos(a), R * std::sin(a), z);
        }
    }

    std::vector<std::vector<std::size_t>> tri(nz + 1, std::vector<std::size_t>(ns));
    std::vector<std::vector<std::size_t>> wall(nz, std::vector<std::size_t>(ns));

    for (std::size_t k = 0; k <= nz; ++k) {
        const std::size_t b = k * p_layer;
        for (std::size_t s = 0; s < ns; ++s) {
            const std::size_t sn = (s + 1) % ns;
            tri[k][s] = m.faces().n_faces();
            if (k == 0)
                m.faces().push_face({b, b + 1 + sn, b + 1 + s});
            else
                m.faces().push_face({b, b + 1 + s, b + 1 + sn});
        }
    }

    for (std::size_t k = 0; k < nz; ++k) {
        const std::size_t b0 = k * p_layer;
        const std::size_t b1 = (k + 1) * p_layer;
        for (std::size_t s = 0; s < ns; ++s) {
            const std::size_t sn = (s + 1) % ns;
            wall[k][s] = m.faces().n_faces();
            m.faces().push_face({b0 + 1 + s, b0 + 1 + sn,
                                 b1 + 1 + sn, b1 + 1 + s});
        }
    }

    m.ownership().resize(m.n_faces());
    for (std::size_t k = 0; k <= nz; ++k) {
        for (std::size_t s = 0; s < ns; ++s) {
            const auto f = tri[k][s];
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
            m.ownership().set_owner(f, k == 0 ? 0 : k - 1);
            if (k > 0 && k < nz)
                m.ownership().set_neighbour(f, static_cast<int>(k));
        }
    }
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t s = 0; s < ns; ++s) {
            m.ownership().set_owner(wall[k][s], k);
            m.ownership().set_neighbour(wall[k][s], FaceOwnership::BOUNDARY);
        }
    }

    for (std::size_t k = 0; k < nz; ++k) {
        std::vector<std::size_t> faces = {
            tri[k][0], tri[k + 1][0], wall[k][0], wall[k][1], wall[k][ns - 1]
        };
        m.cells().push_cell(faces);
        // Each triangular prism is represented by the first sector's two end
        // triangles plus its two neighbouring wall sectors. The remaining
        // wall faces are attached below by adding all cells sector-wise.
        m.cells().pop_last_face_placeholder();
    }
    throw std::runtime_error("VMFL005: cell connectivity construction placeholder");
}

} // namespace

int main()
{
    std::cout << "VMFL005_VALIDATION: NOT_IMPLEMENTED\n";
    return 1;
}
