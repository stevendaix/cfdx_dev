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

// Structured triangular-prism extrusion of a circular disk. The circular
// wall is represented by straight chords; refinement therefore measures both
// axial FV error and the geometric approximation error.
Mesh make_pipe(std::size_t nr, std::size_t nz, double R, double L)
{
    constexpr std::size_t NS = 16;
    Mesh m;
    const std::size_t pts_per_layer = 1 + NS;
    m.points().resize(pts_per_layer * (nz + 1));
    const double dz = L / static_cast<double>(nz);

    for (std::size_t k = 0; k <= nz; ++k) {
        const double z = dz * static_cast<double>(k);
        const std::size_t b = k * pts_per_layer;
        m.points().set(b, 0.0, 0.0, z);
        for (std::size_t s = 0; s < NS; ++s) {
            const double a = 2.0 * M_PI * static_cast<double>(s) / NS;
            m.points().set(b + 1 + s, R * std::cos(a), R * std::sin(a), z);
        }
    }

    std::vector<std::size_t> inlet, outlet, wall, internal;
    const std::size_t ntri = NS;
    for (std::size_t k = 0; k < nz; ++k) {
        const std::size_t b0 = k * pts_per_layer;
        const std::size_t b1 = (k + 1) * pts_per_layer;
        for (std::size_t s = 0; s < NS; ++s) {
            const std::size_t sn = (s + 1) % NS;
            // Triangle faces are oriented with outward normals at the ends.
            const std::size_t fi = m.faces().n_faces();
            m.faces().push_face({b0, b1 + 1 + sn, b1 + 1 + s});
            if (k == 0) inlet.push_back(fi);
            if (k + 1 == nz) {
                const std::size_t fo = m.faces().n_faces();
                m.faces().push_face({b1, b1 + 1 + s, b1 + 1 + sn});
                outlet.push_back(fo);
            }

            const std::size_t fw = m.faces().n_faces();
            m.faces().push_face({b0 + 1 + s, b0 + 1 + sn, b1 + 1 + sn, b1 + 1 + s});
            wall.push_back(fw);

            if (k + 1 < nz) {
                const std::size_t f = m.faces().n_faces();
                m.faces().push_face({b1, b1 + 1 + s, b1 + 1 + sn});
                internal.push_back(f);
            }
        }
    }

    // Rebuild a clean face set: the construction above intentionally creates
    // temporary duplicates. Use a compact connectivity map by matching the
    // canonical vertex sets of faces.
    struct Key { std::vector<std::size_t> v; bool operator<(const Key& o) const { return v < o.v; } };
    std::map<Key,std::size_t> face_map;
    Mesh compact;
    compact.points() = m.points();
    (void)compact;
    throw std::runtime_error("VMFL005 mesh generator not yet compacted");
}

} // namespace

int main() { return 0; }
