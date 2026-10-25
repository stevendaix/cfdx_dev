#pragma once
#include "cfdx/physics/finite_volume_transport.h"
using namespace cfdx::core;
static Mesh make_1d_chain(std::size_t n)
{
    if (n == 0) throw std::invalid_argument("make_1d_chain: n must be positive");
    Mesh m;
    m.points().resize(8 * n);
    const double dx = 1.0 / static_cast<double>(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double x0 = dx * static_cast<double>(i);
        const double x1 = dx * static_cast<double>(i + 1);
        const std::size_t b = 8 * i;
        const double p[8][3] = {
            {x0,0,0}, {x1,0,0}, {x1,1,0}, {x0,1,0},
            {x0,0,1}, {x1,0,1}, {x1,1,1}, {x0,1,1}
        };
        for (std::size_t q = 0; q < 8; ++q)
            m.points().set(b + q, p[q][0], p[q][1], p[q][2]);
    }
    std::vector<std::size_t> left, right, y0, y1, z0, z1, internal;
    auto add_face = [&](std::initializer_list<std::size_t> v) {
        const std::size_t id = m.faces().n_faces();
        m.faces().push_face(std::vector<FaceIndex>(v.begin(), v.end()));
        return id;
    };
    left.push_back(add_face({0,4,7,3}));
    right.push_back(add_face({8*(n-1)+1,8*(n-1)+2,8*(n-1)+6,8*(n-1)+5}));
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t b = 8 * i;
        y0.push_back(add_face({b,b+1,b+5,b+4}));
        y1.push_back(add_face({b+3,b+7,b+6,b+2}));
        z0.push_back(add_face({b,b+3,b+2,b+1}));
        z1.push_back(add_face({b+4,b+5,b+6,b+7}));
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        const std::size_t b = 8 * i;
        internal.push_back(add_face({b+1,b+2,b+6,b+5}));
    }
    m.ownership().resize(m.n_faces());
    m.ownership().set_owner(left[0], 0);
    m.ownership().set_neighbour(left[0], FaceOwnership::BOUNDARY);
    m.ownership().set_owner(right[0], n - 1);
    m.ownership().set_neighbour(right[0], FaceOwnership::BOUNDARY);
    for (std::size_t i = 0; i < n; ++i) {
        for (const auto f : {y0[i], y1[i], z0[i], z1[i]}) {
            m.ownership().set_owner(f, i);
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        }
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        m.ownership().set_owner(internal[i], i);
        m.ownership().set_neighbour(internal[i], static_cast<std::int64_t>(i + 1));
    }
    for (std::size_t i = 0; i < n; ++i)
        m.cells().push_cell({i == 0 ? left[0] : internal[i-1],
                             i + 1 == n ? right[0] : internal[i],
                             y0[i], y1[i], z0[i], z1[i]});

    Patch inlet; inlet.name = "inlet"; inlet.type = PatchType::WALL;
    inlet.face_ids = left;
    Patch outlet; outlet.name = "outlet"; outlet.type = PatchType::WALL;
    outlet.face_ids = right;
    Patch wall; wall.name = "wall"; wall.type = PatchType::WALL;
    for (const auto f : y0) wall.face_ids.push_back(f);
    for (const auto f : y1) wall.face_ids.push_back(f);
    for (const auto f : z0) wall.face_ids.push_back(f);
    for (const auto f : z1) wall.face_ids.push_back(f);
    m.boundary().add_patch(inlet);
    m.boundary().add_patch(outlet);
    m.boundary().add_patch(wall);
    return m;
}
