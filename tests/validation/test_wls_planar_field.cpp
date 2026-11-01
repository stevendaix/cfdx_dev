#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/mesh/mesh.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using namespace cfdx::core;
namespace {
Mesh make_channel_mesh(std::size_t nx, std::size_t ny)
{
    if (nx < 2 || ny < 4) throw std::invalid_argument("channel mesh too small");

    Mesh mesh;
    const std::size_t plane = (nx + 1) * (ny + 1);
    mesh.points().resize(2 * plane);
    const auto id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };

    for (std::size_t j = 0; j <= ny; ++j) {
        for (std::size_t i = 0; i <= nx; ++i) {
            const double x = static_cast<double>(i) / static_cast<double>(nx);
            const double y = static_cast<double>(j) / static_cast<double>(ny);
            mesh.points().set(id(i, j, 0), x, y, 0.0);
            mesh.points().set(id(i, j, 1), x, y, 1.0);
        }
    }

    std::map<std::vector<std::size_t>, std::size_t> face_map;
    std::vector<std::vector<std::size_t>> cell_faces(nx * ny);

    auto add_face = [&](std::initializer_list<std::size_t> vertices,
                        std::size_t cell) {
        std::vector<std::size_t> key(vertices);
        std::sort(key.begin(), key.end());
        const auto found = face_map.find(key);
        if (found != face_map.end()) {
            mesh.ownership().set_neighbour(found->second, static_cast<int>(cell));
            return found->second;
        }

        const std::size_t f = mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        mesh.ownership().set_owner(f, cell);
        mesh.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        face_map.emplace(std::move(key), f);
        return f;
    };

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t c = j * nx + i;
            const auto a = id(i, j, 0);
            const auto b = id(i + 1, j, 0);
            const auto c0 = id(i + 1, j + 1, 0);
            const auto d = id(i, j + 1, 0);
            const auto e = id(i, j, 1);
            const auto f = id(i + 1, j, 1);
            const auto g = id(i + 1, j + 1, 1);
            const auto h = id(i, j + 1, 1);
            cell_faces[c] = {
                add_face({a, d, c0, b}, c),
                add_face({e, f, g, h}, c),
                add_face({a, b, f, e}, c),
                add_face({d, h, g, c0}, c),
                add_face({a, e, h, d}, c),
                add_face({b, c0, g, f}, c)};
        }
    }

    for (const auto& faces : cell_faces) mesh.cells().push_cell(faces);

    Patch inlet{"inlet", PatchType::INLET, {}};
    Patch outlet{"outlet", PatchType::OUTLET, {}};
    Patch bottom{"bottom", PatchType::WALL, {}};
    Patch top{"top", PatchType::WALL, {}};
    Patch front{"front", PatchType::EMPTY, {}};
    Patch back{"back", PatchType::EMPTY, {}};

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) >= 0) continue;

        const auto& vertices = mesh.faces().vertices();
        const auto begin = vertices.begin() +
            static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto end = begin +
            static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));

        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        for (auto it = begin; it != end; ++it) {
            x += mesh.points().x(*it);
            y += mesh.points().y(*it);
            z += mesh.points().z(*it);
        }
        const double n = static_cast<double>(mesh.faces().face_size(f));
        x /= n; y /= n; z /= n;

        constexpr double tol = 1e-12;
        if (std::abs(x) < tol) inlet.face_ids.push_back(f);
        else if (std::abs(x - 1.0) < tol) outlet.face_ids.push_back(f);
        else if (std::abs(y) < tol) bottom.face_ids.push_back(f);
        else if (std::abs(y - 1.0) < tol) top.face_ids.push_back(f);
        else if (std::abs(z) < tol) front.face_ids.push_back(f);
        else if (std::abs(z - 1.0) < tol) back.face_ids.push_back(f);
        else throw std::runtime_error("channel boundary face is unclassified");
    }

    mesh.boundary().add_patch(inlet);
    mesh.boundary().add_patch(outlet);
    mesh.boundary().add_patch(bottom);
    mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(front);
    mesh.boundary().add_patch(back);
    return mesh;
}


}
int main() {
    try {
        const auto mesh = make_channel_mesh(4,4);
        const auto geom = make_geometry_cache(mesh);
        ScalarCellField field(mesh.n_cells(), "affine", "1", 1);
        for (std::size_t c=0;c<mesh.n_cells();++c) {
            const auto pos=geom.cell_centres[c];
            field(c)=3.0+2.0*pos.x-3.0*pos.y;
        }
        const auto gradient=compute_gradient_weighted_least_squares(
            field,mesh,GradientWeighting::INVERSE_DISTANCE_SQUARED,
            1.e5,BoundaryGradientPolicy::EXCLUDE_BOUNDARY,nullptr,2,true);
        for (std::size_t c=0;c<mesh.n_cells();++c) {
            if (std::abs(gradient(c,0)-2.0)>1.e-10 ||
                std::abs(gradient(c,1)+3.0)>1.e-10 ||
                std::abs(gradient(c,2))>1.e-10)
                throw std::runtime_error("field WLS planar affine consistency failed");
        }
        bool rejected=false;
        try {
            (void)compute_gradient_weighted_least_squares(
                field,mesh,GradientWeighting::INVERSE_DISTANCE_SQUARED,
                1.e5,BoundaryGradientPolicy::EXCLUDE_BOUNDARY,nullptr,3,true);
        } catch (const std::runtime_error&) { rejected=true; }
        if (!rejected) throw std::runtime_error("3-D field mode accepted planar stencil");
        std::cout << "WLS_PLANAR_FIELD: PASS cells="<<mesh.n_cells()<<'\n';
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "WLS_PLANAR_FIELD: FAIL "<<e.what()<<'\n';
        return 1;
    }
}
