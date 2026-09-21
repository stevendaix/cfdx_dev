#include "convection_assembly.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::core::numerics;

namespace {

struct GeometryData {
    std::vector<Vec3> face_centres;
    std::vector<Vec3> face_area_vectors;
    std::vector<Vec3> cell_centres;
};

GeometryData build_geometry(const Mesh& mesh)
{
    GeometryData g;
    g.face_centres.resize(mesh.n_faces());
    g.face_area_vectors.resize(mesh.n_faces());
    g.cell_centres.resize(mesh.n_cells());

    for (std::size_t f=0; f<mesh.n_faces(); ++f) {
        const auto off=mesh.faces().offsets_data()[f];
        const auto n=mesh.faces().offsets_data()[f+1]-off;
        const auto fg=compute_face_geometry(
            mesh.points().x_data(),mesh.points().y_data(),mesh.points().z_data(),
            mesh.faces().vertices_data(),off,n);
        g.face_centres[f]=fg.centre;
        g.face_area_vectors[f]=fg.Sf;
    }

    compute_area_weighted_cell_centres(
        mesh, g.face_centres.data(), g.face_area_vectors.data(), g.cell_centres.data());
    orient_mesh_face_vectors(mesh, g.face_centres, g.cell_centres, g.face_area_vectors);
    return g;
}

double scheme_weight(const std::string& scheme,double F)
{
    if(scheme=="upwind") return F>=0.0 ? 1.0 : 0.0;
    if(scheme=="linear") return 0.5;
    if(scheme=="limited") {
        // Bounded blend between upwind and central interpolation.
        return F>=0.0 ? 0.75 : 0.25;
    }
    throw std::invalid_argument("assembleConvectionCSR: unsupported scheme: "+scheme);
}
