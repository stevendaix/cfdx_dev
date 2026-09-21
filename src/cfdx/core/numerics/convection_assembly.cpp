#include "convection_assembly.h"
#include "cfdx/core/geometry/face_geometry.h"
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

    for (std::size_t c=0; c<mesh.n_cells(); ++c) {
        const auto off=mesh.cells().offsets_data()[c];
        const auto n=mesh.cells().offsets_data()[c+1]-off;
        Vec3 centre;
        double weight=0.0;
        for (std::size_t k=0;k<n;++k) {
            const auto f=mesh.cells().faces_data()[off+k];
            const double a=g.face_area_vectors[f].mag();
            centre=centre+g.face_centres[f]*a;
            weight+=a;
        }
        if (!(weight>0.0)) throw std::runtime_error("convection assembly: invalid cell geometry");
        g.cell_centres[c]=centre*(1.0/weight);
    }
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
