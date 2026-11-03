#pragma once

#include "cfdx/physics/forces.h"
#include "cfdx/physics/wall_quadratic_stencil.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

// Convert a mesh WALL patch into the solver-agnostic force integration
// contract. The wall pressure and velocity are explicit face-indexed inputs;
// this adapter never guesses boundary values or falls back to cell gradients.
inline std::vector<forces::SurfaceSample> reconstruct_quadratic_wall_surface(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    std::size_t wall_patch,
    const std::vector<double>& cell_pressure,
    const std::vector<cfdx::core::Vec3>& cell_velocity,
    const std::vector<double>& wall_pressure_by_face,
    const std::vector<cfdx::core::Vec3>& wall_velocity_by_face,
    double dynamic_viscosity,
    std::size_t max_rings = 3)
{
    using cfdx::core::PatchType;
    using cfdx::core::Vec3;
    if (wall_patch >= mesh.boundary().n_patches() ||
        mesh.boundary().patch(wall_patch).type != PatchType::WALL ||
        geometry.face_centres.size() != mesh.n_faces() ||
        geometry.face_area_vectors.size() != mesh.n_faces() ||
        wall_pressure_by_face.size() != mesh.n_faces() ||
        wall_velocity_by_face.size() != mesh.n_faces() ||
        !(dynamic_viscosity >= 0.0) || !std::isfinite(dynamic_viscosity))
        throw std::invalid_argument("quadratic wall forces: invalid patch, fields or viscosity");

    const auto& patch = mesh.boundary().patch(wall_patch);
    std::vector<forces::SurfaceSample> samples;
    samples.reserve(patch.face_ids.size());
    for (const auto face : patch.face_ids) {
        if (face >= mesh.n_faces())
            throw std::invalid_argument("quadratic wall forces: patch face is out of range");
        const auto reconstructed = reconstruct_wall_face_quadratic(
            build_wall_quadratic_stencil(
                mesh, geometry, face, cell_pressure, cell_velocity,
                wall_velocity_by_face[face], wall_pressure_by_face[face], max_rings),
            geometry.face_centres[face], geometry.face_area_vectors[face],
            dynamic_viscosity);

        const Vec3 fluid_outward =
            geometry.face_area_vectors[face].normalized();
        const Vec3 body_outward = fluid_outward * -1.0;
        const Vec3 viscous_on_body =
            reconstructed.force_per_area - fluid_outward * reconstructed.pressure;
        const double area = geometry.face_area_vectors[face].mag();
        forces::SurfaceSample sample;
        sample.position = {geometry.face_centres[face].x,
                           geometry.face_centres[face].y,
                           geometry.face_centres[face].z};
        sample.normal = {body_outward.x, body_outward.y, body_outward.z};
        sample.area = area;
        sample.pressure = reconstructed.pressure;
        sample.viscous_traction = {viscous_on_body.x, viscous_on_body.y,
                                   viscous_on_body.z};
        samples.push_back(sample);
    }
    return samples;
}

} // namespace cfdx::physics
