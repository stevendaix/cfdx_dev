#pragma once

#include "cfdx/physics/steady_incompressible_solver.h"
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

namespace cfdx::validation {

inline std::shared_ptr<const cfdx::core::ValueProvider> constant_provider(double value)
{
    return std::make_shared<cfdx::core::ConstantValueProvider>(value);
}

inline void add_velocity_dirichlet(
    cfdx::physics::BoundaryConstraintMap& bc,
    const std::string& patch,
    const cfdx::core::Vec3& value)
{
    bc[patch].emplace_back("U.x", cfdx::core::Dirichlet{constant_provider(value.x)});
    bc[patch].emplace_back("U.y", cfdx::core::Dirichlet{constant_provider(value.y)});
    bc[patch].emplace_back("U.z", cfdx::core::Dirichlet{constant_provider(value.z)});
}

inline void add_velocity_flux_dependent(
    cfdx::physics::BoundaryConstraintMap& bc,
    const std::string& patch,
    const cfdx::core::Vec3& inflow,
    const cfdx::core::Vec3& outflow_gradient)
{
    bc[patch].emplace_back(
        "U.x",
        cfdx::core::FluxDependent{
            constant_provider(inflow.x), constant_provider(outflow_gradient.x)});
    bc[patch].emplace_back(
        "U.y",
        cfdx::core::FluxDependent{
            constant_provider(inflow.y), constant_provider(outflow_gradient.y)});
    bc[patch].emplace_back(
        "U.z",
        cfdx::core::FluxDependent{
            constant_provider(inflow.z), constant_provider(outflow_gradient.z)});
}

inline void add_pressure_dirichlet(
    cfdx::physics::BoundaryConstraintMap& bc,
    const std::string& patch,
    double value)
{
    bc[patch].emplace_back("p", cfdx::core::Dirichlet{constant_provider(value)});
}

inline void add_pressure_neumann(
    cfdx::physics::BoundaryConstraintMap& bc,
    const std::string& patch,
    double gradient)
{
    bc[patch].emplace_back("p", cfdx::core::Neumann{constant_provider(gradient)});
}

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_validation_face_flux(
    const cfdx::core::Mesh& mesh,
    const cfdx::physics::FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U)
{
    if (U.dimension() != 3 || U.size() != mesh.n_cells())
        throw std::invalid_argument("validation BC preflight: invalid velocity field");

    cfdx::core::Field<double, cfdx::core::Location::FACE>
        flux(mesh.n_faces(), "validation_phi", "kg/s", 1);

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const auto& own = mesh.ownership();
        const std::size_t o = own.owner(f);
        cfdx::core::Vec3 uf;
        U.get(o, uf.x, uf.y, uf.z);

        const int n = own.neighbour(f);
        if (n >= 0) {
            cfdx::core::Vec3 un;
            U.get(static_cast<std::size_t>(n), un.x, un.y, un.z);
            uf = (uf + un) * 0.5;
        } else {
            const std::size_t patch = geometry.face_patch[f];
            if (patch < mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch).type == cfdx::core::PatchType::EMPTY) {
                flux(f) = 0.0;
                continue;
            }
        }
        flux(f) = uf.dot(geometry.face_area_vectors[f]);
    }
    return flux;
}

inline void exercise_new_velocity_bc_contract(
    const cfdx::core::Mesh& mesh,
    const cfdx::physics::FvGeometry& geometry,
    const cfdx::physics::BoundaryConstraintMap& constraints,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const std::string& label)
{
    if (mass_flux.size() != mesh.n_faces() || mass_flux.dimension() != 1)
        throw std::invalid_argument(label + ": invalid validation face flux");

    cfdx::core::Field<double, cfdx::core::Location::CELL>
        grad_p(mesh.n_cells(), "bc_preflight_gradp", "Pa/m", 3);
    cfdx::core::Field<double, cfdx::core::Location::CELL>
        body(mesh.n_cells(), "bc_preflight_body", "N/m3", 1);
    grad_p.fill(0.0);
    body.fill(0.0);

    for (std::size_t component = 0; component < 3; ++component) {
        const std::string field =
            "U." + std::string(1, static_cast<char>('x' + component));
        const auto resolved = cfdx::physics::resolve_scalar_boundary_constraints(
            mesh, geometry.face_centres, constraints, field, 0.0, &mass_flux);

        for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
            const int n = mesh.ownership().neighbour(face);
            if (n >= 0) continue;
            const std::size_t patch = geometry.face_patch[face];
            if (patch >= mesh.boundary().n_patches()) continue;
            if (mesh.boundary().patch(patch).type == cfdx::core::PatchType::EMPTY) continue;
            if (!resolved.has(face))
                throw std::runtime_error(
                    label + ": new BC resolver left boundary face unspecified");
        }

        const auto equation = cfdx::physics::assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body, 1.0,
            constraints, component, false, 0.0,
            cfdx::physics::ConvectionScheme::UPWIND, nullptr);
        if (equation.matrix.n_rows() != mesh.n_cells())
            throw std::runtime_error(label + ": new BC FVM assembly returned invalid matrix");
    }
}

} // namespace cfdx::validation
