#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/io/restart/dat_restart.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "cfdx/physics/solver_control.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <functional>

namespace cfdx::physics {

struct VelocityBoundaryCondition {
    enum class Type { FIXED_VALUE, ZERO_GRADIENT };
    Type type = Type::ZERO_GRADIENT;
    cfdx::core::Vec3 value{0.0, 0.0, 0.0};
};

using VelocityBoundaryConditions = std::map<std::string, VelocityBoundaryCondition>;

enum class IncompressibleProbeField {
    PRESSURE,
    U_X,
    U_Y,
    U_Z,
    U_MAGNITUDE,
};

struct IncompressiblePointProbe {
    std::string name;
    cfdx::core::Vec3 location{0.0, 0.0, 0.0};
    IncompressibleProbeField field = IncompressibleProbeField::PRESSURE;
};

struct IncompressibleProbeSample {
    std::string name;
    std::size_t iteration = 0;
    double time = 0.0;
    double value = 0.0;
};

struct IncompressibleSolverControls {
    PressureVelocityAlgorithm algorithm = PressureVelocityAlgorithm::SIMPLE;
    CouplingControls coupling;
    ConvergenceCriteria convergence;
    std::size_t linear_max_iterations = 1000;
    double linear_tolerance = 1e-10;
    double density = 1.0;
    double kinematic_viscosity = 1.0e-3;
    double turbulent_viscosity = 0.0;
    cfdx::core::Vec3 body_force{0.0, 0.0, 0.0};
    std::size_t pressure_reference_cell = 0;
    double pressure_reference_value = 0.0;
    bool use_bounded_convection = true;
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND;
    std::vector<IncompressiblePointProbe> probes;
    std::function<void(const IncompressibleProbeSample&)> probe_callback;
    // Called after each completed nonlinear iteration with the authoritative solver state.
    // Return false to stop the solve after the current iteration.
    using IterationOutputCallback = std::function<bool(
        std::size_t,
        double,
        const cfdx::core::Mesh&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&)>;
    IterationOutputCallback iteration_output_callback;
};

struct IncompressibleIteration {
    std::size_t iteration = 0;
    double momentum_residual = std::numeric_limits<double>::infinity();
    double pressure_residual = std::numeric_limits<double>::infinity();
    double continuity_l1 = std::numeric_limits<double>::infinity();
    double continuity_linf = std::numeric_limits<double>::infinity();
    double continuity_normalized = std::numeric_limits<double>::infinity();
    double momentum_equation_residual = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_relative = std::numeric_limits<double>::infinity();
    double velocity_change_inf = std::numeric_limits<double>::infinity();
    double pressure_change_inf = std::numeric_limits<double>::infinity();
    std::size_t momentum_linear_iterations = 0;
    std::size_t pressure_linear_iterations = 0;
    // Conservative pressure-corrected flux continuity at the end of the iteration.
    double corrected_flux_continuity_linf = std::numeric_limits<double>::infinity();
    // Continuity obtained by reconstructing the face flux from the corrected cell velocity.
    double reconstructed_velocity_continuity_linf = std::numeric_limits<double>::infinity();
    // Maximum face-flux discrepancy between the conservative flux and reconstructed U flux.
    double flux_velocity_mismatch_linf = std::numeric_limits<double>::infinity();
    // Component-wise and location-aware diagnostics for independently
    // reassembled momentum equations. These are diagnostic signals only;
    // acceptance thresholds remain owned by the validation tests.
    std::array<double, 3> momentum_equation_residual_components{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    double momentum_equation_residual_internal = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_boundary = std::numeric_limits<double>::infinity();
    double pressure_gradient_linf = std::numeric_limits<double>::infinity();
    double pressure_gradient_l2 = std::numeric_limits<double>::infinity();
    std::size_t momentum_residual_cell = 0;
    double momentum_residual_no_pressure = std::numeric_limits<double>::infinity();
    double momentum_pressure_contribution = std::numeric_limits<double>::infinity();
    std::string momentum_residual_patch;
};

struct IncompressibleSolveResult {
    bool converged = false;
    std::size_t iterations = 0;
    std::vector<IncompressibleIteration> history;
};

inline void validate_incompressible_controls(
    const IncompressibleSolverControls& c, std::size_t n_cells)
{
    validate_coupling_controls(c.coupling);
    validate_convergence_criteria(c.convergence);
    if (!(c.density > 0.0) || !(c.kinematic_viscosity >= 0.0) ||
        !(c.turbulent_viscosity >= 0.0))
        throw std::invalid_argument("invalid incompressible material properties");
    if (c.pressure_reference_cell >= n_cells)
        throw std::invalid_argument("pressure reference cell out of range");
    if (c.linear_max_iterations == 0 || !(c.linear_tolerance > 0.0))
        throw std::invalid_argument("invalid linear solver controls");
}

inline cfdx::core::Field<double, cfdx::core::Location::CELL>
gauss_gradient_with_boundary(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& field,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarBoundaryConditions& bcs)
{
    using namespace cfdx::core;
    if (field.dimension() != 1 || field.size() != mesh.n_cells())
        throw std::invalid_argument("gauss_gradient_with_boundary: invalid field");

    Field<double, Location::CELL> grad(
        mesh.n_cells(), field.name() + "_grad_bc", field.metadata().unit + "/m", 3);
    const auto& own = mesh.ownership();

    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        Vec3 sum{0.0, 0.0, 0.0};
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        for (Offset k = 0; k < count; ++k) {
            const std::size_t f = mesh.cells().faces_data()[off + k];
            const bool owner = own.owner(f) == c;
            Vec3 Sf = owner ? geometry.face_area_vectors[f] : geometry.face_area_vectors[f] * -1.0;
            double vf = field(c);

            if (own.neighbour(f) >= 0) {
                const std::size_t n = static_cast<std::size_t>(own.neighbour(f));
                vf = 0.5 * (field(c) + field(n));
            } else {
                const std::size_t p = geometry.face_patch[f];
                if (p < mesh.boundary().n_patches() &&
                    mesh.boundary().patch(p).type == PatchType::EMPTY)
                    continue;
                if (p < mesh.boundary().n_patches()) {
                    const auto& name = mesh.boundary().patch(p).name;
                    const auto it = bcs.find(name);
                    if (it != bcs.end() && it->second.type == ScalarBoundaryType::FIXED_VALUE)
                        vf = it->second.value;
                }
            }
            sum += Sf * vf;
        }
        const double invV = 1.0 / geometry.cell_volumes[c];
        grad.component_data(0)[c] = sum.x * invV;
        grad.component_data(1)[c] = sum.y * invV;
        grad.component_data(2)[c] = sum.z * invV;
    }
    return grad;
}

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs = {})
{
    using namespace cfdx::core;
    (void)pressure_bcs;
    Field<double, Location::FACE> flux(mesh.n_faces(), "phi", "kg/s", 1);
    const auto& own = mesh.ownership();

    for (std::size_t c = 0; c < U.size(); ++c) {
        const double ux = U.component_data(0)[c];
        const double uy = U.component_data(1)[c];
        const double uz = U.component_data(2)[c];
        if (!std::isfinite(ux) || !std::isfinite(uy) || !std::isfinite(uz))
            throw std::runtime_error("make_mass_flux: non-finite cell velocity");
    }
    for (std::size_t f=0; f<mesh.n_faces(); ++f) {
        const std::size_t o=own.owner(f);
        Vec3 Uf;
        U.get(o,Uf.x,Uf.y,Uf.z);
        if (own.neighbour(f)>=0) {
            const std::size_t n=static_cast<std::size_t>(own.neighbour(f));
            Vec3 Un;
            U.get(n,Un.x,Un.y,Un.z);
            Uf=(Uf+Un)*0.5;
        } else {
            const std::size_t p=geometry.face_patch[f];
            if (p < mesh.boundary().n_patches() &&
                mesh.boundary().patch(p).type == PatchType::EMPTY) {
                flux(f) = 0.0;
                continue;
            }
            if(p<mesh.boundary().n_patches()) {
                const auto& name=mesh.boundary().patch(p).name;
                const auto it=bcs.find(name);
                if(it!=bcs.end() && it->second.type==VelocityBoundaryCondition::Type::FIXED_VALUE)
                    Uf=it->second.value;
            }
        }
        flux(f)=rho*Uf.dot(geometry.face_area_vectors[f]);
    }
    return flux;
}

inline double rhie_chow_pressure_flux_internal(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    std::size_t face,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& grad_p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho)
{
    using namespace cfdx::core;
    const auto nr = mesh.ownership().neighbour(face);
    if (nr < 0) return 0.0;
    const std::size_t o = mesh.ownership().owner(face);
    const std::size_t n = static_cast<std::size_t>(nr);
    const Vec3 Sf = geometry.face_area_vectors[face];
    const double area = Sf.mag();
    if (!(area > 0.0))
        throw std::runtime_error("rhie_chow_pressure_flux_internal: degenerate face");
    const Vec3 dvec = geometry.cell_centres[n] - geometry.cell_centres[o];
    const double d = dvec.mag();
    if (!(d > 0.0))
        throw std::runtime_error("rhie_chow_pressure_flux_internal: degenerate centre distance");
    const Vec3 e{dvec.x/d, dvec.y/d, dvec.z/d};
    // rAU is the inverse integrated momentum diagonal (1/A_P). The volume
    // factor is introduced here, at the velocity/pressure projection stage,
    // because the pressure gradient is a volumetric momentum source. Thus
    // cx/cy/cz are dAU = V*rAU, not a redefinition of rAU itself.
    const double cx = 0.5*(geometry.cell_volumes[o] * rAU[0][o] +
                           geometry.cell_volumes[n] * rAU[0][n]);
    const double cy = 0.5*(geometry.cell_volumes[o] * rAU[1][o] +
                           geometry.cell_volumes[n] * rAU[1][n]);
    const double cz = 0.5*(geometry.cell_volumes[o] * rAU[2][o] +
                           geometry.cell_volumes[n] * rAU[2][n]);
    const double rfn = cx*e.x*e.x + cy*e.y*e.y + cz*e.z*e.z;
    const double orthogonal_area = Sf.dot(e);
    const Vec3 Snon{
        Sf.x - orthogonal_area*e.x,
        Sf.y - orthogonal_area*e.y,
        Sf.z - orthogonal_area*e.z};
    const Vec3 gp{
        0.5*(grad_p.component_data(0)[o] + grad_p.component_data(0)[n]),
        0.5*(grad_p.component_data(1)[o] + grad_p.component_data(1)[n]),
        0.5*(grad_p.component_data(2)[o] + grad_p.component_data(2)[n])};
    // The face response uses the same dAU = V/A_P operator as the cell
    // velocity reconstruction. Keeping the orthogonal pressure difference
    // implicit and the non-orthogonal remainder deferred is important: the
    // continuity matrix and the conservative Rhie-Chow flux must represent
    // the same discrete operator, otherwise a residual floor can survive
    // even when the momentum equations appear converged.
    // Orthogonal pressure difference is implicit in the pressure matrix.
    // The non-orthogonal remainder is deferred explicitly using the same
    // Gauss gradient. This keeps the face flux and pressure equation on the
    // same discrete operator on arbitrary meshes.
    const double pressure_gradient_flux =
        (p(n) - p(o))/d * orthogonal_area + gp.dot(Snon);
    return rho * rfn * pressure_gradient_flux;
}

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_rhie_chow_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs = {})
;

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_rhie_chow_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& p,
    const std::vector<double>& rAU,
    double rho,
    const VelocityBoundaryConditions& bcs)
{
    std::array<std::vector<double>, 3> directional_rAU{rAU, rAU, rAU};
    return make_rhie_chow_mass_flux(
        mesh, geometry, U, p, directional_rAU, rho, bcs, {});
}

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_rhie_chow_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs)
{
    using namespace cfdx::core;
    const auto& own = mesh.ownership();
    if (p.size() != mesh.n_cells())
        throw std::invalid_argument("make_rhie_chow_mass_flux: field size mismatch");
    // Boundary pressure is passed explicitly to the gradient operator. A
    // boundary face has no neighbour coefficient: any dAU/rAU interpolation
    // must therefore use the owner-cell value rather than a ghost/zero entry.
    const auto grad_p = gauss_gradient_with_boundary(
        p, mesh, geometry, pressure_bcs);
    for (const auto& component : rAU) {
        if (component.size() != mesh.n_cells())
            throw std::invalid_argument("make_rhie_chow_mass_flux: inverse momentum diagonal size mismatch");
        for (const double value : component) {
            if (!std::isfinite(value) || value <= 0.0)
                throw std::invalid_argument(
                    "make_rhie_chow_mass_flux: inverse momentum diagonal must be finite and positive");
        }
    }

    // Momentum-weighted (Rhie-Chow) face flux:
    //
    //     phi_f = rho [ H/A |_f . Sf
    //                    - d_f (p_N - p_P)/|C_N-C_P| |Sf| ]
    //
    // The first term is linearly interpolated from the cell-centred predictor.
    // The pressure difference is evaluated directly on the face-normal line.
    // This is deliberately the same pressure operator used to build the
    // segregated pressure equation; it must never be replaced by an unrelated
    // least-squares fit of the cell velocity.
    //
    // For a fixed-pressure boundary the boundary pressure is known. For a
    // zero-gradient pressure boundary there is no normal pressure correction.
    Field<double, Location::FACE> flux(mesh.n_faces(), "phi", "kg/s", 1);
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const std::size_t o = own.owner(f);
        Vec3 Hf;
        U.get(o, Hf.x, Hf.y, Hf.z);
        const auto nr = own.neighbour(f);
        double correction = 0.0;
        if (nr >= 0) {
            const std::size_t n = static_cast<std::size_t>(nr);
            Vec3 Hn;
            U.get(n, Hn.x, Hn.y, Hn.z);
            Hf = (Hf + Hn) * 0.5;
            correction = rhie_chow_pressure_flux_internal(
                mesh, geometry, f, p, grad_p, rAU, rho);
        } else {
            const std::size_t patch = geometry.face_patch[f];
            if (patch < mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch).type == PatchType::EMPTY) {
                flux(f) = 0.0;
                continue;
            }
            if (patch < mesh.boundary().n_patches()) {
                const auto& name = mesh.boundary().patch(patch).name;
                const auto ubc = bcs.find(name);
                if (ubc != bcs.end() &&
                    ubc->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                    Hf = ubc->second.value;

                const auto pbc = pressure_bcs.find(name);
                if (pbc != pressure_bcs.end() &&
                    pbc->second.type == ScalarBoundaryType::FIXED_VALUE) {
                    const std::size_t o = own.owner(f);
                    const Vec3 Sf = geometry.face_area_vectors[f];
                    const double area = Sf.mag();
                    const Vec3 dvec =
                        geometry.face_centres[f] - geometry.cell_centres[o];
                    const double d = dvec.mag();
                    if (!(d > 0.0))
                        throw std::runtime_error("make_rhie_chow_mass_flux: degenerate boundary distance");
                    const Vec3 e{dvec.x/d, dvec.y/d, dvec.z/d};
                    const double rfn =
                        geometry.cell_volumes[o] * (
                            rAU[0][o]*e.x*e.x +
                            rAU[1][o]*e.y*e.y +
                            rAU[2][o]*e.z*e.z);
                    // For a prescribed boundary pressure the normal pressure
                    // difference is implicit. The non-orthogonal remainder
                    // requires a boundary gradient model and is therefore not
                    // fabricated here.
                    correction = rho * rfn * (pbc->second.value - p(o)) *
                        Sf.dot(e) / d;
                }
            }
        }
        flux(f) = rho * Hf.dot(geometry.face_area_vectors[f]) - correction;
    }
    return flux;
}

inline ScalarEquation assemble_momentum_component(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure_gradient_component,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& body_component,
    double effective_dynamic_viscosity,
    const ScalarBoundaryConditions& velocity_bcs,
    std::size_t component,
    bool bounded,
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* convected_field = nullptr)
{
    cfdx::core::Field<double, cfdx::core::Location::CELL> source(
        mesh.n_cells(), "momentum_source", "N/m3", 1);
    cfdx::core::Field<double, cfdx::core::Location::CELL> sp(
        mesh.n_cells(), "momentum_sp", "kg/m3/s", 1);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        if (pressure_gradient_component.dimension() != 3 ||
            pressure_gradient_component.size() != mesh.n_cells())
            throw std::invalid_argument("assemble_momentum_component: pressure gradient must be a 3-component cell field");
        source(c) = body_component(c) - pressure_gradient_component.component_data(component)[c];
        sp(c) = 0.0;
    }
    if (component >= 3)
        throw std::invalid_argument("assemble_momentum_component: component index out of range");
    return assemble_scalar_equation(
        mesh, geometry, mass_flux, effective_dynamic_viscosity,
        source, sp, velocity_bcs, bounded, nullptr, nullptr, nullptr, nullptr,
        convection_scheme, convected_field);
}

inline void apply_velocity_boundary_conditions(
    const cfdx::core::Mesh& mesh,
    const VelocityBoundaryConditions& bcs,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U)
{
    // Cell-centred finite-volume fields have no independent patch storage.
    // Fixed-value conditions are therefore enforced through the matrix and
    // the face flux; no cell overwrite is performed here.
    (void)mesh;
    (void)bcs;
    (void)U;
}

inline double sample_incompressible_probe(
    const IncompressiblePointProbe& probe,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& p)
{
    if (mesh.n_cells() == 0)
        throw std::invalid_argument("sample_incompressible_probe: mesh has no cells");
    std::size_t nearest = 0;
    double best_distance = std::numeric_limits<double>::infinity();
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        const auto delta = geometry.cell_centres[c] - probe.location;
        const double distance = delta.mag2();
        if (distance < best_distance) {
            best_distance = distance;
            nearest = c;
        }
    }
    switch (probe.field) {
    case IncompressibleProbeField::PRESSURE: return p(nearest);
    case IncompressibleProbeField::U_X: { cfdx::core::Vec3 v; U.get(nearest, v.x, v.y, v.z); return v.x; }
    case IncompressibleProbeField::U_Y: { cfdx::core::Vec3 v; U.get(nearest, v.x, v.y, v.z); return v.y; }
    case IncompressibleProbeField::U_Z: { cfdx::core::Vec3 v; U.get(nearest, v.x, v.y, v.z); return v.z; }
    case IncompressibleProbeField::U_MAGNITUDE: {
        cfdx::core::Vec3 velocity;
        U.get(nearest, velocity.x, velocity.y, velocity.z);
        return velocity.mag();
    }
    }
    throw std::invalid_argument("sample_incompressible_probe: unsupported field");
}

inline void relax_momentum_equation(
    ScalarEquation& equation,
    const cfdx::core::Vector& old_field,
    double alpha)
{
    if (!(alpha > 0.0 && alpha <= 1.0) ||
        old_field.size() != equation.rhs.size())
        throw std::invalid_argument("relax_momentum_equation: invalid relaxation");

    if (alpha == 1.0) return;

    const auto* rows = equation.matrix.row_offsets_data();
    const auto* cols = equation.matrix.columns_data();
    auto* values = equation.matrix.values_data();

    // Standard equation relaxation:
    //
    //   A u = b
    //
    // becomes
    //
    //   [A_off + A_P/alpha] u
    //       = b + (1-alpha)/alpha A_P u_old.
    //
    // This preserves the exact algebraic relation used later to construct
    // HbyA and rAU. A field blend after the solve would break that relation.
    for (std::size_t c = 0; c < old_field.size(); ++c) {
        const double a = equation.diagonal[c];
        if (!(a > 0.0) || !std::isfinite(a))
            throw std::runtime_error("relax_momentum_equation: invalid diagonal");

        double* diagonal_value = nullptr;
        for (std::uint32_t k = rows[c]; k < rows[c + 1]; ++k) {
            if (cols[k] == c) {
                diagonal_value = &values[k];
                break;
            }
        }
        if (!diagonal_value)
            throw std::runtime_error("relax_momentum_equation: missing diagonal entry");

        const double added = a * (1.0 / alpha - 1.0);
        *diagonal_value += added;
        equation.diagonal[c] += added;
        equation.rhs(c) += added * old_field(c);
    }
}

inline cfdx::core::SolverResult solve_coupled_momentum_continuity(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    ScalarEquation ex,
    ScalarEquation ey,
    ScalarEquation ez,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U_old,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& p_old,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    double rho,
    std::size_t reference_cell,
    double reference_value,
    std::size_t max_iterations,
    double tolerance,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    const std::size_t nv = 3 * nc;
    const std::size_t np = nc;
    const std::size_t n = nv + np;

    SparseMatrix A(n, n);
    Vector b(n, 0.0);

    // The momentum equations received the current pressure explicitly during
    // segregated assembly. Undo that explicit contribution here: the coupled
    // matrix contains the pressure gradient implicitly instead.
    const auto grad_p_old =
        gauss_gradient_with_boundary(p_old, mesh, geometry, pressure_bcs);
    auto add_momentum_block = [&](const ScalarEquation& eq, std::size_t component) {
        const std::size_t row_base = component * nc;
        const auto* ro = eq.matrix.row_offsets_data();
        const auto* co = eq.matrix.columns_data();
        const auto* va = eq.matrix.values_data();
        for (std::size_t r = 0; r < nc; ++r) {
            for (std::uint32_t k = ro[r]; k < ro[r + 1]; ++k)
                A.push_back(row_base + r, component * nc + co[k], va[k]);
            b(row_base + r) =
                eq.rhs(r) + grad_p_old.component_data(component)[r] * geometry.cell_volumes[r];
        }
    };
    add_momentum_block(ex, 0);
    add_momentum_block(ey, 1);
    add_momentum_block(ez, 2);

    // Move the current explicit pressure term from b back into the linear
    // system. The FV Gauss gradient is linear in the pressure field, so its
    // coefficients can be assembled face-by-face.
    for (std::size_t c = 0; c < nc; ++c) {
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        for (Offset k = 0; k < count; ++k) {
            const std::size_t f = mesh.cells().faces_data()[off + k];
            const bool owner = mesh.ownership().owner(f) == c;
            const double sign = owner ? 1.0 : -1.0;
            const Vec3 Sf = geometry.face_area_vectors[f] * sign;
            const auto nr = mesh.ownership().neighbour(f);

            if (nr >= 0) {
                const std::size_t ncell = static_cast<std::size_t>(nr);
                // Internal Gauss face pressure is the arithmetic average.
                const double coeff = 0.5;
                A.push_back(c, 3*nc + c, Sf.x * coeff);
                A.push_back(c, 3*nc + ncell, Sf.x * coeff);
                A.push_back(nc + c, 3*nc + c, Sf.y * coeff);
                A.push_back(nc + c, 3*nc + ncell, Sf.y * coeff);
                A.push_back(2*nc + c, 3*nc + c, Sf.z * coeff);
                A.push_back(2*nc + c, 3*nc + ncell, Sf.z * coeff);
            } else {
                const std::size_t patch = geometry.face_patch[f];
                const ScalarBoundaryCondition* bc = nullptr;
                ScalarBoundaryCondition local_bc;
                if (patch < mesh.boundary().n_patches()) {
                    const auto& name = mesh.boundary().patch(patch).name;
                    const auto it = pressure_bcs.find(name);
                    if (it != pressure_bcs.end()) {
                        local_bc = it->second;
                        bc = &local_bc;
                    }
                }
                if (bc && bc->type == ScalarBoundaryType::FIXED_VALUE) {
                    // Known boundary pressure is moved to the momentum RHS.
                    b(c) -= Sf.x * bc->value;
                    b(nc + c) -= Sf.y * bc->value;
                    b(2*nc + c) -= Sf.z * bc->value;
                } else {
                    // The transport-layer Gauss gradient treats zero-gradient
                    // pressure as owner-cell extrapolation. Retain that exact
                    // operator here so the coupled and segregated residuals
                    // remain comparable.
                    A.push_back(c, 3*nc + c, Sf.x);
                    A.push_back(nc + c, 3*nc + c, Sf.y);
                    A.push_back(2*nc + c, 3*nc + c, Sf.z);
                }
            }
        }
    }

    // Continuity is assembled from the same face velocity interpolation used
    // by the pressure-based solver. The pressure block is the implicit
    // momentum-weighted (Rhie-Chow) face-flux derivative.
    for (std::size_t c = 0; c < nc; ++c) {
        const std::size_t row = nv + c;
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        for (Offset k = 0; k < count; ++k) {
            const std::size_t f = mesh.cells().faces_data()[off + k];
            const bool owner = mesh.ownership().owner(f) == c;
            const double sign = owner ? 1.0 : -1.0;
            const Vec3 Sf = geometry.face_area_vectors[f] * sign;
            const auto nr = mesh.ownership().neighbour(f);

            if (nr >= 0) {
                const std::size_t ncell = static_cast<std::size_t>(nr);
                const double half_rho = 0.5 * rho;
                A.push_back(row, c, half_rho * Sf.x);
                A.push_back(row, ncell, half_rho * Sf.x);
                A.push_back(row, nc + c, half_rho * Sf.y);
                A.push_back(row, nc + ncell, half_rho * Sf.y);
                A.push_back(row, 2*nc + c, half_rho * Sf.z);
                A.push_back(row, 2*nc + ncell, half_rho * Sf.z);

                const Vec3 rawSf = geometry.face_area_vectors[f];
                const double area = rawSf.mag();
                const double d = (geometry.cell_centres[ncell] -
                                  geometry.cell_centres[c]).mag();
                const double nx = rawSf.x / area;
                const double ny = rawSf.y / area;
                const double nz = rawSf.z / area;
                const double aox = ex.diagonal[c], anx = ex.diagonal[ncell];
                const double aoy = ey.diagonal[c], any = ey.diagonal[ncell];
                const double aoz = ez.diagonal[c], anz = ez.diagonal[ncell];
                if (!(aox > 0.0) || !(anx > 0.0) ||
                    !(aoy > 0.0) || !(any > 0.0) ||
                    !(aoz > 0.0) || !(anz > 0.0) ||
                    !std::isfinite(aox) || !std::isfinite(anx) ||
                    !std::isfinite(aoy) || !std::isfinite(any) ||
                    !std::isfinite(aoz) || !std::isfinite(anz))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid momentum diagonal on face " +
                        std::to_string(f) + " cells " + std::to_string(c) + "/" +
                        std::to_string(ncell));
                // The coupled Schur complement uses the same dAU = V/A_P
                // pressure response as the segregated Rhie-Chow operator.
                const double rAUx_f = 0.5 * (
                    geometry.cell_volumes[c] / aox +
                    geometry.cell_volumes[ncell] / anx);
                const double rAUy_f = 0.5 * (
                    geometry.cell_volumes[c] / aoy +
                    geometry.cell_volumes[ncell] / any);
                const double rAUz_f = 0.5 * (
                    geometry.cell_volumes[c] / aoz +
                    geometry.cell_volumes[ncell] / anz);
                const double rfn = rAUx_f*nx*nx + rAUy_f*ny*ny + rAUz_f*nz*nz;
                if (!(rfn > 0.0) || !std::isfinite(rfn))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face rAU on face " +
                        std::to_string(f));
                const double D = rho * rfn * area / d;
                if (!(d > 0.0) || !(area > 0.0) ||
                    !std::isfinite(d) || !std::isfinite(area) ||
                    !(D > 0.0) || !std::isfinite(D)) {
                    const char* debug = std::getenv("CFDX_DEBUG_COUPLED");
                    if (debug && *debug) {
                        std::cerr << "COUPLED_FACE_DEBUG face=" << f
                                  << " cells=" << c << "/" << ncell
                                  << " area=" << area
                                  << " distance=" << d
                                  << " a=(" << aox << "," << anx
                                  << ";" << aoy << "," << any
                                  << ";" << aoz << "," << anz << ")"
                                  << " rAU=(" << rAUx_f << "," << rAUy_f
                                  << "," << rAUz_f << ")"
                                  << " rfn=" << rfn
                                  << " D=" << D << "\n";
                    }
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid pressure coefficient on face " +
                        std::to_string(f) + " cells " + std::to_string(c) + "/" +
                        std::to_string(ncell));
                }
                // phi_p = D (p_P - p_N). Only this orthogonal derivative is
                // implicit in the coupled matrix. The non-orthogonal remainder
                // is a deferred correction evaluated from p_old below, exactly
                // as in the segregated Rhie-Chow operator.
                A.push_back(row, nv + c, D);
                A.push_back(row, nv + ncell, -D);

                const Vec3 dvec = geometry.cell_centres[ncell] -
                                  geometry.cell_centres[c];
                const Vec3 e{dvec.x/d, dvec.y/d, dvec.z/d};
                const double orthogonal_area = rawSf.dot(e);
                const Vec3 Snon{
                    rawSf.x - orthogonal_area*e.x,
                    rawSf.y - orthogonal_area*e.y,
                    rawSf.z - orthogonal_area*e.z};
                const Vec3 gp{
                    0.5*(grad_p_old.component_data(0)[c] +
                         grad_p_old.component_data(0)[ncell]),
                    0.5*(grad_p_old.component_data(1)[c] +
                         grad_p_old.component_data(1)[ncell]),
                    0.5*(grad_p_old.component_data(2)[c] +
                         grad_p_old.component_data(2)[ncell])};
                const double nonorth_flux = rho * rfn * gp.dot(Snon);
                // The continuity equation is div(phi)=0 and the pressure
                // correction is written as phi = phi_U - D grad(p).
                // Therefore the explicit non-orthogonal pressure flux is
                // subtracted from the current continuity residual.
                b(row) += nonorth_flux;
            } else {
                const std::size_t patch = geometry.face_patch[f];
                if (patch < mesh.boundary().n_patches() &&
                    mesh.boundary().patch(patch).type == PatchType::EMPTY)
                    continue;
                const VelocityBoundaryCondition* vbc = nullptr;
                VelocityBoundaryCondition local_vbc;
                const ScalarBoundaryCondition* pbc = nullptr;
                ScalarBoundaryCondition local_pbc;
                if (patch < mesh.boundary().n_patches()) {
                    const auto& name = mesh.boundary().patch(patch).name;
                    const auto vit = velocity_bcs.find(name);
                    if (vit != velocity_bcs.end()) {
                        local_vbc = vit->second;
                        vbc = &local_vbc;
                    }
                    const auto pit = pressure_bcs.find(name);
                    if (pit != pressure_bcs.end()) {
                        local_pbc = pit->second;
                        pbc = &local_pbc;
                    }
                }

                if (vbc && vbc->type == VelocityBoundaryCondition::Type::FIXED_VALUE) {
                    b(row) -= rho * vbc->value.dot(Sf);
                } else {
                    A.push_back(row, c, rho * Sf.x);
                    A.push_back(row, nc + c, rho * Sf.y);
                    A.push_back(row, 2*nc + c, rho * Sf.z);
                }

                if (pbc && pbc->type == ScalarBoundaryType::FIXED_VALUE) {
                    const Vec3 rawSf = geometry.face_area_vectors[f];
                    const double area = rawSf.mag();
                    const double d = (geometry.face_centres[f] -
                                      geometry.cell_centres[c]).mag();
                    const double nx = rawSf.x / area;
                    const double ny = rawSf.y / area;
                    const double nz = rawSf.z / area;
                    const double a_x = ex.diagonal[c];
                    const double a_y = ey.diagonal[c];
                    const double a_z = ez.diagonal[c];
                    if (!(a_x > 0.0) || !(a_y > 0.0) || !(a_z > 0.0) ||
                        !std::isfinite(a_x) || !std::isfinite(a_y) || !std::isfinite(a_z))
                        throw std::runtime_error(
                            "solve_coupled_momentum_continuity: invalid boundary momentum diagonal on face " +
                            std::to_string(f) + " cell " + std::to_string(c));
                    const double rAUx = geometry.cell_volumes[c] / a_x;
                    const double rAUy = geometry.cell_volumes[c] / a_y;
                    const double rAUz = geometry.cell_volumes[c] / a_z;
                    const double rfn = rAUx*nx*nx + rAUy*ny*ny + rAUz*nz*nz;
                    if (!(rfn > 0.0) || !std::isfinite(rfn))
                        throw std::runtime_error(
                            "solve_coupled_momentum_continuity: invalid boundary face rAU on face " +
                            std::to_string(f));
                    const double D = rho * rfn * area / d;
                    if (!(D > 0.0) || !std::isfinite(D))
                        throw std::runtime_error(
                            "solve_coupled_momentum_continuity: invalid boundary pressure coefficient on face " +
                            std::to_string(f));
                    A.push_back(row, nv + c, D);
                    b(row) += D * pbc->value;
                }
            }
        }
    }

    // A pure-Neumann pressure field has one gauge null mode. Eliminate the
    // reference pressure column and replace the reference continuity row by
    // p_ref = reference_value. This keeps the coupled system deterministic.
    bool has_fixed_pressure = false;
    for (const auto& [name, bc] : pressure_bcs) {
        (void)name;
        has_fixed_pressure = has_fixed_pressure ||
            bc.type == ScalarBoundaryType::FIXED_VALUE;
    }

    if (!has_fixed_pressure) {
        // Rebuild the matrix while eliminating the known gauge column.
        SparseMatrix reduced(n, n);
        const auto* ro = A.row_offsets_data();
        const auto* co = A.columns_data();
        const auto* va = A.values_data();
        for (std::size_t r = 0; r < n; ++r) {
            if (r == nv + reference_cell) {
                reduced.push_back(r, nv + reference_cell, 1.0);
                b(r) = reference_value;
                continue;
            }
            for (std::uint32_t k = ro[r]; k < ro[r + 1]; ++k) {
                const std::size_t col = co[k];
                if (col == nv + reference_cell) {
                    b(r) -= va[k] * reference_value;
                } else {
                    reduced.push_back(r, col, va[k]);
                }
            }
        }
        reduced.finalize();
        A = std::move(reduced);
    } else {
        A.finalize();
    }

    Vector x(n, 0.0);
    for (std::size_t c = 0; c < nc; ++c) {
        x(c) = U_old.component_data(0)[c];
        x(nc + c) = U_old.component_data(1)[c];
        x(2*nc + c) = U_old.component_data(2)[c];
        x(nv + c) = p_old(c);
    }

    // Reject non-finite coefficients before entering GMRES. A bad local
    // coefficient must be diagnosed at assembly time rather than surfacing as
    // an opaque Krylov failure.
    {
        const auto* values = A.values_data();
        for (std::size_t k = 0; k < A.nnz(); ++k) {
            if (!std::isfinite(values[k]))
                throw std::runtime_error(
                    "solve_coupled_momentum_continuity: non-finite matrix coefficient at nnz " +
                    std::to_string(k));
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (!std::isfinite(b(i)))
                throw std::runtime_error(
                    "solve_coupled_momentum_continuity: non-finite RHS at row " +
                    std::to_string(i));
        }
    }

    // The coupled matrix has strongly different momentum and pressure
    // scales. CellBlockJacobi is a local exact 4x4 block preconditioner;
    // block-Schur/AMG remains a separate roadmap item.
    CellBlockJacobiPreconditioner coupled_preconditioner(nc);
    auto result = solve_gmres(
        A, b, x, 64, max_iterations, tolerance, &coupled_preconditioner);
    if (result.status != SolverStatus::CONVERGED)
        return result;

    for (std::size_t c = 0; c < nc; ++c) {
        U.component_data(0)[c] = x(c);
        U.component_data(1)[c] = x(nc + c);
        U.component_data(2)[c] = x(2*nc + c);
        p(c) = x(nv + c);
    }
    return result;
}

inline IncompressibleSolveResult solve_steady_incompressible(
    const cfdx::core::Mesh& mesh,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    const IncompressibleSolverControls& controls = {},
    const std::string& restart_path = {},
    cfdx::io::DatRestartFields restart_fields = {})
{
    using namespace cfdx::core;
    if (U.dimension() != 3 || U.size() != mesh.n_cells() ||
        p.dimension() != 1 || p.size() != mesh.n_cells())
        throw std::invalid_argument("solve_steady_incompressible: invalid fields");

    if (!restart_path.empty())
        (void)cfdx::io::read_dat_restart_fields(restart_path, mesh, U, p, restart_fields);

    validate_incompressible_controls(controls, mesh.n_cells());
    const FvGeometry geometry = build_fv_geometry(mesh);
    const double mu_eff = controls.density *
        (controls.kinematic_viscosity + controls.turbulent_viscosity);

    bool has_fixed_pressure_boundary = false;
    for (const auto& [name, bc] : pressure_bcs) {
        (void)name;
        if (bc.type == ScalarBoundaryType::FIXED_VALUE) {
            has_fixed_pressure_boundary = true;
            break;
        }
    }

    const std::size_t pcorr =
        (controls.algorithm == PressureVelocityAlgorithm::PISO ||
         controls.algorithm == PressureVelocityAlgorithm::PIMPLE)
            ? static_cast<std::size_t>(controls.coupling.n_pressure_correctors)
            : (controls.algorithm == PressureVelocityAlgorithm::FRACTIONAL_STEP
                ? static_cast<std::size_t>(controls.coupling.n_fractional_steps)
                : 1u);
    const std::size_t minimum_outer_correctors =
        controls.algorithm == PressureVelocityAlgorithm::PIMPLE
            ? static_cast<std::size_t>(controls.coupling.n_outer_correctors)
            : 1u;

    IncompressibleSolveResult result;

    auto build_hbya = [&](const ScalarEquation& eq,
                          const Vector& field,
                          const Field<double, Location::CELL>& gradp,
                          std::size_t component,
                          const std::vector<double>& rAU,
                          const std::vector<double>& rAtU) {
        std::vector<double> hbya(mesh.n_cells(), 0.0);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            double h = eq.rhs(c) + gradp.component_data(component)[c] * geometry.cell_volumes[c];
            const auto begin = eq.matrix.row_offsets_data()[c];
            const auto end = eq.matrix.row_offsets_data()[c + 1];
            for (std::uint32_t k = begin; k < end; ++k) {
                const auto col = eq.matrix.columns_data()[k];
                if (col != c)
                    h -= eq.matrix.values_data()[k] * field(col);
            }
            // SIMPLEC follows the established consistent formulation:
            // HbyA = rAU*H - (rAU-rAtU)*grad(p).
            hbyA[c] = rAU[c] * h;
            if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC)
                hbyA[c] -= (rAU[c] - rAtU[c]) *
                           gradp.component_data(component)[c] * geometry.cell_volumes[c];
            if (!std::isfinite(hbyA[c]))
                throw std::runtime_error("solve_steady_incompressible: non-finite HbyA");
        }
        return hbya;
    };

    auto make_hbya_field = [&](const std::array<std::vector<double>,3>& hbya) {
        Field<double, Location::CELL> field(
            mesh.n_cells(), "HbyA", "m/s", 3);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            field.component_data(0)[c] = hbyA[0][c];
            field.component_data(1)[c] = hbyA[1][c];
            field.component_data(2)[c] = hbyA[2][c];
        }
        return field;
    };

    // Face pressure projection coefficient.
    //
    // Momentum diagonals are integrated finite-volume coefficients A_P [N s/m].
    // rAU = 1/A_P therefore has units of m/(N s), while the velocity response
    // to a pressure gradient is dAU = V*rAU.  The same dAU operator must be
    // used by cell velocity reconstruction, Rhie-Chow face fluxes and the
    // pressure-correction matrix.
    auto directional_face_coefficient =
        [&](const std::array<std::vector<double>,3>& coeff,
            std::size_t o, std::size_t n,
            const Vec3& Sf) {
            const double area = Sf.mag();
            if (!(area > 0.0))
                throw std::runtime_error("solve_steady_incompressible: degenerate face area");
            const Vec3 nf{Sf.x/area, Sf.y/area, Sf.z/area};
            const double cx = 0.5 * (
                geometry.cell_volumes[o] * coeff[0][o] +
                geometry.cell_volumes[n] * coeff[0][n]);
            const double cy = 0.5 * (
                geometry.cell_volumes[o] * coeff[1][o] +
                geometry.cell_volumes[n] * coeff[1][n]);
            const double cz = 0.5 * (
                geometry.cell_volumes[o] * coeff[2][o] +
                geometry.cell_volumes[n] * coeff[2][n]);
            return cx*nf.x*nf.x + cy*nf.y*nf.y + cz*nf.z*nf.z;
        };

    Field<double, Location::FACE> mass_flux =
        make_mass_flux(mesh, geometry, U, controls.density, velocity_bcs);

    const char* freeze_state_env = std::getenv("CFDX_FREEZE_STATE");
    const bool freeze_state = freeze_state_env && std::string(freeze_state_env) == "1";
    Field<double, Location::FACE> phi_used_in_momentum;
    std::array<std::vector<double>, 3> frozen_rAU;
    std::array<std::vector<double>, 3> frozen_hbyA;
    Field<double, Location::CELL> frozen_grad_p;
    bool frozen_state_valid = false;

    for (std::size_t iter = 1; iter <= controls.convergence.max_iterations; ++iter) {
        const auto U_old = U;
        const auto p_old = p;

        // Exact conservative flux consumed by the momentum assembly.
        if (freeze_state)
            phi_used_in_momentum = mass_flux;

        auto grad_p = gauss_gradient_with_boundary(p, mesh, geometry, pressure_bcs);

        Field<double, Location::CELL> body_x(mesh.n_cells(), "body_x", "N/m3", 1);
        Field<double, Location::CELL> body_y(mesh.n_cells(), "body_y", "N/m3", 1);
        Field<double, Location::CELL> body_z(mesh.n_cells(), "body_z", "N/m3", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            body_x(c) = controls.density * controls.body_force.x;
            body_y(c) = controls.density * controls.body_force.y;
            body_z(c) = controls.density * controls.body_force.z;
        }

        ScalarBoundaryConditions ubc_x, ubc_y, ubc_z;
        Field<double, Location::CELL> u_x_field(mesh.n_cells(), "Ux_reconstruction", "m/s", 1);
        Field<double, Location::CELL> u_y_field(mesh.n_cells(), "Uy_reconstruction", "m/s", 1);
        Field<double, Location::CELL> u_z_field(mesh.n_cells(), "Uz_reconstruction", "m/s", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            u_x_field(c) = U.component_data(0)[c];
            u_y_field(c) = U.component_data(1)[c];
            u_z_field(c) = U.component_data(2)[c];
        }
        for (const auto& [name, bc] : velocity_bcs) {
            ubc_x[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                               ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                           bc.value.x, 0.0};
            ubc_y[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                               ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                           bc.value.y, 0.0};
            ubc_z[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                               ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                           bc.value.z, 0.0};
        }

        // Momentum predictor uses the conservative corrected flux from the
        // previous outer iteration. This preserves the SIMPLE phi state.

        auto ex = assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body_x, mu_eff, ubc_x, 0,
            controls.use_bounded_convection, controls.convection_scheme, &u_x_field);
        auto ey = assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body_y, mu_eff, ubc_y, 1,
            controls.use_bounded_convection, controls.convection_scheme, &u_y_field);
        auto ez = assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body_z, mu_eff, ubc_z, 2,
            controls.use_bounded_convection, controls.convection_scheme, &u_z_field);

        // Preserve the pre-solve velocity. HbyA is the explicit momentum
        // predictor H/A and must be reconstructed from the velocity state that
        // was used to assemble the current matrix, not from the already solved
        // momentum state. Using the post-solve vector here double-counts the
        // momentum solve and leaves the independent momentum residual O(1).
        Vector ux(mesh.n_cells()), uy(mesh.n_cells()), uz(mesh.n_cells());
        Vector ux_old(mesh.n_cells()), uy_old(mesh.n_cells()), uz_old(mesh.n_cells());
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            ux(c) = U.component_data(0)[c];
            uy(c) = U.component_data(1)[c];
            uz(c) = U.component_data(2)[c];
            ux_old(c) = ux(c);
            uy_old(c) = uy(c);
            uz_old(c) = uz(c);
        }

        cfdx::core::SolverResult rx{}, ry{}, rz{};
        double pressure_residual = std::numeric_limits<double>::infinity();
        std::size_t pressure_iterations = 0;
        if (controls.algorithm == PressureVelocityAlgorithm::COUPLED) {
            const auto coupled_result = solve_coupled_momentum_continuity(
                mesh, geometry, ex, ey, ez, U_old, p_old,
                velocity_bcs, pressure_bcs, controls.density,
                controls.pressure_reference_cell,
                controls.pressure_reference_value,
                controls.coupling.coupled_max_iterations,
                controls.coupling.coupled_linear_tolerance,
                U, p);
            if (coupled_result.status != cfdx::core::SolverStatus::CONVERGED) {
                throw std::runtime_error(
                    "solve_steady_incompressible: coupled momentum-continuity solve did not converge "
                    "(status=" + std::to_string(static_cast<int>(coupled_result.status)) +
                    ", iterations=" + std::to_string(coupled_result.iterations) +
                    ", residual=" + std::to_string(coupled_result.residual) +
                    ", relative=" + std::to_string(coupled_result.residual_relative) + ")");
            }

            // The coupled linear system is a Picard/Newton linearization of the
            // current nonlinear iterate. Keep the nonlinear field update explicit
            // and bounded; there is no hidden scalar fallback.
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                for (std::size_t d = 0; d < 3; ++d)
                    U.component_data(d)[c] =
                        U_old.component_data(d)[c] +
                        controls.coupling.alpha_u *
                        (U.component_data(d)[c] - U_old.component_data(d)[c]);
                p(c) = p_old(c) +
                    controls.coupling.alpha_p * (p(c) - p_old(c));
            }

            // The coupled solve has already enforced the discrete continuity
            // equation. Rebuild the authoritative flux from the same RC
            // interpolation used by the block matrix for diagnostics and for
            // the next nonlinear iteration.
            mass_flux = make_rhie_chow_mass_flux(
                mesh, geometry, U, p,
                std::array<std::vector<double>,3>{
                    [&] { std::vector<double> v(mesh.n_cells()); for (std::size_t c=0;c<mesh.n_cells();++c) v[c]=geometry.cell_volumes[c]/std::max(ex.diagonal[c],1e-30); return v; }(),
                    [&] { std::vector<double> v(mesh.n_cells()); for (std::size_t c=0;c<mesh.n_cells();++c) v[c]=geometry.cell_volumes[c]/std::max(ey.diagonal[c],1e-30); return v; }(),
                    [&] { std::vector<double> v(mesh.n_cells()); for (std::size_t c=0;c<mesh.n_cells();++c) v[c]=geometry.cell_volumes[c]/std::max(ez.diagonal[c],1e-30); return v; }()},
                controls.density, velocity_bcs, pressure_bcs);

            // A coupled iteration has no segregated momentum/pressure Krylov
            // sub-solves. The block residual is the authoritative linear metric.
            // The independent nonlinear momentum residual below remains the
            // final physical acceptance metric.
        } else {
        relax_momentum_equation(ex, ux_old, controls.coupling.alpha_u);
        relax_momentum_equation(ey, uy_old, controls.coupling.alpha_u);
        relax_momentum_equation(ez, uz_old, controls.coupling.alpha_u);

        // Equation relaxation is already part of the matrix. The Krylov solve
        // therefore returns the solution of the relaxed equation directly;
        // applying alpha_u again here would double-relax the predictor.
        rx = solve_scalar_equation(ex, ux, {
            controls.linear_max_iterations, controls.linear_tolerance, 1.0});
        ry = solve_scalar_equation(ey, uy, {
            controls.linear_max_iterations, controls.linear_tolerance, 1.0});
        rz = solve_scalar_equation(ez, uz, {
            controls.linear_max_iterations, controls.linear_tolerance, 1.0});

        auto require_linear_convergence = [](const char* component, const auto& solve) {
            if (solve.status != cfdx::core::SolverStatus::CONVERGED)
                throw std::runtime_error(
                    std::string("solve_steady_incompressible: ") + component +
                    " momentum solve did not converge (status=" +
                    std::to_string(static_cast<int>(solve.status)) +
                    ", iterations=" + std::to_string(solve.iterations) +
                    ", residual=" + std::to_string(solve.residual) +
                    ", relative=" + std::to_string(solve.residual_relative) + ")");
        };
        require_linear_convergence("Ux", rx);
        require_linear_convergence("Uy", ry);
        require_linear_convergence("Uz", rz);

        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            U.component_data(0)[c] = ux(c);
            U.component_data(1)[c] = uy(c);
            U.component_data(2)[c] = uz(c);
        }

        std::array<std::vector<double>,3> rAU, rAtU, hbya;
        for (std::size_t d = 0; d < 3; ++d) {
            rAU[d].resize(mesh.n_cells());
            rAtU[d].resize(mesh.n_cells());
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                const ScalarEquation* eqs[3] = {&ex, &ey, &ez};
                const double a = eqs[d]->diagonal[c];
                if (!(a > 0.0) || !std::isfinite(a))
                    throw std::runtime_error(
                        "solve_steady_incompressible: invalid momentum diagonal");
                rAU[d][c] = 1.0 / a;
                rAtU[d][c] = rAU[d][c];
                if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC) {
                    double h1 = 0.0;
                    const auto begin = eqs[d]->matrix.row_offsets_data()[c];
                    const auto end = eqs[d]->matrix.row_offsets_data()[c + 1];
                    for (std::uint32_t k = begin; k < end; ++k) {
                        const auto col = eqs[d]->matrix.columns_data()[k];
                        // Boundary contributions are assembled into the owner
                        // diagonal/source; off-diagonal entries are therefore
                        // genuine internal neighbour coefficients.
                        if (col != c)
                            h1 -= eqs[d]->matrix.values_data()[k];
                    }
                    const double denom = 1.0/rAU[d][c] - h1;
                    if (!(denom > 0.0) || !std::isfinite(denom))
                        throw std::runtime_error(
                            "solve_steady_incompressible: invalid SIMPLEC consistent diagonal");
                    rAtU[d][c] = 1.0 / denom;
                }
            }
        }

        hbyA[0] = build_hbya(ex, ux, grad_p, 0, rAU[0], rAtU[0]);
        hbyA[1] = build_hbya(ey, uy, grad_p, 1, rAU[1], rAtU[1]);
        hbyA[2] = build_hbya(ez, uz, grad_p, 2, rAU[2], rAtU[2]);
        auto HbyA = make_hbya_field(hbya);

        if (freeze_state) {
            frozen_rAU = rAU;
            frozen_hbyA = hbya;
            frozen_grad_p = grad_p;
            frozen_state_valid = true;
        }

        // rAU is the inverse integrated momentum diagonal (1/A_P).
        // Pressure-gradient terms in the integrated momentum equation therefore
        // carry the cell volume explicitly when reconstructing cell velocity.
        // phiHbyA is the authoritative predictor flux. The Rhie-Chow term
        // uses the same directional inverse momentum coefficient as the
        // pressure equation, avoiding a second, incompatible face operator.
        auto phiHbyA = make_mass_flux(
            mesh, geometry, HbyA, controls.density, velocity_bcs);

        if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC) {
            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                const auto nr = mesh.ownership().neighbour(f);
                if (nr < 0) continue;
                const std::size_t o = mesh.ownership().owner(f);
                const std::size_t n = static_cast<std::size_t>(nr);
                const auto Sf = geometry.face_area_vectors[f];
                const double area = Sf.mag();
                const double d = (geometry.cell_centres[n] - geometry.cell_centres[o]).mag();
                const Vec3 nf{Sf.x/area, Sf.y/area, Sf.z/area};
                const double dx = 0.5*(geometry.cell_volumes[o]*rAtU[0][o] +
                                       geometry.cell_volumes[n]*rAtU[0][n]) -
                                  0.5*(geometry.cell_volumes[o]*rAU[0][o] +
                                       geometry.cell_volumes[n]*rAU[0][n]);
                const double dy = 0.5*(geometry.cell_volumes[o]*rAtU[1][o] +
                                       geometry.cell_volumes[n]*rAtU[1][n]) -
                                  0.5*(geometry.cell_volumes[o]*rAU[1][o] +
                                       geometry.cell_volumes[n]*rAU[1][n]);
                const double dz = 0.5*(geometry.cell_volumes[o]*rAtU[2][o] +
                                       geometry.cell_volumes[n]*rAtU[2][n]) -
                                  0.5*(geometry.cell_volumes[o]*rAU[2][o] +
                                       geometry.cell_volumes[n]*rAU[2][n]);
                const double drn = dx*nf.x*nf.x + dy*nf.y*nf.y + dz*nf.z*nf.z;
                phiHbyA(f) += controls.density * drn *
                    (p(n)-p(o))/d * area;
            }
        }

        for (std::size_t corr = 0; corr < pcorr; ++corr) {
            const std::size_t nc = mesh.n_cells();
            SparseMatrix A(nc, nc);
            Vector b(nc, 0.0);
            std::vector<std::map<std::size_t, double>> rows(nc);
            std::vector<double> diag(nc, 0.0);
            std::vector<double> continuity(nc, 0.0);
            const auto current_grad_p =
                gauss_gradient_with_boundary(p, mesh, geometry, pressure_bcs);

            const auto* cell_faces = mesh.cells().faces_data();
            const auto* cell_offsets = mesh.cells().offsets_data();
            for (std::size_t c = 0; c < nc; ++c) {
                const Offset off = cell_offsets[c];
                const Offset count = cell_offsets[c + 1] - off;
                for (Offset k = 0; k < count; ++k) {
                    const std::size_t f = cell_faces[off+k];
                    continuity[c] += mesh.ownership().owner(f) == c
                        ? phiHbyA(f) : -phiHbyA(f);
                    if (mesh.ownership().neighbour(f) >= 0) {
                        const double phi_nonorth = rhie_chow_pressure_flux_internal(
                            mesh, geometry, f, p, current_grad_p,
                            controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU,
                            controls.density);
                        continuity[c] += mesh.ownership().owner(f) == c
                            ? -phi_nonorth : phi_nonorth;
                    }
                }
            }

            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                const auto nr = mesh.ownership().neighbour(f);
                if (nr < 0) continue;
                const std::size_t o = mesh.ownership().owner(f);
                const std::size_t n = static_cast<std::size_t>(nr);
                const Vec3 Sf = geometry.face_area_vectors[f];
                const double area = Sf.mag();
                const double d = (geometry.cell_centres[n]-geometry.cell_centres[o]).mag();
                const auto& pressure_coeff =
                    controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU;
                const double rfn = directional_face_coefficient(
                    pressure_coeff, o, n, Sf);
                const double coeff = controls.density * rfn * area / d;
                diag[o] += coeff;
                diag[n] += coeff;
                rows[o][n] -= coeff;
                rows[n][o] -= coeff;
            }

            for (std::size_t c = 0; c < nc; ++c) {
                rows[c][c] += diag[c];
                b(c) = -continuity[c];
            }

            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                if (mesh.ownership().neighbour(f) >= 0) continue;
                const std::size_t patch = geometry.face_patch[f];
                if (patch >= mesh.boundary().n_patches()) continue;
                const auto& name = mesh.boundary().patch(patch).name;
                const auto it = pressure_bcs.find(name);
                if (it == pressure_bcs.end() ||
                    it->second.type != ScalarBoundaryType::FIXED_VALUE) continue;
                const std::size_t o = mesh.ownership().owner(f);
                const double distance =
                    (geometry.face_centres[f]-geometry.cell_centres[o]).mag();
                const Vec3 Sf = geometry.face_area_vectors[f];
                const double area = Sf.mag();
                const auto& pressure_coeff =
                    controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU;
                const double rfn =
                    pressure_coeff[0][o]*(Sf.x/area)*(Sf.x/area) +
                    pressure_coeff[1][o]*(Sf.y/area)*(Sf.y/area) +
                    pressure_coeff[2][o]*(Sf.z/area)*(Sf.z/area);
                rows[o][o] += controls.density * rfn * area / distance;
            }

            if (!has_fixed_pressure_boundary) {
                const std::size_t ref = controls.pressure_reference_cell;
                for (std::size_t row = 0; row < nc; ++row) {
                    if (row == ref) continue;
                    // Keep the column out of the reduced equations. The
                    // reference row itself is the gauge equation.
                    rows[row].erase(ref);
                }
                rows[ref].clear();
                rows[ref][ref] = 1.0;
                b(ref) = 0.0;
            }

            for (std::size_t row = 0; row < nc; ++row)
                for (const auto& [col, value] : rows[row])
                    A.push_back(row, col, value);
            A.finalize();

            Vector p_corr(nc, 0.0);
            const auto rp = has_fixed_pressure_boundary
                ? solve_cg(A, b, p_corr, controls.linear_max_iterations, controls.linear_tolerance)
                : solve_gmres(A, b, p_corr, 64,
                              controls.linear_max_iterations, controls.linear_tolerance);
            pressure_residual = rp.residual_relative;
            pressure_iterations = rp.iterations;
            if (rp.status != cfdx::core::SolverStatus::CONVERGED)
                throw std::runtime_error(
                    "solve_steady_incompressible: pressure-correction solve did not converge "
                    "(status=" + std::to_string(static_cast<int>(rp.status)) +
                    ", iterations=" + std::to_string(rp.iterations) +
                    ", residual=" + std::to_string(rp.residual) +
                    ", relative=" + std::to_string(rp.residual_relative) + ")");

            // Relax the physical pressure exactly once. The reference is a
            // gauge: enforce it by a uniform shift, never by overwriting one
            // cell and creating an artificial local pressure jump.
            const double alpha_p = controls.coupling.alpha_p;
            for (std::size_t c = 0; c < nc; ++c)
                p(c) += alpha_p * p_corr(c);
            if (!has_fixed_pressure_boundary) {
                const double shift =
                    controls.pressure_reference_value - p(controls.pressure_reference_cell);
                for (std::size_t c = 0; c < nc; ++c)
                    p(c) += shift;
            }

            auto corrected_grad_p =
                gauss_gradient_with_boundary(p, mesh, geometry, pressure_bcs);

            // Rebuild the corrected velocity from the same HbyA/rAtU
            // operator used by the pressure equation. This is the momentum
            // correction; no least-squares flux fitting is permitted.
            for (std::size_t c = 0; c < nc; ++c) {
                U.component_data(0)[c] =
                    HbyA.component_data(0)[c] - rAtU[0][c] * corrected_grad_p.component_data(0)[c] * geometry.cell_volumes[c];
                U.component_data(1)[c] =
                    HbyA.component_data(1)[c] - rAtU[1][c] * corrected_grad_p.component_data(1)[c] * geometry.cell_volumes[c];
                U.component_data(2)[c] =
                    HbyA.component_data(2)[c] - rAtU[2][c] * corrected_grad_p.component_data(2)[c] * geometry.cell_volumes[c];
            }

            // The conservative face flux is rebuilt from the same
            // momentum-weighted interpolation used by the pressure equation.
            // This is the single authoritative phi state for continuity.
            mass_flux = make_rhie_chow_mass_flux(
                mesh, geometry, HbyA, p,
                controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU,
                controls.density, velocity_bcs, pressure_bcs);

            if (corr + 1 < pcorr) {
                // PISO's next pressure correction is driven by the current
                // conservative face flux, not by the original predictor.
                // The momentum matrix is intentionally frozen inside the
                // inner PISO loop; the updated flux is the split-operator
                // correction that carries the first pressure solve into the
                // next one. Rebuilding HbyA here would incorrectly restart
                // the correction sequence from the same predictor.
                phiHbyA = mass_flux;
            }
        }

        }
        
        // Independent final momentum residual: evaluate the actual corrected
        // state against freshly assembled equations. Do not use a U reconstructed
        // from the conservative flux as the residual state.
        Field<double, Location::CELL> final_ux_field(mesh.n_cells(), "Ux_final", "m/s", 1);
        Field<double, Location::CELL> final_uy_field(mesh.n_cells(), "Uy_final", "m/s", 1);
        Field<double, Location::CELL> final_uz_field(mesh.n_cells(), "Uz_final", "m/s", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            final_ux_field(c) = U.component_data(0)[c];
            final_uy_field(c) = U.component_data(1)[c];
            final_uz_field(c) = U.component_data(2)[c];
        }
        auto final_grad_p = gauss_gradient_with_boundary(p, mesh, geometry, pressure_bcs);
        auto final_ex = assemble_momentum_component(
            mesh, geometry, mass_flux, final_grad_p, body_x, mu_eff, ubc_x, 0,
            controls.use_bounded_convection, controls.convection_scheme, &final_ux_field);
        auto final_ey = assemble_momentum_component(
            mesh, geometry, mass_flux, final_grad_p, body_y, mu_eff, ubc_y, 1,
            controls.use_bounded_convection, controls.convection_scheme, &final_uy_field);
        auto final_ez = assemble_momentum_component(
            mesh, geometry, mass_flux, final_grad_p, body_z, mu_eff, ubc_z, 2,
            controls.use_bounded_convection, controls.convection_scheme, &final_uz_field);

        Vector final_ux(mesh.n_cells()), final_uy(mesh.n_cells()), final_uz(mesh.n_cells());
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            final_ux(c) = U.component_data(0)[c];
            final_uy(c) = U.component_data(1)[c];
            final_uz(c) = U.component_data(2)[c];
        }
        const double final_momentum_residual = std::max({
            scalar_equation_residual_inf(final_ex, final_ux),
            scalar_equation_residual_inf(final_ey, final_uy),
            scalar_equation_residual_inf(final_ez, final_uz)
        });

        // Recompute the cell-wise residuals explicitly so a large independent
        // residual can be localized to a component and to an interior/boundary
        // cell. This is intentionally redundant with scalar_equation_residual_inf:
        // the purpose is forensic diagnostics, not a second numerical path.
        const auto cell_has_boundary_face = [&](std::size_t cell) {
            const Offset off = mesh.cells().offsets_data()[cell];
            const Offset count = mesh.cells().offsets_data()[cell + 1] - off;
            for (Offset k = 0; k < count; ++k) {
                const std::size_t face = mesh.cells().faces_data()[off + k];
                if (mesh.ownership().neighbour(face) < 0) {
                    const std::size_t patch = geometry.face_patch[face];
                    if (patch < mesh.boundary().n_patches() &&
                        mesh.boundary().patch(patch).type == PatchType::EMPTY)
                        continue;
                    return true;
                }
            }
            return false;
        };
        struct ResidualDiagnostic {
            double global = 0.0;
            double interior = 0.0;
            double boundary = 0.0;
            std::size_t max_cell = 0;
            double no_pressure = 0.0;
            double pressure_contribution = 0.0;
            std::string patch;
        };
        const auto residual_diagnostics =
            [&](const ScalarEquation& equation,
                const Vector& solution,
                std::size_t component) {
                ResidualDiagnostic diagnostic;
                for (std::size_t row = 0; row < mesh.n_cells(); ++row) {
                    double ax = 0.0;
                    const auto begin = equation.matrix.row_offsets_data()[row];
                    const auto end = equation.matrix.row_offsets_data()[row + 1];
                    for (std::uint32_t k = begin; k < end; ++k)
                        ax += equation.matrix.values_data()[k] *
                              solution(equation.matrix.columns_data()[k]);
                    const double full_signed = ax - equation.rhs(row);
                    const double residual = std::abs(full_signed);
                    if (residual > diagnostic.global) {
                        diagnostic.global = residual;
                        diagnostic.max_cell = row;
                    }
                    if (cell_has_boundary_face(row))
                        diagnostic.boundary = std::max(diagnostic.boundary, residual);
                    else
                        diagnostic.interior = std::max(diagnostic.interior, residual);

                    // final_ex is assembled with source = body - grad(p).
                    // Remove the pressure source only for diagnosis, never for
                    // the acceptance metric, to distinguish pressure/transport
                    // inconsistency from boundary transport inconsistency.
                    const double pressure = final_grad_p.component_data(component)[row] *
                                            geometry.cell_volumes[row];
                    const double transport_signed = full_signed - pressure;
                    diagnostic.no_pressure =
                        std::max(diagnostic.no_pressure, std::abs(transport_signed));
                    diagnostic.pressure_contribution =
                        std::max(diagnostic.pressure_contribution, std::abs(pressure));
                }
                const std::size_t row = diagnostic.max_cell;
                const Offset off = mesh.cells().offsets_data()[row];
                const Offset count = mesh.cells().offsets_data()[row + 1] - off;
                for (Offset k = 0; k < count; ++k) {
                    const std::size_t face = mesh.cells().faces_data()[off + k];
                    if (mesh.ownership().neighbour(face) < 0) {
                        const std::size_t patch = geometry.face_patch[face];
                        if (patch < mesh.boundary().n_patches()) {
                            diagnostic.patch = mesh.boundary().patch(patch).name;
                            break;
                        }
                    }
                }
                return diagnostic;
            };
        const auto rx_diag = residual_diagnostics(final_ex, final_ux, 0);
        const auto ry_diag = residual_diagnostics(final_ey, final_uy, 1);
        const auto rz_diag = residual_diagnostics(final_ez, final_uz, 2);

        double pressure_gradient_linf = 0.0;
        double pressure_gradient_l2_sum = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            const double gx = final_grad_p.component_data(0)[c];
            const double gy = final_grad_p.component_data(1)[c];
            const double gz = final_grad_p.component_data(2)[c];
            const double gp2 = gx*gx + gy*gy + gz*gz;
            pressure_gradient_linf = std::max(pressure_gradient_linf, std::sqrt(gp2));
            pressure_gradient_l2_sum += gp2;
        }
        const double pressure_gradient_l2 =
            std::sqrt(pressure_gradient_l2_sum /
                      std::max<std::size_t>(mesh.n_cells(), 1));

        double l1 = 0.0, linf = 0.0;
        const auto* cell_faces = mesh.cells().faces_data();
        const auto* cell_offsets = mesh.cells().offsets_data();
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            double div = 0.0;
            const Offset off = cell_offsets[c];
            const Offset count = cell_offsets[c+1] - off;
            for (Offset k = 0; k < count; ++k) {
                const std::size_t f = cell_faces[off+k];
                div += mesh.ownership().owner(f) == c ? mass_flux(f) : -mass_flux(f);
            }
            l1 += std::abs(div);
            linf = std::max(linf, std::abs(div));
        }

        double velocity_change_inf = 0.0, pressure_change_inf = 0.0;
        double velocity_scale = 1.0, pressure_scale = 1.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            for (std::size_t d = 0; d < 3; ++d) {
                velocity_change_inf = std::max(
                    velocity_change_inf,
                    std::abs(U.component_data(d)[c]-U_old.component_data(d)[c]));
                velocity_scale = std::max(
                    velocity_scale, std::abs(U.component_data(d)[c]));
            }
            pressure_change_inf = std::max(
                pressure_change_inf, std::abs(p(c)-p_old(c)));
            pressure_scale = std::max(pressure_scale, std::abs(p(c)));
        }
        velocity_change_inf /= velocity_scale;
        pressure_change_inf /= pressure_scale;

        double momentum_rhs_scale = 1.0;
        for (const auto* eq : {&final_ex, &final_ey, &final_ez})
            for (std::size_t c = 0; c < mesh.n_cells(); ++c)
                momentum_rhs_scale = std::max(momentum_rhs_scale, std::abs(eq->rhs(c)));

        IncompressibleIteration h;
        h.iteration = iter;
        h.momentum_residual =
            controls.algorithm == PressureVelocityAlgorithm::COUPLED
                ? final_momentum_residual / momentum_rhs_scale
                : std::max({rx.residual_relative, ry.residual_relative, rz.residual_relative});
        h.pressure_residual = pressure_residual;
        h.continuity_l1 = l1;
        h.continuity_linf = linf;
        h.momentum_equation_residual = final_momentum_residual;
        h.momentum_equation_residual_relative =
            final_momentum_residual / momentum_rhs_scale;
        const double domain_volume =
            std::accumulate(geometry.cell_volumes.begin(), geometry.cell_volumes.end(), 0.0);
        const double characteristic_area =
            std::max(std::pow(domain_volume, 2.0/3.0), 1e-30);
        h.continuity_normalized =
            linf / std::max(controls.density*velocity_scale*characteristic_area, 1e-30);
        if (controls.algorithm == PressureVelocityAlgorithm::COUPLED)
            h.pressure_residual = h.continuity_normalized;
        h.velocity_change_inf = velocity_change_inf;
        h.pressure_change_inf = pressure_change_inf;
        h.momentum_linear_iterations =
            controls.algorithm == PressureVelocityAlgorithm::COUPLED
                ? 0 : std::max({rx.iterations, ry.iterations, rz.iterations});
        h.pressure_linear_iterations =
            controls.algorithm == PressureVelocityAlgorithm::COUPLED ? 0 : pressure_iterations;
        h.corrected_flux_continuity_linf = linf;
        auto reconstructed_flux =
            make_mass_flux(mesh, geometry, U, controls.density, velocity_bcs);
        double reconstructed_linf = 0.0;
        double flux_mismatch_linf = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            double div_reconstructed = 0.0;
            const Offset off = cell_offsets[c];
            const Offset count = cell_offsets[c+1] - off;
            for (Offset k = 0; k < count; ++k) {
                const std::size_t f = cell_faces[off+k];
                div_reconstructed +=
                    mesh.ownership().owner(f) == c ? reconstructed_flux(f) : -reconstructed_flux(f);
                flux_mismatch_linf = std::max(
                    flux_mismatch_linf,
                    std::abs(mass_flux(f) - reconstructed_flux(f)));
            }
            reconstructed_linf = std::max(reconstructed_linf, std::abs(div_reconstructed));
        }
        h.reconstructed_velocity_continuity_linf = reconstructed_linf;
        h.flux_velocity_mismatch_linf = flux_mismatch_linf;
        h.momentum_equation_residual_components =
            std::array<double, 3>{rx_diag.global, ry_diag.global, rz_diag.global};
        h.momentum_equation_residual_internal =
            std::max({rx_diag.interior, ry_diag.interior, rz_diag.interior});
        h.momentum_equation_residual_boundary =
            std::max({rx_diag.boundary, ry_diag.boundary, rz_diag.boundary});
        h.pressure_gradient_linf = pressure_gradient_linf;
        h.pressure_gradient_l2 = pressure_gradient_l2;
        const ResidualDiagnostic* worst_diag = &rx_diag;
        if (ry_diag.global > worst_diag->global) worst_diag = &ry_diag;
        if (rz_diag.global > worst_diag->global) worst_diag = &rz_diag;
        h.momentum_residual_cell = worst_diag->max_cell;
        h.momentum_residual_no_pressure = worst_diag->no_pressure;
        h.momentum_pressure_contribution = worst_diag->pressure_contribution;
        h.momentum_residual_patch = worst_diag->patch;

        // Optional forensic microscope for one cell. This is deliberately
        // disabled unless CFDX_DEBUG_CELL is set, so production/CI output is
        // unchanged. Use CFDX_DEBUG_CELL=<index> for a fixed cell or
        // CFDX_DEBUG_CELL=auto to inspect the independently reassembled
        // worst-residual cell. The dump is emitted once per algorithm at
        // convergence or at the configured iteration limit.
        const char* debug_cell_env = std::getenv("CFDX_DEBUG_CELL");
        bool debug_cell_enabled = false;
        bool debug_cell_auto = false;
        std::size_t debug_cell = 0;
        if (debug_cell_env && *debug_cell_env) {
            try {
                const std::string value(debug_cell_env);
                if (value == "auto") {
                    debug_cell_auto = true;
                    debug_cell_enabled = true;
                } else {
                    std::size_t parsed = 0;
                    debug_cell = std::stoull(value, &parsed);
                    debug_cell_enabled =
                        parsed == value.size() && debug_cell < mesh.n_cells();
                }
            } catch (...) {
                debug_cell_enabled = false;
            }
        }
        if (debug_cell_auto)
            debug_cell = h.momentum_residual_cell;

        const bool converged_now =
            iter >= minimum_outer_correctors && iter > 1 &&
            h.momentum_residual <= controls.convergence.relative_tolerance &&
            h.momentum_equation_residual_relative <= controls.convergence.relative_tolerance &&
            h.pressure_residual <= controls.convergence.relative_tolerance &&
            h.continuity_linf <= controls.convergence.continuity_tolerance &&
            h.velocity_change_inf <= controls.convergence.relative_tolerance &&
            h.pressure_change_inf <= controls.convergence.relative_tolerance;
        if (debug_cell_enabled && (converged_now || iter == controls.convergence.max_iterations)) {
            std::cerr << "\n=== CFDX MOMENTUM MICROSCOPE cell=" << debug_cell
                      << " cell_source=" << (debug_cell_auto ? "auto_worst" : "explicit")
                      << " algorithm=" << static_cast<int>(controls.algorithm)
                      << " iteration=" << iter << " ===\n";
            const auto dump_equation = [&](const char* name,
                                           const ScalarEquation& eq,
                                           const Vector& solution) {
                double ax = 0.0;
                const auto begin = eq.matrix.row_offsets_data()[debug_cell];
                const auto end = eq.matrix.row_offsets_data()[debug_cell + 1];
                std::cerr << name << " matrix row: rhs=" << eq.rhs(debug_cell)
                          << " diagonal=" << eq.diagonal[debug_cell]
                          << " nnz=" << (end - begin) << "\n";
                for (std::uint32_t k = begin; k < end; ++k) {
                    const auto col = eq.matrix.columns_data()[k];
                    const auto value = eq.matrix.values_data()[k];
                    ax += value * solution(col);
                    std::cerr << "  col=" << col << " a=" << value
                              << " x=" << solution(col)
                              << " ax=" << value * solution(col) << "\n";
                }
                const double residual = ax - eq.rhs(debug_cell);
                std::cerr << "  exact_row_ax=" << ax
                          << " exact_row_residual=" << residual
                          << " abs=" << std::abs(residual) << "\n";
            };
            dump_equation("Ux", final_ex, final_ux);
            dump_equation("Uy", final_ey, final_uy);
            dump_equation("Uz", final_ez, final_uz);

            // Reconstruct the same momentum predictor used by the pressure-
            // velocity coupling from the final, independently assembled
            // momentum rows. This is diagnostic only; it does not alter the
            // converged state. dAU = V/A_P is the exact pressure-response
            // coefficient used by the segregated Rhie-Chow operator.
            std::array<std::vector<double>, 3> diag_rAU;
            std::array<std::vector<double>, 3> hbyA;
            const ScalarEquation* final_eqs[3] = {&final_ex, &final_ey, &final_ez};
            const Vector* final_u[3] = {&final_ux, &final_uy, &final_uz};
            for (std::size_t d = 0; d < 3; ++d) {
                diag_rAU[d].resize(mesh.n_cells());
                hbyA[d].resize(mesh.n_cells());
                for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                    const double ap = final_eqs[d]->diagonal[c];
                    if (!(ap > 0.0) || !std::isfinite(ap))
                        throw std::runtime_error(
                            "CFDX momentum microscope: invalid final momentum diagonal");
                    diag_rAU[d][c] = 1.0 / ap;
                    double h_explicit = final_eqs[d]->rhs(c) +
                        final_grad_p.component_data(d)[c] * geometry.cell_volumes[c];
                    const auto begin = final_eqs[d]->matrix.row_offsets_data()[c];
                    const auto end = final_eqs[d]->matrix.row_offsets_data()[c + 1];
                    for (std::uint32_t k = begin; k < end; ++k) {
                        const auto col = final_eqs[d]->matrix.columns_data()[k];
                        if (col != c)
                            h_explicit -=
                                final_eqs[d]->matrix.values_data()[k] * (*final_u[d])(col);
                    }
                    hbyA[d][c] = diag_rAU[d][c] * h_explicit;
                }
            }

            const auto reconstructed_flux =
                make_mass_flux(mesh, geometry, U, controls.density, velocity_bcs);
            double local_div_authoritative = 0.0;
            double local_div_reconstructed = 0.0;
            double max_delta_phi = 0.0;
            double max_delta_grad_p = 0.0;
            double max_delta_u_face = 0.0;
            double max_momentum_pressure_flux = 0.0;
            double max_correction_pressure_flux = 0.0;
            double max_rhie_chow_actual = 0.0;
            double max_rhie_chow_expected = 0.0;
            double max_rhie_chow_mismatch = 0.0;
            if (freeze_state && frozen_state_valid) {
                double phi_delta_linf = 0.0;
                double hbyA_delta_linf = 0.0;
                std::size_t hbyA_delta_cell = 0;
                for (std::size_t f = 0; f < mesh.n_faces(); ++f)
                    phi_delta_linf = std::max(phi_delta_linf,
                        std::abs(mass_flux(f) - phi_used_in_momentum(f)));
                for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                    double delta = 0.0;
                    for (std::size_t d = 0; d < 3; ++d)
                        delta = std::max(delta, std::abs(hbyA[d][c] - frozen_hbyA[d][c]));
                    if (delta > hbyA_delta_linf) {
                        hbyA_delta_linf = delta;
                        hbyA_delta_cell = c;
                    }
                }
                const std::size_t comparison_cell =
                    debug_cell == 33 && mesh.n_cells() > 70 ? 70 : 33;
                auto dump_state_cell = [&](std::size_t cell, const char* label) {
                    if (cell >= mesh.n_cells()) return;
                    double phi_delta_l1 = 0.0;
                    double phi_used_l1 = 0.0;
                    double phi_reassembled_l1 = 0.0;
                    double hbyA_delta = 0.0;
                    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                        const auto owner = mesh.ownership().owner(f);
                        const auto neighbour = mesh.ownership().neighbour(f);
                        if (owner != cell && (neighbour < 0 || static_cast<std::size_t>(neighbour) != cell)) continue;
                        const double sign = owner == cell ? 1.0 : -1.0;
                        const double used = sign * phi_used_in_momentum(f);
                        const double current = sign * mass_flux(f);
                        phi_used_l1 += std::abs(used);
                        phi_reassembled_l1 += std::abs(current);
                        phi_delta_l1 += std::abs(used-current);
                    }
                    for (std::size_t d = 0; d < 3; ++d)
                        hbyA_delta = std::max(hbyA_delta,
                            std::abs(hbyA[d][cell] - frozen_hbyA[d][cell]));
                    std::cerr << "FROZEN_STATE_CELL label=" << label
                              << " cell=" << cell
                              << " phi_used_in_momentum_assembly=" << phi_used_l1
                              << " phi_used_in_reassembly=" << phi_reassembled_l1
                              << " phi_delta_l1=" << phi_delta_l1
                              << " HbyA_delta_linf=" << hbyA_delta
                              << " frozen_rAU=(" << frozen_rAU[0][cell] << ","
                              << frozen_rAU[1][cell] << "," << frozen_rAU[2][cell] << ")"
                              << " final_rAU=(" << diag_rAU[0][cell] << ","
                              << diag_rAU[1][cell] << "," << diag_rAU[2][cell] << ")\n";
                };
                std::cerr << "FROZEN_STATE_COMPARE phi_global_linf=" << phi_delta_linf
                          << " HbyA_global_linf=" << hbyA_delta_linf
                          << " HbyA_max_cell=" << hbyA_delta_cell << "\n";
                dump_state_cell(debug_cell, "debug");
                if (comparison_cell != debug_cell)
                    dump_state_cell(comparison_cell, "comparison");
            }
            const Offset off = mesh.cells().offsets_data()[debug_cell];
            const Offset count = mesh.cells().offsets_data()[debug_cell + 1] - off;
            std::cerr << "faces=" << count << "\n";
            for (Offset k = 0; k < count; ++k) {
                const std::size_t face = mesh.cells().faces_data()[off + k];
                const auto owner = mesh.ownership().owner(face);
                const auto neighbour = mesh.ownership().neighbour(face);
                const bool owner_side = owner == debug_cell;
                const double sign = owner_side ? 1.0 : -1.0;
                const double phi_auth = sign * mass_flux(face);
                const double phi_reconstructed = sign * reconstructed_flux(face);
                local_div_authoritative += phi_auth;
                local_div_reconstructed += phi_reconstructed;
                const Vec3 Sf = geometry.face_area_vectors[face] * sign;
                std::size_t patch = geometry.face_patch[face];
                std::string patch_name = "INTERNAL";
                int patch_type = -1;
                if (neighbour < 0 && patch < mesh.boundary().n_patches()) {
                    patch_name = mesh.boundary().patch(patch).name;
                    patch_type = static_cast<int>(mesh.boundary().patch(patch).type);
                }

                Vec3 u_face_auth{0.0, 0.0, 0.0};
                Vec3 u_face_reconstructed{0.0, 0.0, 0.0};
                Vec3 grad_face_momentum{0.0, 0.0, 0.0};
                Vec3 grad_face_correction{0.0, 0.0, 0.0};
                Vec3 hbyA_face{0.0, 0.0, 0.0};
                double dAU_n = 0.0;
                double pressure_flux_correction = 0.0;
                double rhie_chow_expected = 0.0;
                double rhie_chow_actual = 0.0;

                if (neighbour >= 0) {
                    const std::size_t ncell = static_cast<std::size_t>(neighbour);
                    u_face_auth = Vec3{
                        0.5 * (final_ux(debug_cell) + final_ux(ncell)),
                        0.5 * (final_uy(debug_cell) + final_uy(ncell)),
                        0.5 * (final_uz(debug_cell) + final_uz(ncell))};
                    hbyA_face = Vec3{
                        0.5 * (hbyA[0][debug_cell] + hbyA[0][ncell]),
                        0.5 * (hbyA[1][debug_cell] + hbyA[1][ncell]),
                        0.5 * (hbyA[2][debug_cell] + hbyA[2][ncell])};
                    // Build the local face frame from the debug-cell side.
                    // For a neighbour-side debug cell, the global owner-oriented
                    // face vector and the centre-to-centre vector have opposite
                    // orientations. Mixing those frames was the source of the
                    // NaN/incorrect pressure-flux microscope records and could
                    // make a valid internal face look non-commuting.
                    const std::size_t other_cell = owner_side
                        ? ncell : owner;
                    const Vec3 dvec = geometry.cell_centres[other_cell] -
                                      geometry.cell_centres[debug_cell];
                    const double d = dvec.mag();
                    if (!(d > 0.0) || !std::isfinite(d))
                        throw std::runtime_error(
                            "CFDX momentum microscope: degenerate internal face distance");
                    const Vec3 e{dvec.x/d, dvec.y/d, dvec.z/d};
                    const double cx = 0.5 * (
                        geometry.cell_volumes[debug_cell] * diag_rAU[0][debug_cell] +
                        geometry.cell_volumes[other_cell] * diag_rAU[0][other_cell]);
                    const double cy = 0.5 * (
                        geometry.cell_volumes[debug_cell] * diag_rAU[1][debug_cell] +
                        geometry.cell_volumes[other_cell] * diag_rAU[1][other_cell]);
                    const double cz = 0.5 * (
                        geometry.cell_volumes[debug_cell] * diag_rAU[2][debug_cell] +
                        geometry.cell_volumes[other_cell] * diag_rAU[2][other_cell]);
                    dAU_n = cx*e.x*e.x + cy*e.y*e.y + cz*e.z*e.z;
                    const Vec3 gp{
                        0.5 * (final_grad_p.component_data(0)[debug_cell] +
                               final_grad_p.component_data(0)[other_cell]),
                        0.5 * (final_grad_p.component_data(1)[debug_cell] +
                               final_grad_p.component_data(1)[other_cell]),
                        0.5 * (final_grad_p.component_data(2)[debug_cell] +
                               final_grad_p.component_data(2)[other_cell])};
                    // Compare the authoritative face-flux correction against
                    // the exact Rhie-Chow operator used by the pressure equation.
                    // This is deliberately evaluated from HbyA, not the final
                    // velocity, because phi = HbyA_f.Sf - RC(p) is the
                    // conservative pressure-velocity coupling operator.
                    // Use the exact inverse momentum diagonals reconstructed
                    // above from the final momentum rows.  The microscope is
                    // diagnostic only and must not introduce a separate
                    // SIMPLEC/SIMPLE coefficient path.
                    rhie_chow_expected = rhie_chow_pressure_flux_internal(
                        mesh, geometry, face, p, final_grad_p,
                        diag_rAU, controls.density);
                    const double hbyA_flux = controls.density * hbyA_face.dot(Sf);
                    rhie_chow_actual = hbyA_flux - phi_auth;
                    max_rhie_chow_actual = std::max(
                        max_rhie_chow_actual, std::abs(rhie_chow_actual));
                    max_rhie_chow_expected = std::max(
                        max_rhie_chow_expected, std::abs(rhie_chow_expected));
                    max_rhie_chow_mismatch = std::max(
                        max_rhie_chow_mismatch,
                        std::abs(rhie_chow_actual - rhie_chow_expected));
                    const Vec3 localSf = Sf;
                    const double orthogonal_area = localSf.dot(e);
                    const Vec3 Snon{
                        localSf.x - orthogonal_area*e.x,
                        localSf.y - orthogonal_area*e.y,
                        localSf.z - orthogonal_area*e.z};
                    const double normal_dp = (p(other_cell) - p(debug_cell)) / d;
                    grad_face_correction =
                        e * normal_dp + gp - e * gp.dot(e);
                    grad_face_momentum = gp;
                    pressure_flux_correction =
                        controls.density * dAU_n *
                        (normal_dp * orthogonal_area + gp.dot(Snon));
                    u_face_reconstructed = hbyA_face -
                        Vec3{
                            dAU_n * gp.x,
                            dAU_n * gp.y,
                            dAU_n * gp.z};
                } else {
                    if (patch < mesh.boundary().n_patches()) {
                        const auto& name = mesh.boundary().patch(patch).name;
                        const auto pit = pressure_bcs.find(name);
                        const auto uit = velocity_bcs.find(name);
                        if (uit != velocity_bcs.end() &&
                            uit->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                            u_face_auth = uit->second.value;
                        else
                            u_face_auth = Vec3{
                                final_ux(debug_cell), final_uy(debug_cell), final_uz(debug_cell)};
                        if (pit != pressure_bcs.end() &&
                            pit->second.type == ScalarBoundaryType::FIXED_VALUE) {
                            const Vec3 Sf_raw = geometry.face_area_vectors[face];
                            const Vec3 dvec = geometry.face_centres[face] -
                                              geometry.cell_centres[debug_cell];
                            const double d = dvec.mag();
                            const Vec3 e{dvec.x/d, dvec.y/d, dvec.z/d};
                            dAU_n =
                                geometry.cell_volumes[debug_cell] * (
                                    diag_rAU[0][debug_cell]*e.x*e.x +
                                    diag_rAU[1][debug_cell]*e.y*e.y +
                                    diag_rAU[2][debug_cell]*e.z*e.z);
                            const double normal_dp =
                                (pit->second.value - p(debug_cell)) / d;
                            pressure_flux_correction =
                                controls.density * dAU_n * normal_dp * Sf_raw.dot(e);
                            grad_face_correction = e * normal_dp;
                        } else {
                            grad_face_correction = Vec3{
                                final_grad_p.component_data(0)[debug_cell],
                                final_grad_p.component_data(1)[debug_cell],
                                final_grad_p.component_data(2)[debug_cell]};
                        }
                        grad_face_momentum = Vec3{
                            final_grad_p.component_data(0)[debug_cell],
                            final_grad_p.component_data(1)[debug_cell],
                            final_grad_p.component_data(2)[debug_cell]};
                        hbyA_face = Vec3{
                            hbyA[0][debug_cell],
                            hbyA[1][debug_cell],
                            hbyA[2][debug_cell]};
                        const Vec3 gp = grad_face_momentum;
                        u_face_reconstructed = hbyA_face -
                            Vec3{dAU_n*gp.x, dAU_n*gp.y, dAU_n*gp.z};
                    }
                }

                const double momentum_pressure_flux = grad_face_momentum.dot(Sf);
                max_delta_phi = std::max(
                    max_delta_phi, std::abs(phi_auth - phi_reconstructed));
                max_delta_grad_p = std::max(
                    max_delta_grad_p,
                    (grad_face_momentum - grad_face_correction).mag());
                max_delta_u_face = std::max(
                    max_delta_u_face, (u_face_auth - u_face_reconstructed).mag());
                max_momentum_pressure_flux = std::max(
                    max_momentum_pressure_flux, std::abs(momentum_pressure_flux));
                max_correction_pressure_flux = std::max(
                    max_correction_pressure_flux, std::abs(pressure_flux_correction));

                std::cerr << "  face=" << face
                          << " owner=" << owner
                          << " neighbour=" << neighbour
                          << " patch=" << patch_name
                          << " patch_type=" << patch_type
                          << " Sf=(" << Sf.x << "," << Sf.y << "," << Sf.z << ")"
                          << " phi_authoritative=" << phi_auth
                          << " phi_reconstructed_from_Uface=" << phi_reconstructed
                          << " delta_phi=" << (phi_auth - phi_reconstructed)
                          << " U_face_authoritative=(" << u_face_auth.x << ","
                          << u_face_auth.y << "," << u_face_auth.z << ")"
                          << " U_face_reconstructed_HbyA_dAUgradp=("
                          << u_face_reconstructed.x << "," << u_face_reconstructed.y
                          << "," << u_face_reconstructed.z << ")"
                          << " delta_U_face=" << (u_face_auth-u_face_reconstructed).mag()
                          << " dAU=" << dAU_n
                          << " gradp_face_momentum=(" << grad_face_momentum.x << ","
                          << grad_face_momentum.y << "," << grad_face_momentum.z << ")"
                          << " gradp_face_correction=(" << grad_face_correction.x << ","
                          << grad_face_correction.y << "," << grad_face_correction.z << ")"
                          << " delta_gradp="
                          << (grad_face_momentum-grad_face_correction).mag()
                          << " momentum_pressure_flux=" << momentum_pressure_flux
                          << " correction_pressure_flux=" << pressure_flux_correction
                          << " rhie_chow_expected=" << rhie_chow_expected
                          << " rhie_chow_actual=" << rhie_chow_actual
                          << " delta_rhie_chow="
                          << (rhie_chow_actual-rhie_chow_expected)
                          << "\n";
            }

            std::cerr << "local_mass_flux_divergence=" << local_div_authoritative
                      << " local_reconstructed_flux_divergence=" << local_div_reconstructed
                      << " reconstructed_velocity_continuity_linf="
                      << h.reconstructed_velocity_continuity_linf
                      << " delta_phi_linf=" << max_delta_phi
                      << " delta_u_face_linf=" << max_delta_u_face
                      << " delta_gradp_linf=" << max_delta_grad_p
                      << " momentum_pressure_flux_linf=" << max_momentum_pressure_flux
                      << " correction_pressure_flux_linf=" << max_correction_pressure_flux
                      << " rhie_chow_expected_linf=" << max_rhie_chow_expected
                      << " rhie_chow_actual_linf=" << max_rhie_chow_actual
                      << " delta_rhie_chow_linf=" << max_rhie_chow_mismatch
                      << "\n";
            std::cerr << "diagnostic_global=(" << h.momentum_equation_residual_components[0]
                      << "," << h.momentum_equation_residual_components[1]
                      << "," << h.momentum_equation_residual_components[2] << ")"
                      << " internal=" << h.momentum_equation_residual_internal
                      << " boundary=" << h.momentum_equation_residual_boundary
                      << " worst_cell=" << h.momentum_residual_cell
                      << " worst_patch=" << h.momentum_residual_patch << "\n";
            std::cerr << "no_pressure=" << h.momentum_residual_no_pressure
                      << " pressure_contribution=" << h.momentum_pressure_contribution
                      << " gradp_linf=" << h.pressure_gradient_linf
                      << "\n";
            std::cerr << "=== END CFDX MOMENTUM MICROSCOPE ===\n";
        }
        result.history.push_back(h);

        if (controls.probe_callback) {
            for (const auto& probe : controls.probes) {
                if (probe.name.empty())
                    throw std::invalid_argument("solve_steady_incompressible: probe name must not be empty");
                const double value = sample_incompressible_probe(probe, mesh, geometry, U, p);
                if (!std::isfinite(value))
                    throw std::runtime_error("solve_steady_incompressible: probe value is not finite");
                controls.probe_callback(IncompressibleProbeSample{probe.name, iter, 0.0, value});
            }
        }

        if (controls.iteration_output_callback &&
            !controls.iteration_output_callback(iter, 0.0, mesh, U, p)) {
            result.iterations = iter;
            break;
        }

        if (converged_now) {
            result.converged = true;
            result.iterations = iter;
            break;
        }
        result.iterations = iter;
    }

    return result;
}

} // namespace cfdx::physics