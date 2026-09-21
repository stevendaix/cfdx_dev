#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "cfdx/physics/solver_control.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct VelocityBoundaryCondition {
    enum class Type { FIXED_VALUE, ZERO_GRADIENT };
    Type type = Type::ZERO_GRADIENT;
    cfdx::core::Vec3 value{0.0, 0.0, 0.0};
};

using VelocityBoundaryConditions = std::map<std::string, VelocityBoundaryCondition>;

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
};

struct IncompressibleIteration {
    std::size_t iteration = 0;
    double momentum_residual = std::numeric_limits<double>::infinity();
    double pressure_residual = std::numeric_limits<double>::infinity();
    double continuity_l1 = std::numeric_limits<double>::infinity();
    double continuity_linf = std::numeric_limits<double>::infinity();
    double velocity_change_inf = std::numeric_limits<double>::infinity();
    double pressure_change_inf = std::numeric_limits<double>::infinity();
    double momentum_equation_residual = std::numeric_limits<double>::infinity();
    double continuity_normalized = std::numeric_limits<double>::infinity();
    std::size_t momentum_linear_iterations = 0;
    std::size_t pressure_linear_iterations = 0;
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
    const VelocityBoundaryConditions& bcs)
{
    using namespace cfdx::core;
    Field<double, Location::FACE> flux(mesh.n_faces(), "phi", "kg/s", 1);
    const auto& own = mesh.ownership();

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
    using namespace cfdx::core;
    if(p.size()!=mesh.n_cells() || rAU.size()!=mesh.n_cells())
        throw std::invalid_argument("make_rhie_chow_mass_flux: field size mismatch");
    auto flux=make_mass_flux(mesh,geometry,U,rho,bcs);
    auto gradp=gauss_gradient_with_boundary(p,mesh,geometry,{});
    for(std::size_t f=0;f<mesh.n_faces();++f) {
        const auto nr=mesh.ownership().neighbour(f);
        if(nr<0) continue;
        const std::size_t o=mesh.ownership().owner(f);
        const std::size_t n=static_cast<std::size_t>(nr);
        const double d=(geometry.cell_centres[n]-geometry.cell_centres[o]).mag();
        if(d<=0.0) throw std::runtime_error("make_rhie_chow_mass_flux: degenerate face");
        const double rface=0.5*(rAU[o]+rAU[n]);
        const double dpdn=(p(n)-p(o))/d;
        const cfdx::core::Vec3 gpface={
            0.5*(gradp.component_data(0)[o]+gradp.component_data(0)[n]),
            0.5*(gradp.component_data(1)[o]+gradp.component_data(1)[n]),
            0.5*(gradp.component_data(2)[o]+gradp.component_data(2)[n])};
        const cfdx::core::Vec3 Sf=geometry.face_area_vectors[f];
        const double gradface=gpface.dot(Sf);
        const double orth=dpdn*Sf.mag();
        flux(f)-=rho*rface*(orth-gradface);
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
    bool bounded)
{
    cfdx::core::Field<double, cfdx::core::Location::CELL> source(
        mesh.n_cells(), "momentum_source", "N/m3", 1);
    cfdx::core::Field<double, cfdx::core::Location::CELL> sp(
        mesh.n_cells(), "momentum_sp", "kg/m3/s", 1);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        source(c) = body_component(c) - pressure_gradient_component(c);
        sp(c) = 0.0;
    }
    (void)component;
    return assemble_scalar_equation(
        mesh, geometry, mass_flux, effective_dynamic_viscosity,
        source, sp, velocity_bcs, bounded);
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

inline IncompressibleSolveResult solve_steady_incompressible(
    const cfdx::core::Mesh& mesh,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    const IncompressibleSolverControls& controls = {})
{
    using namespace cfdx::core;
    if (U.dimension() != 3 || U.size() != mesh.n_cells() ||
        p.dimension() != 1 || p.size() != mesh.n_cells())
        throw std::invalid_argument("solve_steady_incompressible: invalid fields");

    validate_incompressible_controls(controls, mesh.n_cells());
    const FvGeometry geometry = build_fv_geometry(mesh);
    const double mu_eff = controls.density *
        (controls.kinematic_viscosity + controls.turbulent_viscosity);

    IncompressibleSolveResult result;
    // The steady solver's outer loop is the nonlinear convergence loop. The
    // algorithm controls how many pressure corrections are performed inside
    // each nonlinear iteration, mirroring the SIMPLE/PISO/PIMPLE structure.
    const std::size_t pcorr =
        controls.algorithm == PressureVelocityAlgorithm::SIMPLE
            ? 1u : static_cast<std::size_t>(controls.coupling.n_pressure_correctors);

    for (std::size_t iter = 1; iter <= controls.convergence.max_iterations; ++iter) {
        const auto U_old = U;
        const auto p_old = p;
        auto mass_flux = make_mass_flux(mesh, geometry, U, controls.density, velocity_bcs);
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

        auto ex = assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body_x, mu_eff, ubc_x, 0,
            controls.use_bounded_convection);
        auto ey = assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body_y, mu_eff, ubc_y, 1,
            controls.use_bounded_convection);
        auto ez = assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body_z, mu_eff, ubc_z, 2,
            controls.use_bounded_convection);

        Vector ux(mesh.n_cells()), uy(mesh.n_cells()), uz(mesh.n_cells());
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            ux(c) = U.component_data(0)[c];
            uy(c) = U.component_data(1)[c];
            uz(c) = U.component_data(2)[c];
        }

        const auto rx = solve_scalar_equation(ex, ux, {controls.linear_max_iterations,
                                                         controls.linear_tolerance,
                                                         controls.coupling.alpha_u});
        const auto ry = solve_scalar_equation(ey, uy, {controls.linear_max_iterations,
                                                         controls.linear_tolerance,
                                                         controls.coupling.alpha_u});
        const auto rz = solve_scalar_equation(ez, uz, {controls.linear_max_iterations,
                                                         controls.linear_tolerance,
                                                         controls.coupling.alpha_u});

        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            U.component_data(0)[c] = ux(c);
            U.component_data(1)[c] = uy(c);
            U.component_data(2)[c] = uz(c);
        }

        std::vector<double> rAU(mesh.n_cells());
        for (std::size_t c=0;c<mesh.n_cells();++c) {
            const double ax=std::max(ex.diagonal[c],1e-30);
            const double ay=std::max(ey.diagonal[c],1e-30);
            const double az=std::max(ez.diagonal[c],1e-30);
            rAU[c]=geometry.cell_volumes[c]/
                std::max((ax+ay+az)/3.0,1e-30);
        }
        mass_flux=make_rhie_chow_mass_flux(
            mesh,geometry,U,p,rAU,controls.density,velocity_bcs);

        double pressure_residual = std::numeric_limits<double>::infinity();
        std::size_t pressure_iterations = 0;

        for (std::size_t corr = 0; corr < pcorr; ++corr) {
            // Assemble the pressure-correction Laplacian with the same
            // face coefficient used for the velocity correction.
            const std::size_t nc = mesh.n_cells();
            SparseMatrix A(nc, nc);
            Vector b(nc, 0.0);
            std::vector<std::map<std::size_t, double>> rows(nc);
            std::vector<double> diag(nc, 0.0);

            std::vector<double> continuity(nc, 0.0);
            const auto* cell_faces = mesh.cells().faces_data();
            const auto* cell_offsets = mesh.cells().offsets_data();
            for (std::size_t c = 0; c < nc; ++c) {
                const Offset off = cell_offsets[c];
                const Offset count = cell_offsets[c + 1] - off;
                for (Offset k = 0; k < count; ++k) {
                    const std::size_t f = cell_faces[off + k];
                    const double sf = mass_flux(f);
                    continuity[c] += mesh.ownership().owner(f) == c ? sf : -sf;
                }
            }

            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                const auto nraw = mesh.ownership().neighbour(f);
                if (nraw < 0) continue;
                const std::size_t o = mesh.ownership().owner(f);
                const std::size_t n = static_cast<std::size_t>(nraw);
                const double d = (geometry.cell_centres[n] - geometry.cell_centres[o]).mag();
                const double area = geometry.face_area_vectors[f].mag();
                const double dface = 0.5 * (rAU[o] + rAU[n]) * area / d;
                const double coeff = controls.density * dface;
                diag[o] += coeff;
                diag[n] += coeff;
                rows[o][n] -= coeff;
                rows[n][o] -= coeff;
            }

            for (std::size_t c = 0; c < nc; ++c) {
                rows[c][c] += diag[c];
                b(c) = -continuity[c];
            }

            // Pressure has a gauge freedom whenever only Neumann-type
            // conditions are present. A single reference cell removes it.
            const std::size_t ref = controls.pressure_reference_cell;
            rows[ref].clear();
            rows[ref][ref] = 1.0;
            b(ref) = 0.0;

            for (std::size_t row = 0; row < nc; ++row) {
                for (const auto& [col, value] : rows[row]) A.push_back(row, col, value);
            }
            A.finalize();

            Vector p_corr(nc, 0.0);
            const auto rp = solve_bicgstab(
                A, b, p_corr, controls.linear_max_iterations, controls.linear_tolerance);
            pressure_residual = rp.residual_relative;
            pressure_iterations = rp.iterations;

            for (std::size_t c = 0; c < nc; ++c)
                p(c) += controls.coupling.alpha_p * p_corr(c);
            // Explicitly enforce the selected pressure gauge after relaxation.
            p(controls.pressure_reference_cell) = controls.pressure_reference_value;

            Field<double, Location::CELL> p_corr_field(
                nc, "p_corr", "Pa", 1);
            for (std::size_t c = 0; c < nc; ++c) p_corr_field(c) = p_corr(c);
            ScalarBoundaryConditions pressure_correction_bcs;
            for (const auto& [name, bc] : pressure_bcs) {
                pressure_correction_bcs[name] = bc;
                if (bc.type == ScalarBoundaryType::FIXED_VALUE)
                    pressure_correction_bcs[name].value = 0.0;
            }
            auto grad_pc = gauss_gradient_with_boundary(
                p_corr_field, mesh, geometry, pressure_correction_bcs);
            for (std::size_t c = 0; c < nc; ++c) {
                U.component_data(0)[c] -= rAU[c] * grad_pc.component_data(0)[c];
                U.component_data(1)[c] -= rAU[c] * grad_pc.component_data(1)[c];
                U.component_data(2)[c] -= rAU[c] * grad_pc.component_data(2)[c];
            }

            // Consistent face-flux correction. Fixed-value velocity patches
            // are left constrained by their prescribed face velocity.
            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                const auto nraw = mesh.ownership().neighbour(f);
                if (nraw < 0) continue;
                const std::size_t o = mesh.ownership().owner(f);
                const std::size_t n = static_cast<std::size_t>(nraw);
                const double d = (geometry.cell_centres[n] - geometry.cell_centres[o]).mag();
                const double area = geometry.face_area_vectors[f].mag();
                const double dface = 0.5 * (rAU[o] + rAU[n]) * area / d;
                mass_flux(f) -= controls.density * dface *
                                (p_corr(n) - p_corr(o));
            }
        }

        double l1 = 0.0;
        double linf = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            double div = 0.0;
            const Offset off = mesh.cells().offsets_data()[c];
            const Offset count = mesh.cells().offsets_data()[c + 1] - off;
            for (Offset k = 0; k < count; ++k) {
                const std::size_t f = mesh.cells().faces_data()[off + k];
                div += mesh.ownership().owner(f) == c ? mass_flux(f) : -mass_flux(f);
            }
            l1 += std::abs(div);
            linf = std::max(linf, std::abs(div));
        }

        // Reassemble the converged discrete momentum equations using the
        // corrected U/p fields. The linear-solver residuals above describe the
        // intermediate momentum solves before pressure correction; they are not
        // a physical residual of the final SIMPLE state.
        auto final_mass_flux = make_mass_flux(
            mesh, geometry, U, controls.density, velocity_bcs);
        auto final_grad_p = gauss_gradient_with_boundary(
            p, mesh, geometry, pressure_bcs);
        auto final_ex = assemble_momentum_component(
            mesh, geometry, final_mass_flux, final_grad_p, body_x, mu_eff,
            ubc_x, 0, controls.use_bounded_convection);
        auto final_ey = assemble_momentum_component(
            mesh, geometry, final_mass_flux, final_grad_p, body_y, mu_eff,
            ubc_y, 1, controls.use_bounded_convection);
        auto final_ez = assemble_momentum_component(
            mesh, geometry, final_mass_flux, final_grad_p, body_z, mu_eff,
            ubc_z, 2, controls.use_bounded_convection);
        Vector final_ux(mesh.n_cells()), final_uy(mesh.n_cells()), final_uz(mesh.n_cells());
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            final_ux(c) = U.component_data(0)[c];
            final_uy(c) = U.component_data(1)[c];
            final_uz(c) = U.component_data(2)[c];
        }
        const double final_momentum_residual = std::max({
            scalar_equation_residual_inf(final_ex, final_ux),
            scalar_equation_residual_inf(final_ey, final_uy),
            scalar_equation_residual_inf(final_ez, final_uz)});
        double velocity_scale = 1.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            velocity_scale = std::max(velocity_scale, std::sqrt(
                std::pow(U.component_data(0)[c], 2) +
                std::pow(U.component_data(1)[c], 2) +
                std::pow(U.component_data(2)[c], 2)));
        const double domain_volume = std::accumulate(
            geometry.cell_volumes.begin(), geometry.cell_volumes.end(), 0.0);
        const double characteristic_area = std::max(
            std::pow(domain_volume, 2.0 / 3.0), 1e-30);
        const double continuity_normalized = linf /
            std::max(controls.density * velocity_scale * characteristic_area, 1e-30);

        double velocity_change_inf = 0.0;
        double pressure_change_inf = 0.0;
        double velocity_scale = 1.0;
        double pressure_scale = 1.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            for (std::size_t d = 0; d < 3; ++d) {
                velocity_change_inf = std::max(
                    velocity_change_inf,
                    std::abs(U.component_data(d)[c] - U_old.component_data(d)[c]));
                velocity_scale = std::max(
                    velocity_scale, std::abs(U.component_data(d)[c]));
            }
            pressure_change_inf = std::max(
                pressure_change_inf, std::abs(p(c) - p_old(c)));
            pressure_scale = std::max(pressure_scale, std::abs(p(c)));
        }
        velocity_change_inf /= velocity_scale;
        pressure_change_inf /= pressure_scale;

        IncompressibleIteration h;
        h.iteration = iter;
        h.momentum_residual = std::max({rx.residual_relative, ry.residual_relative,
                                        rz.residual_relative});
        h.pressure_residual = pressure_residual;
        h.continuity_l1 = l1;
        h.continuity_linf = linf;
        h.velocity_change_inf = velocity_change_inf;
        h.pressure_change_inf = pressure_change_inf;
        h.momentum_equation_residual = final_momentum_residual;
        h.continuity_normalized = continuity_normalized;
        h.momentum_linear_iterations = std::max({rx.iterations, ry.iterations, rz.iterations});
        h.pressure_linear_iterations = pressure_iterations;
        result.history.push_back(h);

        if (iter > 1 &&
            std::isfinite(h.momentum_residual) && std::isfinite(h.pressure_residual) &&
            h.momentum_residual <= controls.convergence.relative_tolerance &&
            h.pressure_residual <= controls.convergence.relative_tolerance &&
            h.momentum_equation_residual <= controls.convergence.relative_tolerance &&
            h.continuity_linf <= controls.convergence.continuity_tolerance &&
            h.velocity_change_inf <= controls.convergence.relative_tolerance &&
            h.pressure_change_inf <= controls.convergence.relative_tolerance) {
            result.converged = true;
            result.iterations = iter;
            break;
        }
        result.iterations = iter;
    }

    return result;
}

} // namespace cfdx::physics
