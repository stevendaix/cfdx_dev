#include "cfdx/physics/continuation.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

Mesh make_cavity_mesh(std::size_t nx, std::size_t ny)
{
    Mesh mesh;
    mesh.points().resize((nx + 1) * (ny + 1) * 2);
    const auto id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };
    for (std::size_t j = 0; j <= ny; ++j)
        for (std::size_t i = 0; i <= nx; ++i) {
            const double x = static_cast<double>(i) / static_cast<double>(nx);
            const double y = static_cast<double>(j) / static_cast<double>(ny);
            mesh.points().set(id(i,j,0), x, y, 0.0);
            mesh.points().set(id(i,j,1), x, y, 1.0);
        }

    std::map<std::vector<std::size_t>, std::size_t> faces;
    std::vector<std::vector<std::size_t>> cell_faces(nx * ny);
    auto add_face = [&](std::initializer_list<std::size_t> vertices, std::size_t cell) {
        std::vector<std::size_t> key(vertices);
        std::sort(key.begin(), key.end());
        const auto found = faces.find(key);
        if (found != faces.end()) {
            mesh.ownership().set_neighbour(found->second, static_cast<int>(cell));
            return found->second;
        }
        const std::size_t face = mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        mesh.ownership().set_owner(face, cell);
        mesh.ownership().set_neighbour(face, FaceOwnership::BOUNDARY);
        faces.emplace(std::move(key), face);
        return face;
    };

    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t c = j * nx + i;
            const auto a=id(i,j,0), b=id(i+1,j,0), cc=id(i+1,j+1,0), d=id(i,j+1,0);
            const auto e=id(i,j,1), f=id(i+1,j,1), g=id(i+1,j+1,1), h=id(i,j+1,1);
            cell_faces[c] = {
                add_face({a,d,cc,b},c), add_face({e,f,g,h},c),
                add_face({a,b,f,e},c), add_face({d,h,g,cc},c),
                add_face({a,e,h,d},c), add_face({b,cc,g,f},c)};
        }
    for (const auto& cf : cell_faces)
        mesh.cells().push_cell(cf);

    Patch bottom{"bottom", PatchType::WALL, {}};
    Patch top{"top", PatchType::WALL, {}};
    Patch left{"left", PatchType::WALL, {}};
    Patch right{"right", PatchType::WALL, {}};
    Patch front{"front", PatchType::WALL, {}};
    Patch back{"back", PatchType::WALL, {}};

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) >= 0)
            continue;
        const auto& vertices = mesh.faces().vertices();
        const auto begin = vertices.begin() + static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto end = begin + static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));
        double x=0.0, y=0.0, z=0.0;
        for (auto p=begin; p!=end; ++p) {
            x += mesh.points().x(*p); y += mesh.points().y(*p); z += mesh.points().z(*p);
        }
        const double n = static_cast<double>(mesh.faces().face_size(f));
        x /= n; y /= n; z /= n;
        constexpr double tol = 1e-12;
        if (std::abs(y) < tol) bottom.face_ids.push_back(f);
        else if (std::abs(y-1.0) < tol) top.face_ids.push_back(f);
        else if (std::abs(x) < tol) left.face_ids.push_back(f);
        else if (std::abs(x-1.0) < tol) right.face_ids.push_back(f);
        else if (std::abs(z) < tol) front.face_ids.push_back(f);
        else if (std::abs(z-1.0) < tol) back.face_ids.push_back(f);
        else throw std::runtime_error("N7 cavity: unclassified boundary face");
    }
    mesh.boundary().add_patch(bottom); mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(left); mesh.boundary().add_patch(right);
    mesh.boundary().add_patch(front); mesh.boundary().add_patch(back);
    return mesh;
}

struct Solution {
    Field<double,Location::CELL> U_direct;
    Field<double,Location::CELL> U_continuation;
    IncompressibleSolveResult direct;
    ContinuationSolveResult continuation;
    FvGeometry geometry;
};

Solution solve_n7_case()
{
    constexpr std::size_t n = 32;
    Mesh mesh = make_cavity_mesh(n,n);
    Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    Field<double,Location::CELL> p(mesh.n_cells(),"p","Pa",1);
    U.fill(0.0);
    p.fill(0.0);

    VelocityBoundaryConditions ubc;
    ubc["bottom"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["left"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["right"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["front"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["back"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["top"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};

    ScalarBoundaryConditions pbc;
    for (const char* name : {"bottom","top","left","right","front","back"})
        pbc[name]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

    IncompressibleSolverControls controls;
    controls.algorithm=PressureVelocityAlgorithm::SIMPLE;
    controls.density=1.0;
    controls.kinematic_viscosity=0.01; // Re=100, U=1, L=1.
    controls.linear_max_iterations=10000;
    controls.linear_tolerance=1e-8;
    controls.pressure_reference_cell=(n/2)*n+(n/2);
    controls.pressure_reference_value=0.0;
    controls.use_bounded_convection=true;
    controls.convection_scheme=ConvectionScheme::UPWIND;
    controls.coupling.alpha_u=0.7;
    controls.coupling.alpha_p=0.3;
    controls.convergence.max_iterations=2500;
    controls.convergence.relative_tolerance=1e-8;
    controls.convergence.continuity_tolerance=1e-8;

    const auto direct = solve_steady_incompressible(mesh,U,p,ubc,pbc,controls);
    if (!direct.converged)
        throw std::runtime_error("N7 direct reference solve did not converge");

    Field<double,Location::CELL> U_cont(mesh.n_cells(),"U","m/s",3);
    Field<double,Location::CELL> p_cont(mesh.n_cells(),"p","Pa",1);
    U_cont.fill(0.0);
    p_cont.fill(0.0);

    controls.adaptive_relaxation.enabled=true;
    controls.adaptive_relaxation.min_alpha_u=0.2;
    controls.adaptive_relaxation.max_alpha_u=0.9;
    controls.adaptive_relaxation.min_alpha_p=0.1;
    controls.adaptive_relaxation.max_alpha_p=0.5;

    ContinuationControls continuation;
    continuation.enabled=true;
    continuation.initial_step=0.25;
    continuation.minimum_step=0.125;
    continuation.maximum_step=0.5;
    continuation.step_growth=1.5;
    continuation.step_reduction=0.5;
    continuation.max_stage_attempts=8;
    continuation.max_stages=16;

    const auto continuation_result =
        solve_steady_incompressible_continuation(
            mesh,U_cont,p_cont,ubc,pbc,controls,continuation);
    if (!continuation_result.converged) {
        std::ostringstream message;
        message << "N7 continuation solve did not converge: "
                << continuation_result.reason
                << " final_parameter=" << continuation_result.final_parameter
                << " attempts=" << continuation_result.total_attempts;
        if (!continuation_result.stages.empty()) {
            const auto& stage = continuation_result.stages.back();
            message << " last_stage=" << stage.stage
                    << " target=" << stage.parameter
                    << " step=" << stage.step
                    << " attempts=" << stage.attempts
                    << " status=" << cfdx::core::to_string(stage.status)
                    << " reason=" << stage.reason;
        }
        throw std::runtime_error(message.str());
    }

    return {std::move(U),std::move(U_cont),direct,std::move(continuation_result),build_fv_geometry(mesh)};
}

double kinetic_energy(const Field<double,Location::CELL>& U,const FvGeometry& geometry)
{
    double energy=0.0;
    const double* ux=U.component_data(0);
    const double* uy=U.component_data(1);
    const double* uz=U.component_data(2);
    for (std::size_t c=0;c<geometry.cell_volumes.size();++c)
        energy += 0.5*(ux[c]*ux[c]+uy[c]*uy[c]+uz[c]*uz[c])*geometry.cell_volumes[c];
    return energy;
}

double centre_component(const Field<double,Location::CELL>& U,bool x_component)
{
    constexpr std::size_t n=32;
    const std::size_t i=n/2-1;
    const std::size_t j=n/2-1;
    return U.component_data(x_component ? 0 : 1)[j*n+i];
}

} // namespace

int main()
{
    try {
        const auto result = solve_n7_case();
        const auto& final = result.continuation.final_solver_result;
        const auto& history = final.history;
        if (history.empty())
            throw std::runtime_error("N7 continuation produced no iteration history");

        const double continuity = history.back().continuity_linf;
        const double momentum = history.back().momentum_equation_residual;
        if (continuity > 1e-7 || momentum > 1e-7)
            throw std::runtime_error("N7 final physical convergence gates failed");

        if (result.continuation.final_parameter != 1.0 ||
            result.continuation.stages.size() < 2 ||
            result.continuation.total_attempts < result.continuation.stages.size())
            throw std::runtime_error("N7 continuation stage contract failed");

        for (const auto& stage : result.continuation.stages) {
            if (!stage.converged ||
                stage.status != cfdx::core::ConvergenceStatus::CONVERGED ||
                !std::isfinite(stage.parameter) ||
                !std::isfinite(stage.step))
                throw std::runtime_error("N7 continuation stage is not qualified");
        }

        const double direct_energy=kinetic_energy(result.U_direct,result.geometry);
        const double continuation_energy=kinetic_energy(result.U_continuation,result.geometry);
        if (!std::isfinite(direct_energy) || !std::isfinite(continuation_energy))
            throw std::runtime_error("N7 independent QoI is non-finite");

        const double direct_u=centre_component(result.U_direct,true);
        const double direct_v=centre_component(result.U_direct,false);
        const double continuation_u=centre_component(result.U_continuation,true);
        const double continuation_v=centre_component(result.U_continuation,false);

        // Independent QoI: the converged Re=100 cavity centre velocity remains
        // close to the published Ghia reference while continuation/adaptive
        // controls reproduce the direct nonlinear solution.
        if (std::abs(continuation_u - direct_u) > 2e-5 ||
            std::abs(continuation_v - direct_v) > 2e-5)
            throw std::runtime_error("N7 continuation changed the direct-solve QoI");

        constexpr double ghia_u=-0.20581;
        constexpr double ghia_v=0.05454;
        if (std::abs(continuation_u-ghia_u) > 0.15 ||
            std::abs(continuation_v-ghia_v) > 0.15)
            throw std::runtime_error("N7 Ghia centre-velocity QoI gate failed");

        std::size_t adaptive_changes=0;
        std::size_t adaptive_history_samples=0;
        for (const auto& stage : result.continuation.stages) {
            for (const auto& h : stage.solver_result.history) {
                ++adaptive_history_samples;
                if (std::abs(h.effective_alpha_u-0.7) > 1e-12 ||
                    std::abs(h.effective_alpha_p-0.3) > 1e-12)
                    ++adaptive_changes;
                if (h.effective_alpha_u < 0.2 || h.effective_alpha_u > 0.9 ||
                    h.effective_alpha_p < 0.1 || h.effective_alpha_p > 0.5 ||
                    !std::isfinite(h.nonlinear_convergence_metric))
                    throw std::runtime_error("N7 adaptive relaxation bounds/history failed");
            }
        }
        if (adaptive_history_samples < 3 || adaptive_changes == 0)
            throw std::runtime_error("N7 adaptive relaxation did not produce qualified history");

        if (std::abs(continuation_energy-direct_energy) > 2e-5)
            throw std::runtime_error("N7 kinetic-energy QoI mismatch");

        return 0;
    } catch (const std::exception& e) {
        return std::fprintf(stderr,"N7_CONVERGENCE_QUALIFICATION: FAIL: %s\n",e.what()), 1;
    }
}
