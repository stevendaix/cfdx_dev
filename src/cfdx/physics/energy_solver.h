#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/finite_volume_transport.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct EnergySolverControls {
    double density = 1.0;
    double cp = 1000.0;
    double conductivity = 1.0;
    double dt = 0.0;                 // <=0 => steady
    double relaxation = 0.9;
    std::size_t max_iterations = 500;
    double tolerance = 1e-10;
};

struct EnergyIteration {
    std::size_t iteration = 0;
    double residual = 0.0;
    double energy_imbalance = 0.0;
};

struct EnergySolveResult {
    bool converged = false;
    std::size_t iterations = 0;
    std::vector<EnergyIteration> history;
};

inline void validate_energy_controls(const EnergySolverControls& c)
{
    if (c.density <= 0.0 || c.cp <= 0.0 || c.conductivity < 0.0 ||
        !(c.relaxation > 0.0 && c.relaxation <= 1.0) ||
        c.max_iterations == 0 || c.tolerance <= 0.0)
        throw std::invalid_argument("invalid energy solver controls");
    if (c.dt < 0.0) throw std::invalid_argument("energy time step must be >= 0");
}

inline ScalarEquation assemble_energy_equation(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& source,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& old_temperature,
    const EnergySolverControls& c,
    const ScalarBoundaryConditions& bcs = {},
    const ScalarBoundaryFaceValues* face_values = nullptr)
{
    validate_energy_controls(c);
    if (source.size()!=mesh.n_cells() || old_temperature.size()!=mesh.n_cells())
        throw std::invalid_argument("energy field size mismatch");

    cfdx::core::Field<double,cfdx::core::Location::CELL> su(
        mesh.n_cells(),"energy_source","W/m3",1);
    cfdx::core::Field<double,cfdx::core::Location::CELL> sp(
        mesh.n_cells(),"energy_sp","W/m3/K",1);
    std::vector<double> transient_diag(mesh.n_cells(),0.0);
    std::vector<double> transient_rhs(mesh.n_cells(),0.0);

    for(std::size_t i=0;i<mesh.n_cells();++i) {
        su(i)=source(i);
        sp(i)=0.0;
        if(c.dt>0.0) {
            transient_diag[i]=c.density*c.cp*geometry.cell_volumes[i]/c.dt;
            transient_rhs[i]=transient_diag[i]*old_temperature(i);
        }
    }

    return assemble_scalar_equation(
        mesh,geometry,mass_flux,c.conductivity,su,sp,bcs,true,
        face_values,&transient_diag,&transient_rhs);
}

inline EnergySolveResult solve_energy(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& source,
    const EnergySolverControls& controls = {},
    const ScalarBoundaryConditions& bcs = {},
    const ScalarBoundaryFaceValues* face_values = nullptr)
{
    validate_energy_controls(controls);
    if(temperature.size()!=mesh.n_cells() || source.size()!=mesh.n_cells())
        throw std::invalid_argument("energy field size mismatch");

    EnergySolveResult result;
    auto old=temperature;
    for(std::size_t iter=1;iter<=controls.max_iterations;++iter) {
        auto eq=assemble_energy_equation(
            mesh,geometry,mass_flux,source,old,controls,bcs,face_values);
        auto candidate=temperature;
        ScalarSolveControls sc;
        sc.max_iterations=2000;
        sc.tolerance=controls.tolerance;
        sc.relaxation=controls.relaxation;
        const auto linear=solve_scalar_equation(eq,candidate,sc);
        double res=scalar_equation_residual_inf(eq,candidate);

        double imbalance=0.0;
        for(std::size_t i=0;i<mesh.n_cells();++i)
            imbalance=std::max(imbalance,std::abs(eq.rhs(i)));

        temperature=candidate;
        result.history.push_back({iter,res,imbalance});
        result.iterations=iter;

        if(linear.status==cfdx::core::SolverStatus::CONVERGED &&
           res<=controls.tolerance) {
            result.converged=true;
            break;
        }
        old=temperature;
    }
    return result;
}

} // namespace cfdx::physics
