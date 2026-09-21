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

inline double energy_balance_relative(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& old_temperature,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& source,
    const EnergySolverControls& controls,
    const ScalarBoundaryConditions& bcs)
{
    double net_flux = 0.0;
    double source_total = 0.0;
    double accumulation = 0.0;
    const auto& own = mesh.ownership();

    for(std::size_t c=0;c<mesh.n_cells();++c) {
        source_total += source(c) * geometry.cell_volumes[c];
        if(controls.dt > 0.0)
            accumulation += controls.density * controls.cp *
                geometry.cell_volumes[c] *
                (temperature(c)-old_temperature(c))/controls.dt;
    }

    for(std::size_t f=0;f<mesh.n_faces();++f) {
        const std::size_t o=own.owner(f);
        const auto nraw=own.neighbour(f);
        if(nraw>=0) continue;

        const std::size_t patch=geometry.face_patch[f];
        ScalarBoundaryCondition bc;
        if(patch<mesh.boundary().n_patches()) {
            const auto& name=mesh.boundary().patch(patch).name;
            auto it=bcs.find(name);
            if(it!=bcs.end()) bc=it->second;
        }

        const auto& Sf0=geometry.face_area_vectors[f];
        const auto dvec=geometry.face_centres[f]-geometry.cell_centres[o];
        const double sign=(Sf0.dot(dvec)>=0.0)?1.0:-1.0;
        const auto Sf=Sf0*sign;
        const double area=Sf.mag();
        const double d=dvec.mag();
        const double F=mass_flux(f);
        double Tf=temperature(o);
        if(bc.type==ScalarBoundaryType::FIXED_VALUE)
            Tf=bc.value;
        else if(bc.type==ScalarBoundaryType::FIXED_GRADIENT)
            Tf=temperature(o)+bc.gradient*d;

        const double k=controls.conductivity;
        double conductive=-k*bc.gradient*area;
        if(bc.type==ScalarBoundaryType::FIXED_VALUE && d>0.0)
            conductive=-k*(Tf-temperature(o))/d*area;
        const double convective=F*(F>=0.0 ? temperature(o) : Tf);
        net_flux += convective + conductive;
    }

    const double imbalance=accumulation + net_flux - source_total;
    const double scale=std::max({1.0,std::abs(accumulation),
                                 std::abs(net_flux),std::abs(source_total)});
    return std::abs(imbalance)/scale;
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
        cfdx::core::Vector candidate(temperature.size(),0.0);
        for(std::size_t i=0;i<temperature.size();++i)
            candidate(i)=temperature(i);
        ScalarSolveControls sc;
        sc.max_iterations=2000;
        sc.tolerance=controls.tolerance;
        sc.relaxation=controls.relaxation;
        const auto linear=solve_scalar_equation(eq,candidate,sc);
        double res=scalar_equation_residual_inf(eq,candidate);

        for(std::size_t i=0;i<temperature.size();++i)
            temperature(i)=candidate(i);
        const double imbalance=energy_balance_relative(
            mesh,geometry,mass_flux,temperature,old,source,controls,bcs);
        result.history.push_back({iter,res,imbalance});
        result.iterations=iter;

        if(linear.status==cfdx::core::SolverStatus::CONVERGED &&
           res<=controls.tolerance) {
            result.converged=true;
            break;
        }
        // For a transient step, "old" is the state at t^n and must remain
        // fixed throughout the nonlinear iterations of the t^(n+1) solve.
        // For steady mode dt<=0 and old is only a linearisation reference.

    }
    return result;
}

} // namespace cfdx::physics
