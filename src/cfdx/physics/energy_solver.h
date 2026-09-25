#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/finite_volume_transport.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <iostream>
#include <limits>
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
    // Legacy common tolerance. Dedicated values override it when > 0.
    double tolerance = 1e-10;
    double linear_tolerance = -1.0;
    double temperature_tolerance = -1.0;
    double energy_balance_tolerance = -1.0;
};

struct EnergyIteration {
    std::size_t iteration = 0;
    double residual = 0.0;
    double backward_error = 0.0;
    double energy_imbalance = 0.0;
};

struct EnergySolveResult {
    bool converged = false;
    std::size_t iterations = 0;
    std::vector<EnergyIteration> history;
};

inline void validate_energy_controls(const EnergySolverControls& c)
{
    if (!std::isfinite(c.density) || !std::isfinite(c.cp) ||
        !std::isfinite(c.conductivity) || !std::isfinite(c.dt) ||
        !std::isfinite(c.relaxation) || !std::isfinite(c.tolerance) ||
        !std::isfinite(c.linear_tolerance) ||
        !std::isfinite(c.temperature_tolerance) ||
        !std::isfinite(c.energy_balance_tolerance) ||
        c.density <= 0.0 || c.cp <= 0.0 || c.conductivity < 0.0 ||
        !(c.relaxation > 0.0 && c.relaxation <= 1.0) ||
        c.max_iterations == 0 || c.tolerance <= 0.0 ||
        (c.linear_tolerance <= 0.0 && c.linear_tolerance != -1.0) ||
        (c.temperature_tolerance <= 0.0 && c.temperature_tolerance != -1.0) ||
        (c.energy_balance_tolerance <= 0.0 && c.energy_balance_tolerance != -1.0))
        throw std::invalid_argument("invalid energy solver controls");
    if (c.dt < 0.0) throw std::invalid_argument("energy time step must be >= 0");
}

inline cfdx::core::Field<double,cfdx::core::Location::FACE> enthalpy_face_flux(
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    double cp)
{
    if (!(cp > 0.0) || !std::isfinite(cp))
        throw std::invalid_argument("enthalpy_face_flux: cp must be positive");
    auto flux=mass_flux;
    for(std::size_t f=0;f<flux.size();++f) flux(f)=mass_flux(f)*cp;
    return flux;
}

inline void validate_energy_implicit_source(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>* source_implicit,
    std::size_t n_cells)
{
    if (!source_implicit) return;
    if (source_implicit->size()!=n_cells)
        throw std::invalid_argument("energy implicit source field size mismatch");
    for(std::size_t i=0;i<n_cells;++i)
        if (!std::isfinite((*source_implicit)(i)) || (*source_implicit)(i)>0.0)
            throw std::invalid_argument(
                "energy implicit source must be finite and non-positive");
}

inline ScalarEquation assemble_energy_equation(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& source,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& old_temperature,
    const EnergySolverControls& c,
    const ScalarBoundaryConditions& bcs = {},
    const ScalarBoundaryFaceValues* face_values = nullptr,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>* source_implicit = nullptr)
{
    validate_energy_controls(c);
    if (source.size()!=mesh.n_cells() || old_temperature.size()!=mesh.n_cells())
        throw std::invalid_argument("energy field size mismatch");
    validate_energy_implicit_source(source_implicit,mesh.n_cells());

    cfdx::core::Field<double,cfdx::core::Location::CELL> su(
        mesh.n_cells(),"energy_source","W/m3",1);
    cfdx::core::Field<double,cfdx::core::Location::CELL> sp(
        mesh.n_cells(),"energy_sp","W/m3/K",1);
    std::vector<double> transient_diag(mesh.n_cells(),0.0);
    std::vector<double> transient_rhs(mesh.n_cells(),0.0);

    for(std::size_t i=0;i<mesh.n_cells();++i) {
        su(i)=source(i);
        sp(i)=source_implicit ? (*source_implicit)(i) : 0.0;
        if(c.dt>0.0) {
            transient_diag[i]=c.density*c.cp*geometry.cell_volumes[i]/c.dt;
            transient_rhs[i]=transient_diag[i]*old_temperature(i);
        }
    }

    const auto enthalpy_flux=enthalpy_face_flux(mass_flux,c.cp);
    return assemble_scalar_equation(
        mesh,geometry,enthalpy_flux,c.conductivity,su,sp,bcs,true,
        face_values,&transient_diag,&transient_rhs);
}

inline double energy_balance_relative_from_equation(
    const ScalarEquation& equation,
    const cfdx::core::Vector& temperature)
{
    if (temperature.size() != equation.rhs.size())
        throw std::invalid_argument("energy balance field size mismatch");
    double imbalance = 0.0;
    double scale = 1.0;
    const auto* row = equation.matrix.row_offsets_data();
    const auto* col = equation.matrix.columns_data();
    const auto* val = equation.matrix.values_data();
    for (std::size_t i = 0; i < temperature.size(); ++i) {
        double ri = -equation.rhs(i);
        double row_scale = std::abs(equation.rhs(i));
        for (std::uint32_t k = row[i]; k < row[i + 1]; ++k) {
            const double term = val[k] * temperature(col[k]);
            ri += term;
            row_scale += std::abs(term);
        }
        imbalance += ri;
        scale = std::max(scale, row_scale);
    }
    return std::abs(imbalance) / scale;
}

inline double energy_balance_relative(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& old_temperature,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& source,
    const EnergySolverControls& controls,
    const ScalarBoundaryConditions& bcs,
    const ScalarBoundaryFaceValues* face_values = nullptr,
    const std::vector<double>* conductivity_values = nullptr)
{
    if (conductivity_values && conductivity_values->size()!=mesh.n_cells())
        throw std::invalid_argument("energy conductivity field size mismatch");
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
        bool face_override_valid=false;
        if(face_values && patch<mesh.boundary().n_patches()) {
            const auto it=face_values->values.find(mesh.boundary().patch(patch).name);
            if(it!=face_values->values.end() && f<it->second.size() &&
               std::isfinite(it->second[f])) {
                Tf=it->second[f];
                face_override_valid=true;
            }
        }
        if(!face_override_valid) {
            if(bc.type==ScalarBoundaryType::FIXED_VALUE)
                Tf=bc.value;
            else if(bc.type==ScalarBoundaryType::FIXED_GRADIENT)
                Tf=temperature(o)+bc.gradient*d;
        }

        const double k=conductivity_values
            ? (*conductivity_values)[o] : controls.conductivity;
        double conductive=-k*bc.gradient*area;
        if(d>0.0 && (bc.type==ScalarBoundaryType::FIXED_VALUE || face_override_valid))
            conductive=-k*(Tf-temperature(o))/d*area;
        const double convective=controls.cp*F*(F>=0.0 ? temperature(o) : Tf);
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
    const ScalarBoundaryFaceValues* face_values = nullptr,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>* source_implicit = nullptr,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>* previous_time_temperature = nullptr)
{
    validate_energy_controls(controls);
    if(temperature.size()!=mesh.n_cells() || source.size()!=mesh.n_cells())
        throw std::invalid_argument("energy field size mismatch");

    validate_energy_implicit_source(source_implicit,mesh.n_cells());
    if (previous_time_temperature &&
        previous_time_temperature->size()!=mesh.n_cells())
        throw std::invalid_argument("energy previous-time field size mismatch");

    EnergySolveResult result;
    // Previous physical time level stays fixed; T_iterate is the accepted
    // nonlinear/Picard state. For a coupled outer iteration, the caller can
    // explicitly provide the physical T^n state.
    auto T_previous_time = previous_time_temperature
        ? *previous_time_temperature : temperature;
    auto T_iterate=temperature;
    const double linear_tolerance =
        controls.linear_tolerance > 0.0 ? controls.linear_tolerance : controls.tolerance;
    const double temperature_tolerance =
        controls.temperature_tolerance > 0.0 ? controls.temperature_tolerance : controls.tolerance;
    const double energy_balance_tolerance =
        controls.energy_balance_tolerance > 0.0 ? controls.energy_balance_tolerance : controls.tolerance;
    for(std::size_t iter=1;iter<=controls.max_iterations;++iter) {
        auto eq=assemble_energy_equation(
            mesh,geometry,mass_flux,source,T_previous_time,controls,bcs,face_values,source_implicit);
        cfdx::core::Vector candidate(temperature.size(),0.0);
        for(std::size_t i=0;i<temperature.size();++i)
            candidate(i)=temperature(i);
        ScalarSolveControls sc;
        sc.max_iterations=2000;
        // Leave margin between the Krylov stopping criterion and the
        // componentwise backward-error gate below. This is deliberately
        // local to energy solves: ScalarSolveControls::tolerance remains a
        // relative Krylov tolerance for the shared momentum/pressure paths.
        sc.tolerance=temperature.size()==1
            ? linear_tolerance
            : std::max(
                10.0*std::numeric_limits<double>::epsilon(),
                0.001*linear_tolerance);
        // solve_scalar_equation returns the fully solved linear predictor;
        // nonlinear temperature relaxation is applied exactly once below.
        sc.relaxation=1.0;
        const auto linear=solve_scalar_equation(eq,candidate,sc);
        cfdx::core::Vector accepted(candidate.size(),0.0);
        for(std::size_t i=0;i<temperature.size();++i)
            accepted(i)=T_iterate(i)+controls.relaxation*(candidate(i)-T_iterate(i));

        const double res=scalar_equation_residual_inf(eq,accepted);

        double temperature_change=0.0;
        double temperature_scale=1.0;
        for(std::size_t i=0;i<temperature.size();++i) {
            temperature_change=std::max(
                temperature_change,
                std::abs(accepted(i)-T_iterate(i)));
            temperature_scale=std::max(
                temperature_scale,
                std::abs(accepted(i)));
        }
        const double relative_temperature_change=
            temperature_change/temperature_scale;

        for(std::size_t i=0;i<temperature.size();++i) {
            temperature(i)=accepted(i);
            T_iterate(i)=accepted(i);
        }

        const double imbalance=energy_balance_relative_from_equation(eq,accepted);
        result.iterations=iter;

        // Validate the *accepted* state, not the unrelaxed predictor.
        // The backward error is ||Ax-b||_inf divided by a componentwise
        // row scale ||A_i,: x||_1 + |b_i|. This must use the same state as
        // the absolute residual; otherwise under-relaxation can report
        // convergence for a state that does not actually solve the equation.
        double linear_backward_error = 0.0;
        const auto* row = eq.matrix.row_offsets_data();
        const auto* col = eq.matrix.columns_data();
        const auto* val = eq.matrix.values_data();
        for (std::size_t i = 0; i < accepted.size(); ++i) {
            double row_residual = -eq.rhs(i);
            double row_scale = std::abs(eq.rhs(i));
            for (std::uint32_t k = row[i]; k < row[i + 1]; ++k) {
                const double term = val[k] * accepted(col[k]);
                row_residual += term;
                row_scale += std::abs(term);
            }
            linear_backward_error = std::max(
                linear_backward_error,
                std::abs(row_residual) / std::max(1.0, row_scale));
        }
        const bool linear_converged =
            linear.status == cfdx::core::SolverStatus::CONVERGED &&
            std::isfinite(res) && std::isfinite(linear_backward_error) &&
            std::isfinite(linear.residual_relative) &&
            linear.residual_relative<=linear_tolerance &&
            linear_backward_error<=linear_tolerance;

        // The energy equation assembled here is linear for a fixed source,
        // conductivity and transient reference state. Therefore a fully
        // converged linear solve is already the exact fixed point when
        // relaxation=1. With under-relaxation, the accepted state is only a
        // predictor and must also stop moving before we declare convergence.
        const bool accepted_state_converged =
            controls.relaxation >= 1.0 - 10.0*std::numeric_limits<double>::epsilon() ||
            relative_temperature_change <= temperature_tolerance;
        result.history.push_back({iter,res,linear_backward_error,imbalance});
        if (mesh.n_cells() <= 64) {
            std::cerr << "THERMAL_RESIDUAL: iteration=" << iter
                      << " linear_status=" << static_cast<int>(linear.status)
                      << " residual=" << res
                      << " backward_error=" << linear_backward_error
                      << " linear_relative=" << linear.residual_relative
                      << " dT_relative=" << relative_temperature_change
                      << " energy_balance=" << imbalance
                      << " accepted_state=" << accepted_state_converged
                      << '\n';
        }
        if(linear_converged &&
           accepted_state_converged &&
           imbalance<=energy_balance_tolerance) {
            result.converged=true;
            break;
        }
        // T_previous_time remains fixed for transient physics. T_iterate is
        // updated above and is the state used for nonlinear/Picard relaxation.
    }
    return result;
}

} // namespace cfdx::physics
