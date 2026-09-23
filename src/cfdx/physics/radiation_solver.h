#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/radiation.h"
#include "cfdx/physics/energy_solver.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <limits>

namespace cfdx::physics {

struct RadiationTransportControls {
    double absorption = 0.0;
    double scattering = 0.0;
    double intensity_relaxation = 1.0;
    std::size_t max_iterations = 100;
    double tolerance = 1e-8;
    std::size_t linear_max_iterations = 2000;
    double linear_tolerance = 1e-10;
};

struct RadiationIteration {
    std::size_t iteration = 0;
    double intensity_residual = 0.0;
    double energy_source_residual = 0.0;
};

struct RadiationSolveResult {
    bool converged = false;
    std::size_t iterations = 0;
    std::vector<RadiationIteration> history;
};

inline void validate_radiation_transport_controls(const RadiationTransportControls& c)
{
    if(c.absorption<0.0 || c.scattering<0.0 ||
       c.max_iterations==0 || c.tolerance<=0.0 ||
       c.linear_max_iterations==0 || c.linear_tolerance<=0.0 ||
       !(c.intensity_relaxation>0.0 && c.intensity_relaxation<=1.0))
        throw std::invalid_argument("invalid radiation transport controls");
}

struct P1RadiationSolveResult {
    bool converged = false;
    std::size_t iterations = 0;
    double residual = 0.0;
};

// Solve the gray isotropic P1 equation
// div(D grad G) - kappa_a G + 4 kappa_a sigma T^4 = 0,
// with D = 1/[3(kappa_a+kappa_s)]. The supplied scalar BCs are
// boundary conditions on irradiation G [W/m2].
inline P1RadiationSolveResult solve_p1_radiation(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& irradiation,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& radiation_source,
    double absorption,
    double scattering,
    const ScalarBoundaryConditions& irradiation_bcs = {},
    std::size_t max_iterations = 2000,
    double tolerance = 1e-10)
{
    if(!std::isfinite(absorption) || !std::isfinite(scattering) ||
       absorption < 0.0 || scattering < 0.0 ||
       absorption + scattering <= 0.0 ||
       max_iterations == 0 || !std::isfinite(tolerance) || tolerance <= 0.0 ||
       temperature.size()!=mesh.n_cells() ||
       irradiation.size()!=mesh.n_cells() ||
       radiation_source.size()!=mesh.n_cells())
        throw std::invalid_argument("invalid P1 radiation solve inputs");

    const double D = 1.0/(3.0*(absorption+scattering));
    cfdx::core::Field<double,cfdx::core::Location::FACE> zero_flux(
        mesh.n_faces(),"p1_flux","m2/s",1);
    zero_flux.fill(0.0);
    cfdx::core::Field<double,cfdx::core::Location::CELL> su(
        mesh.n_cells(),"p1_source","W/m3",1);
    cfdx::core::Field<double,cfdx::core::Location::CELL> sp(
        mesh.n_cells(),"p1_sp","1/m",1);
    for(std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
        const double T=temperature(cell);
        if(!std::isfinite(T) || T<0.0)
            throw std::invalid_argument("P1 temperature must be finite and non-negative");
        su(cell)=4.0*absorption*blackbody_emissive_power(T);
        sp(cell)=-absorption;
    }

    auto eq=assemble_scalar_equation(
        mesh,geometry,zero_flux,D,su,sp,irradiation_bcs,true);
    cfdx::core::Vector G(mesh.n_cells(),0.0);
    for(std::size_t cell=0; cell<mesh.n_cells(); ++cell)
        G(cell)=irradiation(cell);
    ScalarSolveControls controls;
    controls.max_iterations=max_iterations;
    controls.tolerance=tolerance;
    controls.relaxation=1.0;
    const auto linear=solve_scalar_equation(eq,G,controls);

    P1RadiationSolveResult result;
    result.converged=(linear.status==cfdx::core::SolverStatus::CONVERGED);
    result.iterations=linear.iterations;
    result.residual=linear.residual;
    if(!result.converged)
        return result;

    for(std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
        if(!std::isfinite(G(cell)) || G(cell)<-1e-10)
            throw std::runtime_error("P1 produced non-physical irradiation");
        irradiation(cell)=std::max(0.0,G(cell));
        radiation_source(cell)=absorption*
            (4.0*M_PI*blackbody_intensity(temperature(cell))-irradiation(cell));
    }
    return result;
}

inline RadiationSolveResult solve_participating_radiation(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& irradiation,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& radiation_source,
    const std::vector<DiscreteDirection>& directions,
    const RadiationTransportControls& controls = {},
    const ScalarBoundaryConditions& wall_intensity_bcs = {})
{
    validate_radiation_transport_controls(controls);
    validate_discrete_directions(directions);
    if(temperature.size()!=mesh.n_cells() ||
       irradiation.size()!=mesh.n_cells() ||
       radiation_source.size()!=mesh.n_cells())
        throw std::invalid_argument("radiation field size mismatch");

    const std::size_t nc=mesh.n_cells();
    std::vector<cfdx::core::Field<double,cfdx::core::Location::CELL>> intensities;
    intensities.reserve(directions.size());
    for(std::size_t m=0;m<directions.size();++m) {
        intensities.emplace_back(nc,"I_"+std::to_string(m),"W/m2/sr",1);
        intensities.back().fill(0.0);
    }

    RadiationSolveResult result;
    for(std::size_t iter=1;iter<=controls.max_iterations;++iter) {
        auto old_radiation_source=radiation_source;
        std::vector<double> J(nc,0.0);
        for(std::size_t m=0;m<directions.size();++m)
            for(std::size_t c=0;c<nc;++c)
                J[c]+=directions[m].weight*intensities[m](c)/(4.0*M_PI);

        double max_delta=0.0;
        for(std::size_t m=0;m<directions.size();++m) {
            cfdx::core::Field<double,cfdx::core::Location::FACE> directional_flux(
                mesh.n_faces(),"sI","W/m2",1);
            const auto& d=directions[m];
            for(std::size_t f=0;f<mesh.n_faces();++f)
                directional_flux(f)=d.dx*geometry.face_area_vectors[f].x+
                                    d.dy*geometry.face_area_vectors[f].y+
                                    d.dz*geometry.face_area_vectors[f].z;

            cfdx::core::Field<double,cfdx::core::Location::CELL> source(
                nc,"radiation_source","W/m3/sr",1);
            cfdx::core::Field<double,cfdx::core::Location::CELL> sp(
                nc,"radiation_sp","1/m",1);

            for(std::size_t c=0;c<nc;++c) {
                const double Ib=blackbody_intensity(temperature(c));
                source(c)=controls.absorption*Ib +
                          controls.scattering*J[c];
                sp(c)=-(controls.absorption+controls.scattering);
            }

            auto eq=assemble_scalar_equation(
                mesh,geometry,directional_flux,0.0,source,sp,
                wall_intensity_bcs,true);

            auto old=intensities[m];
            cfdx::core::Vector intensity(nc,0.0);
            for(std::size_t c=0;c<nc;++c) intensity(c)=intensities[m](c);
            ScalarSolveControls sc;
            sc.max_iterations=controls.linear_max_iterations;
            sc.tolerance=controls.linear_tolerance;
            sc.relaxation=controls.intensity_relaxation;
            const auto lr=solve_scalar_equation(eq,intensity,sc);
            for(std::size_t c=0;c<nc;++c) intensities[m](c)=intensity(c);
            if(lr.status!=cfdx::core::SolverStatus::CONVERGED) {
                throw std::runtime_error(
                    "radiation intensity linear solve did not converge: residual=" +
                    std::to_string(lr.residual) +
                    " relative=" + std::to_string(lr.residual_relative) +
                    " iterations=" + std::to_string(lr.iterations));
            }

            for(std::size_t c=0;c<nc;++c)
                max_delta=std::max(max_delta,
                    std::abs(intensities[m](c)-old(c)));
        }

        for(std::size_t c=0;c<nc;++c) {
            double G=0.0;
            for(std::size_t m=0;m<directions.size();++m)
                G+=directions[m].weight*intensities[m](c);
            irradiation(c)=G;
            radiation_source(c)=controls.absorption*
                (4.0*M_PI*blackbody_intensity(temperature(c))-G);
        }

        double max_source_delta = 0.0;
        for(std::size_t c=0;c<nc;++c)
            max_source_delta=std::max(max_source_delta,
                std::abs(radiation_source(c)-old_radiation_source(c)));
        result.history.push_back({iter,max_delta,max_source_delta});
        result.iterations=iter;
        double max_relative_delta=0.0;
        double max_relative_source_delta=0.0;
        for(std::size_t c=0;c<nc;++c) {
            double old_scale=1.0, source_scale=1.0;
            for(std::size_t m=0;m<directions.size();++m) {
                old_scale=std::max(old_scale,std::abs(intensities[m](c)));
            }
            old_scale=std::max(old_scale,std::abs(irradiation(c)));
            source_scale=std::max(source_scale,std::abs(radiation_source(c)));
            max_relative_delta=std::max(max_relative_delta,
                max_delta/old_scale);
            max_relative_source_delta=std::max(max_relative_source_delta,
                max_source_delta/source_scale);
        }
        if(max_relative_delta<=controls.tolerance &&
           max_relative_source_delta<=controls.tolerance) {
            result.converged=true;
            break;
        }
    }
    return result;
}

struct RadiationEnergyCouplingControls {
    RadiationTransportControls radiation;
    EnergySolverControls energy;
    std::size_t max_outer_iterations = 100;
    double tolerance = 1e-8;
};

struct RadiationEnergyCouplingResult {
    bool converged = false;
    std::size_t iterations = 0;
    std::vector<double> source_residuals;
    std::vector<double> energy_balance_residuals;
};

inline RadiationEnergyCouplingResult solve_radiation_energy_coupled(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& non_radiative_source,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& irradiation,
    const std::vector<DiscreteDirection>& directions,
    const RadiationEnergyCouplingControls& controls = {},
    const ScalarBoundaryConditions& radiation_bcs = {},
    const ScalarBoundaryConditions& thermal_bcs = {})
{
    if(controls.max_outer_iterations==0 || controls.tolerance<=0.0)
        throw std::invalid_argument("invalid radiation-energy coupling controls");

    const std::size_t nc=mesh.n_cells();
    cfdx::core::Field<double,cfdx::core::Location::CELL> qrad(
        nc,"qrad","W/m3",1);
    cfdx::core::Field<double,cfdx::core::Location::CELL> qtot(
        nc,"qtot","W/m3",1);

    RadiationEnergyCouplingResult result;
    for(std::size_t iter=1;iter<=controls.max_outer_iterations;++iter) {
        cfdx::core::Field<double,cfdx::core::Location::CELL> oldT=temperature;
        cfdx::core::Field<double,cfdx::core::Location::CELL> old_qrad=qrad;
        cfdx::core::Field<double,cfdx::core::Location::CELL> source(
            nc,"radiation_source","W/m3",1);

        auto rr=solve_participating_radiation(
            mesh,geometry,temperature,irradiation,qrad,directions,
            controls.radiation,radiation_bcs);
        if(!rr.converged)
            throw std::runtime_error("radiation inner solve did not converge");

        for(std::size_t c=0;c<nc;++c)
            source(c)=non_radiative_source(c)+qrad(c);

        auto er=solve_energy(
            mesh,geometry,mass_flux,temperature,source,
            controls.energy,thermal_bcs);
        if(!er.converged)
            throw std::runtime_error("energy inner solve did not converge");

        double max_delta=0.0;
        for(std::size_t c=0;c<nc;++c)
            max_delta=std::max(max_delta,std::abs(temperature(c)-oldT(c)));

        double qrad_delta=0.0;
        for(std::size_t c=0;c<nc;++c)
            qrad_delta=std::max(qrad_delta,std::abs(qrad(c)-old_qrad(c)));
        const double energy_balance_residual = er.history.empty()
            ? std::numeric_limits<double>::infinity()
            : er.history.back().energy_imbalance;
        result.source_residuals.push_back(qrad_delta);
        result.energy_balance_residuals.push_back(energy_balance_residual);
        result.iterations=iter;
        if(max_delta<=controls.tolerance &&
           qrad_delta<=controls.tolerance &&
           energy_balance_residual<=controls.tolerance) {
            result.converged=true;
            break;
        }
    }
    return result;
}

} // namespace cfdx::physics
