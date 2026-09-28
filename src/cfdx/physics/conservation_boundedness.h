#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

struct TransportConservationReport {
    double boundary_flux = 0.0;
    double source_sum = 0.0;
    double accumulation = 0.0;
    double residual = 0.0;
    double normalized_residual = 0.0;
    double max_cell_residual = 0.0;
    std::size_t worst_cell = 0;
    std::size_t nonfinite_cells = 0;
    std::size_t nonfinite_sources = 0;

    bool finite() const noexcept {
        return nonfinite_cells == 0 && nonfinite_sources == 0 &&
               std::isfinite(boundary_flux) && std::isfinite(source_sum) &&
               std::isfinite(accumulation) && std::isfinite(residual);
    }
};

inline TransportConservationReport audit_transport_balance(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& face_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& source,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>* accumulation = nullptr)
{
    const std::size_t nf=mesh.n_faces(), nc=mesh.n_cells();
    if(face_flux.size()!=nf || face_flux.dimension()!=1 ||
       source.size()!=nc || source.dimension()!=1)
        throw std::invalid_argument("audit_transport_balance: field dimensions do not match mesh");
    if(accumulation && (accumulation->size()!=nc || accumulation->dimension()!=1))
        throw std::invalid_argument("audit_transport_balance: invalid accumulation field");

    TransportConservationReport r;
    const auto& own=mesh.ownership();
    std::vector<double> cell(nc,0.0);
    for(std::size_t f=0;f<nf;++f) {
        const double phi=face_flux(f);
        if(!std::isfinite(phi)) { ++r.nonfinite_cells; continue; }
        const std::size_t o=own.owner(f);
        if(o>=nc) throw std::runtime_error("audit_transport_balance: invalid owner");
        cell[o]+=phi;
        const auto n=own.neighbour(f);
        if(n>=0) {
            const auto ni=static_cast<std::size_t>(n);
            if(ni>=nc) throw std::runtime_error("audit_transport_balance: invalid neighbour");
            cell[ni]-=phi;
        } else {
            r.boundary_flux+=phi;
        }
    }
    for(std::size_t c=0;c<nc;++c) {
        const double s=source(c);
        if(!std::isfinite(s)) { ++r.nonfinite_sources; continue; }
        r.source_sum+=s;
        const double a=accumulation ? (*accumulation)(c) : 0.0;
        if(!std::isfinite(a)) { ++r.nonfinite_cells; continue; }
        r.accumulation+=a;
        const double rc=cell[c]+s-a;
        r.residual+=rc;
        if(std::abs(rc)>r.max_cell_residual) {
            r.max_cell_residual=std::abs(rc);
            r.worst_cell=c;
        }
    }
    const double scale=std::max({std::abs(r.boundary_flux),
                                 std::abs(r.source_sum),
                                 std::abs(r.accumulation),1.0});
    r.normalized_residual=std::abs(r.residual)/scale;
    return r;
}

struct ScalarBoundReport {
    double minimum=std::numeric_limits<double>::infinity();
    double maximum=-std::numeric_limits<double>::infinity();
    double lower=0.0;
    double upper=std::numeric_limits<double>::infinity();
    std::size_t nonfinite=0;
    std::size_t below=0;
    std::size_t above=0;
    std::size_t worst_cell=0;
    double worst_violation=0.0;

    bool bounded() const noexcept {
        return nonfinite==0 && below==0 && above==0;
    }
};

inline ScalarBoundReport audit_scalar_bounds(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& field,
    double lower, double upper=std::numeric_limits<double>::infinity(),
    double tolerance=0.0)
{
    if(field.dimension()!=1 || field.size()==0)
        throw std::invalid_argument("audit_scalar_bounds: expected non-empty scalar cell field");
    if(!std::isfinite(lower) || (std::isfinite(upper)&&lower>upper) ||
       (std::isinf(upper)&&upper<0.0) || !std::isfinite(tolerance) || tolerance<0.0)
        throw std::invalid_argument("audit_scalar_bounds: invalid bounds");
    ScalarBoundReport r{.lower=lower,.upper=upper};
    const double lo=lower-tolerance, hi=upper+tolerance;
    for(std::size_t c=0;c<field.size();++c) {
        const double v=field(c);
        if(!std::isfinite(v)) { ++r.nonfinite; continue; }
        r.minimum=std::min(r.minimum,v);
        r.maximum=std::max(r.maximum,v);
        const double violation=std::max(lo-v, std::isfinite(hi)?v-hi:0.0);
        if(v<lo) ++r.below;
        if(std::isfinite(hi)&&v>hi) ++r.above;
        if(violation>r.worst_violation) {
            r.worst_violation=violation;
            r.worst_cell=c;
        }
    }
    return r;
}

inline ScalarBoundReport audit_positive_scalar(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& field,
    double floor=0.0)
{
    if(!std::isfinite(floor) || floor<0.0)
        throw std::invalid_argument("audit_positive_scalar: invalid floor");
    return audit_scalar_bounds(field,floor);
}

} // namespace cfdx::physics
