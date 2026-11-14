#pragma once
// CFDX V54 - research/qualification-only independent ordinate-by-ordinate
// transport auditor. Does NOT invoke the V49-V53 HO/LO correction/identity.
// Input is a candidate x=(J, reflected wall R), not an accepted solution.
#include "cfdx/physics/dom_holo_variable_v50.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct DomV54AngularAudit {
    bool accepted=false;
    double max_ordinate_equation=std::numeric_limits<double>::infinity();
    double max_wall_reflection=std::numeric_limits<double>::infinity();
    double max_cell_energy=std::numeric_limits<double>::infinity();
    double global_energy=std::numeric_limits<double>::infinity();
    double max_fixed_moment=std::numeric_limits<double>::infinity();
    double face_partition_error=std::numeric_limits<double>::infinity();
    double maximum_wall_abs_watts=0.;
};

inline DomV54AngularAudit audit_dom_angular_v54(
    DomMomentSchurSystemV44& ho,const std::vector<double>& x,
    double physical_tolerance,std::size_t simulated_partitions=1) {
    if(!(physical_tolerance>0.) || !std::isfinite(physical_tolerance) ||
       simulated_partitions==0 || simulated_partitions>128 ||
       x.size()!=ho.dimension())
        throw std::invalid_argument("V54: invalid angular audit configuration");
    for(const double z:x)if(!(z>=0.) || !std::isfinite(z))
        throw std::invalid_argument("V54: invalid HO candidate");
    const auto& mesh=ho.source_mesh();
    const auto& geo=ho.source_geometry();
    const auto& directions=ho.directions();
    const auto& ka=ho.cell_absorption();
    const auto& ks=ho.cell_scattering();
    const std::size_t nc=ho.n_cells(),nf=mesh.n_faces();
    std::vector<DomMomentSchurSystemV44::CellField> fields;
    const auto projected=ho.project(x,true,&fields);
    if(fields.size()!=directions.size() || projected.size()!=x.size())
        throw std::runtime_error("V54: angular transport returned incomplete state");
    std::vector<const DomMomentSchurSystemV44::WallFace*> faces(nf,nullptr);
    for(const auto& wall:ho.wall_faces()) {
        if(wall.face>=nf || faces[wall.face])
            throw std::runtime_error("V54: wall topology mismatch");
        faces[wall.face]=&wall;
    }
    std::vector<long double> cell_flux(nc,0.),cell_flux_abs(nc,0.);
    std::vector<long double> outgoing(nf,0.);
    std::vector<long double> face_q(nf,0.);
    std::vector<long double> partition_flux(nc,0.);
    std::vector<long double> cell_moment(nc,0.);
    DomV54AngularAudit a;
    a.max_ordinate_equation=0.;a.max_wall_reflection=0.;
    a.max_cell_energy=0.;a.global_energy=0.;
    a.max_fixed_moment=0.;a.face_partition_error=0.;
    constexpr long double fourpi=12.566370614359172953850573533118L;
    for(std::size_t m=0;m<directions.size();++m) {
        const auto& d=directions[m];
        std::vector<long double> angle_div(nc,0.),angle_abs(nc,0.);
        for(std::size_t f=0;f<nf;++f) {
            const auto& sf=geo.face_area_vectors[f];
            const long double mu=static_cast<long double>(d.dx)*sf.x+
                    static_cast<long double>(d.dy)*sf.y+
                    static_cast<long double>(d.dz)*sf.z;
            const std::size_t c=mesh.ownership().owner(f);
            const auto nb=mesh.ownership().neighbour(f);
            long double upwind=0.;
            if(mu>=0.) {
                upwind=fields[m](c);
                if(nb<0)outgoing[f]+=static_cast<long double>(d.weight)*mu*upwind;
            } else if(nb>=0) {
                upwind=fields[m](static_cast<std::size_t>(nb));
            } else {
                const auto* wall=faces[f];
                if(!wall || !(wall->hminus>0.))
                    throw std::runtime_error("V54: angular boundary is missing");
                const long double reflected=wall->reflect_index<ho.n_reflecting_faces()?
                    (1.-wall->emissivity)*(wall->hplus/wall->hminus)*
                    x[nc+wall->reflect_index]:0.;
                upwind=wall->emissivity*wall->ib+reflected;
            }
            const long double q=static_cast<long double>(d.weight)*mu*upwind;
            angle_div[c]+=q;angle_abs[c]+=std::abs(q);
            face_q[f]+=q;
            if(nb>=0) {
                const auto j=static_cast<std::size_t>(nb);
                angle_div[j]-=q;angle_abs[j]+=std::abs(q);
            }
        }
        for(std::size_t c=0;c<nc;++c) {
            const long double value=fields[m](c);
            if(!std::isfinite(static_cast<double>(value)) || value<0.)
                throw std::runtime_error("V54: nonphysical angular intensity");
            cell_moment[c]+=static_cast<long double>(d.weight)*value;
            const long double gas=blackbody_intensity(ho.source_temperature()(c));
            const long double v=geo.cell_volumes[c];
            const long double loss=static_cast<long double>(d.weight)*v*
                ((ka[c]+ks[c])*value-ka[c]*gas-ks[c]*x[c]);
            const long double r=angle_div[c]+loss;
            const long double scale=std::max({1.L,angle_abs[c],std::abs(loss),
                static_cast<long double>(d.weight)*v*(ka[c]+ks[c])*value});
            a.max_ordinate_equation=std::max(a.max_ordinate_equation,
                                              static_cast<double>(std::abs(r)/scale));
        }
    }
    for(std::size_t f=0;f<nf;++f) {
        const auto c=mesh.ownership().owner(f);
        const auto nb=mesh.ownership().neighbour(f);
        cell_flux[c]+=face_q[f];cell_flux_abs[c]+=std::abs(face_q[f]);
        if(nb>=0) {
            const auto j=static_cast<std::size_t>(nb);
            cell_flux[j]-=face_q[f];cell_flux_abs[j]+=std::abs(face_q[f]);
        } else {
            if(!faces[f])throw std::runtime_error("V54: missing boundary face");
            a.maximum_wall_abs_watts=std::max(a.maximum_wall_abs_watts,
                                                static_cast<double>(std::abs(face_q[f])));
            const auto& w=*faces[f];
            if(w.reflect_index<ho.n_reflecting_faces()){
                const double Rout=static_cast<double>(outgoing[f]/(w.hplus*w.area));
                a.max_wall_reflection=std::max(a.max_wall_reflection,
                    std::abs(x[nc+w.reflect_index]-Rout)/
                    std::max({1.,std::abs(Rout),x[nc+w.reflect_index]}));
            }
        }
    }
    // Simulate distinct rank partitions, without claiming an MPI/GPU run.
    // Accumulate each owned face exactly once and reproduce the signed
    // owner/neighbor sum. The partition count changes the accumulation order.
    for(std::size_t rank=0;rank<simulated_partitions;++rank)
        for(std::size_t f=rank;f<nf;f+=simulated_partitions) {
            const auto c=mesh.ownership().owner(f);
            partition_flux[c]+=face_q[f];
            const auto nb=mesh.ownership().neighbour(f);
            if(nb>=0)partition_flux[static_cast<std::size_t>(nb)]-=face_q[f];
        }
    long double total_boundary=0.,total_volume=0.,total_scale=1.;
    for(std::size_t f=0;f<nf;++f)if(mesh.ownership().neighbour(f)<0) {
        total_boundary+=face_q[f];total_scale+=std::abs(face_q[f]);
    }
    for(std::size_t c=0;c<nc;++c) {
        const long double J=cell_moment[c]/fourpi;
        const long double ib=blackbody_intensity(ho.source_temperature()(c));
        const long double emission=fourpi*ka[c]*geo.cell_volumes[c]*(J-ib);
        const long double r=cell_flux[c]+emission;
        const long double scale=std::max({1.L,cell_flux_abs[c],std::abs(emission)});
        a.max_cell_energy=std::max(a.max_cell_energy,
            static_cast<double>(std::abs(r)/scale));
        a.face_partition_error=std::max(a.face_partition_error,
            static_cast<double>(std::abs(cell_flux[c]-partition_flux[c])/scale));
        a.max_fixed_moment=std::max(a.max_fixed_moment,
            std::abs(x[c]-static_cast<double>(J)) /
            std::max({1.,std::abs(x[c]),std::abs(static_cast<double>(J))}));
        total_volume+=emission;total_scale+=std::abs(emission);
    }
    a.global_energy=static_cast<double>(std::abs(total_boundary+total_volume)/total_scale);
    a.accepted=a.max_ordinate_equation<=physical_tolerance &&
               a.max_wall_reflection<=physical_tolerance &&
               a.max_cell_energy<=physical_tolerance &&
               a.global_energy<=physical_tolerance &&
               a.max_fixed_moment<=physical_tolerance &&
               a.face_partition_error<=physical_tolerance;
    return a;
}
} // namespace cfdx::physics
