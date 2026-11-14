#pragma once
// V39 opt-in: conservative gray P1 + diffuse-gray Marshak wall closure.
// The radiative outward flux at a wall is
// q_out = eps/[2(2-eps)] * (G_surface - 4*sigma*T_wall^4).
// Its cell-centre FV conductance is A/(d_normal/D + 1/c), with
// D=1/(3*(kappa_a+kappa_s)), c=eps/[2(2-eps)].
// All fluxes use the SAME coefficient as assembly. No global energy
// correction or physical gate relaxation is applied.
#include "cfdx/physics/radiation_solver.h"
#include <map>
#include <string>
#include <vector>
#include <limits>
#include <numeric>

namespace cfdx::physics {

struct P1MarshakWall {
    double emissivity=1.0;
    double temperature=300.0;
};
using P1MarshakWalls=std::map<std::string,P1MarshakWall>;
// Optional face-resolved wall properties for non-isothermal patches.
// Keys are actual mesh face IDs, not patch-relative indices.
using P1MarshakFaceWalls=std::map<std::size_t,P1MarshakWall>;

struct P1MarshakControls {
    std::size_t max_outer_iterations=40;
    std::size_t linear_max_iterations=4000;
    double linear_tolerance=1.e-12;
    double physical_tolerance=1.e-10;
    double balance_tolerance=1.e-10;
    // V40 experimental: include wall Robin normals in the boundary-aware WLS
    // and its deferred nonorthogonal correction, not a zero-Neumann surrogate.
    bool skew_marshak_robin_v40=false;
    std::size_t robin_gradient_sweeps=12;
};
struct P1MarshakResult {
    bool converged=false;
    std::size_t outer_iterations=0;
    std::size_t linear_iterations=0;
    double maximum_normalized_cell_residual=std::numeric_limits<double>::infinity();
    double relative_global_balance=std::numeric_limits<double>::infinity();
    double total_boundary_outward_W=0.0;
    double volume_emission_minus_absorption_W=0.0;
    std::map<std::string,double> patch_outward_W;
    std::map<std::size_t,double> face_outward_W;
    std::size_t maximum_robin_gradient_sweeps=0;
};

inline P1MarshakResult solve_p1_radiation_marshak(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& irradiation,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& radiation_source,
    const std::vector<double>& absorption,
    const std::vector<double>& scattering,
    const P1MarshakWalls& walls,
    const P1MarshakControls& controls={},
    const P1MarshakFaceWalls* face_walls=nullptr)
{
    using namespace cfdx::core;
    const std::size_t nc=mesh.n_cells(), nf=mesh.n_faces();
    if (nc==0 || nf==0 || temperature.size()!=nc || irradiation.size()!=nc ||
        radiation_source.size()!=nc || absorption.size()!=nc || scattering.size()!=nc ||
        geometry.cell_volumes.size()!=nc || geometry.cell_centres.size()!=nc ||
        geometry.face_centres.size()!=nf || geometry.face_area_vectors.size()!=nf ||
        controls.max_outer_iterations==0 || controls.linear_max_iterations==0 ||
        (controls.skew_marshak_robin_v40 && controls.robin_gradient_sweeps==0) ||
        !(controls.linear_tolerance>0) || !std::isfinite(controls.linear_tolerance) ||
        !(controls.physical_tolerance>0) || !std::isfinite(controls.physical_tolerance) ||
        !(controls.balance_tolerance>0) || !std::isfinite(controls.balance_tolerance))
        throw std::invalid_argument("P1 Marshak: invalid sizes or solve controls");

    std::vector<double> diffusion(nc), driving(nc), staged(nc);
    Field<double,Location::CELL> su(nc,"p1_marshak_su","W/m3",1);
    Field<double,Location::CELL> sp(nc,"p1_marshak_sp","1/m",1);
    Field<double,Location::FACE> face_flux(nf,"p1_marshak_phi","W",1);
    face_flux.fill(0.0);
    bool volume_anchor=false;
    for (std::size_t i=0;i<nc;++i) {
        const double a=absorption[i], s=scattering[i], T=temperature(i), G=irradiation(i);
        if (!std::isfinite(a) || !std::isfinite(s) || a<0 || s<0 || !(a+s>0) ||
            !std::isfinite(T) || T<0 || !std::isfinite(G) || G<0 ||
            !(geometry.cell_volumes[i]>0) || !std::isfinite(geometry.cell_volumes[i]))
            throw std::invalid_argument("P1 Marshak: nonphysical cell state/coefficients");
        diffusion[i]=1.0/(3.0*(a+s));
        driving[i]=4.0*blackbody_emissive_power(T);
        su(i)=a*driving[i];
        sp(i)=-a;
        staged[i]=G;
        volume_anchor=volume_anchor || a>0;
        if (!(diffusion[i]>0) || !std::isfinite(diffusion[i]) ||
            !std::isfinite(driving[i]) || !std::isfinite(su(i)))
            throw std::invalid_argument("P1 Marshak: invalid diffusion/emission");
    }
    if (walls.size()!=mesh.boundary().n_patches())
        throw std::invalid_argument("P1 Marshak: missing/extra boundary patch");
    std::vector<double> extra_diag(nc,0),extra_rhs(nc,0);
    struct WallFace {std::size_t face,owner; double conductance,driving,area,diffusion;
                     std::string patch;};
    std::vector<WallFace> boundary_faces;
    bool boundary_anchor=false;
    std::vector<unsigned char> face_covered(nf,0);
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& patch=mesh.boundary().patch(p);
        auto it=walls.find(patch.name);
        if(it==walls.end()) throw std::invalid_argument("P1 Marshak: unknown wall patch");
        if (!std::isfinite(it->second.emissivity) || it->second.emissivity<0 ||
            it->second.emissivity>1 || !std::isfinite(it->second.temperature) ||
            it->second.temperature<0)
            throw std::invalid_argument("P1 Marshak: invalid patch wall");
        if (patch.type==PatchType::EMPTY)
            throw std::invalid_argument("P1 Marshak: EMPTY patch needs separate reduced-dimension policy");
        for(const auto f:patch.face_ids) {
            if(f>=nf || mesh.ownership().neighbour(f)>=0)
                throw std::invalid_argument("P1 Marshak: invalid wall face");
            if(face_covered[f]++)
                throw std::invalid_argument("P1 Marshak: duplicate boundary face");
            P1MarshakWall fw=it->second;
            if(face_walls) {
                const auto override=face_walls->find(f);
                if(override!=face_walls->end()) fw=override->second;
            }
            if(!std::isfinite(fw.emissivity) || fw.emissivity<0 || fw.emissivity>1 ||
               !std::isfinite(fw.temperature) || fw.temperature<0)
                throw std::invalid_argument("P1 Marshak: invalid face emissivity/temperature");
            const double c=fw.emissivity/(2.0*(2.0-fw.emissivity));
            const double Gwall=4.0*blackbody_emissive_power(fw.temperature);
            if(!std::isfinite(Gwall)) throw std::invalid_argument("P1 Marshak: wall emission overflow");
            const auto owner=static_cast<std::size_t>(mesh.ownership().owner(f));
            if(owner>=nc) throw std::invalid_argument("P1 Marshak: invalid owner");
            const double area=geometry.face_area_vectors[f].mag();
            const auto outward=geometry.face_centres[f]-geometry.cell_centres[owner];
            if(!(geometry.face_area_vectors[f].dot(outward)>0))
                throw std::invalid_argument("P1 Marshak: reversed boundary normal");
            const double d=boundary_normal_distance(geometry,f,owner);
            if(!(area>0) || !std::isfinite(area) || !(d>0) || !std::isfinite(d))
                throw std::invalid_argument("P1 Marshak: invalid boundary geometry");
            const double U = c>0.0 ? (c*diffusion[owner])/(diffusion[owner]+c*d) : 0.0;
            const double conductance=area*U;
            if(!(conductance>=0) || !std::isfinite(conductance))
                throw std::invalid_argument("P1 Marshak: boundary conductance overflow");
            extra_diag[owner]+=conductance;
            extra_rhs[owner]+=conductance*Gwall;
            if(!std::isfinite(extra_diag[owner]) || !std::isfinite(extra_rhs[owner]))
                throw std::invalid_argument("P1 Marshak: boundary assembly overflow");
            boundary_anchor=boundary_anchor || (conductance>0);
            boundary_faces.push_back({f,owner,conductance,Gwall,area,diffusion[owner],patch.name});
        }
    }
    for(std::size_t f=0;f<nf;++f)
        if(mesh.ownership().neighbour(f)<0 && !face_covered[f])
            throw std::invalid_argument("P1 Marshak: unassigned boundary face");
    if(face_walls) for(const auto& [f,unused]:*face_walls) {
        (void)unused;
        if(f>=nf || !face_covered[f])
            throw std::invalid_argument("P1 Marshak: override is not a boundary face");
    }
    if(!volume_anchor && !boundary_anchor)
        throw std::invalid_argument("P1 Marshak: unanchored reflecting, nonabsorbing domain");

    ScalarSolveControls linear_controls;
    linear_controls.max_iterations=controls.linear_max_iterations;
    linear_controls.tolerance=controls.linear_tolerance;
    linear_controls.relaxation=1.0;
    // Always stage updates. A failed linear/outer solve cannot alter committed
    // irradiation or energy source, nor publish partial wall power.
    Field<double,Location::CELL> current(nc,"p1_marshak_candidate","W/m2",1);
    std::vector<double> wall_tangential(nf,0.0);
    std::size_t last_robin_sweeps=0;
    // Recompute all wall corrections from one coherent trial field, including
    // for the independent post-solve residual. A wall is represented once by
    // one implicit Robin conductance and one explicit tangential correction.
    const auto assemble=[&](const std::vector<double>& x){
        for(std::size_t i=0;i<nc;++i) current(i)=x[i];
        if(!controls.skew_marshak_robin_v40)
            return assemble_scalar_equation(mesh,geometry,face_flux,0.0,su,sp,
                ScalarBoundaryConditions{},true,nullptr,&extra_diag,&extra_rhs,
                &diffusion,ConvectionScheme::UPWIND,&current);
        ScalarBoundaryFaceConditions face_conditions(nf);
        for(const auto& wf:boundary_faces) {
            face_conditions.specified[wf.face]=1;
            face_conditions.conditions[wf.face]={
                ScalarBoundaryFaceCondition::Type::FIXED_GRADIENT,0.0};
            wall_tangential[wf.face]=0.0;
        }
        std::vector<double> normal_flux(nf,0.0);
        last_robin_sweeps=0;
        bool robin_gradient_settled=false;
        for(const auto& wf:boundary_faces)
            normal_flux[wf.face]=wf.conductance*(x[wf.owner]-wf.driving);
        // Update the actual Robin normal gradient while WLS reconstructs the
        // tangential gradient from cell neighbours + matching wall flux.
        // On a skewed boundary, phi_wall = phi_owner + tangential - q*d/(D*A).
        // Therefore q = U*A*(phi_owner + tangential - phi_wall_radiative).
        for(std::size_t pass=0;pass<controls.robin_gradient_sweeps;++pass) {
            last_robin_sweeps=pass+1;
            for(const auto& wf:boundary_faces)
                face_conditions.conditions[wf.face].value =
                    -normal_flux[wf.face]/(wf.area*wf.diffusion);
            const auto gradients=compute_boundary_aware_wls_gradient(
                current,mesh,geometry,ScalarBoundaryConditions{},nullptr,
                &face_conditions);
            double maximum_change=0.0;
            for(const auto& wf:boundary_faces) {
                const auto n=geometry.face_area_vectors[wf.face]*(1.0/wf.area);
                const auto r=geometry.face_centres[wf.face]-geometry.cell_centres[wf.owner];
                const auto tangent=r-n*r.dot(n);
                const cfdx::core::Vec3 g{
                    gradients.component_data(0)[wf.owner],
                    gradients.component_data(1)[wf.owner],
                    gradients.component_data(2)[wf.owner]};
                const double correction=g.dot(tangent);
                const double q=wf.conductance*(x[wf.owner]+correction-wf.driving);
                if(!std::isfinite(correction) || !std::isfinite(q))
                    throw std::runtime_error("P1 Marshak: invalid skew Robin reconstruction");
                maximum_change=std::max(maximum_change,std::abs(q-normal_flux[wf.face])/
                    std::max({1.0,std::abs(q),std::abs(normal_flux[wf.face])}));
                normal_flux[wf.face]=q;
                wall_tangential[wf.face]=correction;
            }
            if(maximum_change<=1.e-12){robin_gradient_settled=true;break;}
        }
        if(!robin_gradient_settled)
            throw std::runtime_error("P1 Marshak: boundary Robin gradient did not contract");
        // Need the same normal gradient in the WLS and FV assembly, but its
        // temporary Neumann watts must be cancelled: Robin is already present
        // as one implicit conductance. No doubled boundary energy transfer.
        std::vector<double> robin_rhs(nc,0.0);
        for(std::size_t i=0;i<nc;++i)robin_rhs[i]=extra_rhs[i];
        for(const auto& wf:boundary_faces) {
            const double grad=-normal_flux[wf.face]/(wf.area*wf.diffusion);
            face_conditions.conditions[wf.face].value=grad;
            robin_rhs[wf.owner]-=wf.conductance*wall_tangential[wf.face]
                                   +wf.area*wf.diffusion*grad;
        }
        return assemble_scalar_equation(mesh,geometry,face_flux,0.0,su,sp,
            ScalarBoundaryConditions{},true,nullptr,&extra_diag,&robin_rhs,
            &diffusion,ConvectionScheme::UPWIND,&current,&face_conditions);
    };
    P1MarshakResult result;
    for(std::size_t sweep=0;sweep<controls.max_outer_iterations;++sweep){
        const auto eq=assemble(staged);
        result.maximum_robin_gradient_sweeps=std::max(
            result.maximum_robin_gradient_sweeps,last_robin_sweeps);
        Vector x(nc,0.0);
        for(std::size_t i=0;i<nc;++i)x(i)=staged[i];
        const auto linear=solve_scalar_equation(eq,x,linear_controls);
        result.outer_iterations=sweep+1;
        result.linear_iterations+=linear.iterations;
        if(linear.status!=SolverStatus::CONVERGED) return result;
        for(std::size_t i=0;i<nc;++i){
            if(!std::isfinite(x(i)) || x(i)<0)
                throw std::runtime_error("P1 Marshak: negative/nonfinite irradiation");
            staged[i]=x(i);
        }
        // Reassemble the FV operator on the candidate (not the stale, frozen
        // nonorthogonal correction) and evaluate its independent A*x-b defect.
        const auto reevaluated=assemble(staged);
        result.maximum_robin_gradient_sweeps=std::max(
            result.maximum_robin_gradient_sweeps,last_robin_sweeps);
        long double maximal_relative=0;
        const auto* row=reevaluated.matrix.row_offsets_data();
        const auto* col=reevaluated.matrix.columns_data();
        const auto* val=reevaluated.matrix.values_data();
        for(std::size_t i=0;i<nc;++i){
            long double ax=0,scale=1;
            for(std::size_t j=row[i];j<row[i+1];++j){
                const long double term=static_cast<long double>(val[j])*staged[col[j]];
                ax+=term;
                scale+=std::abs(term);
            }
            const long double defect=ax-reevaluated.rhs(i);
            scale+=std::abs(reevaluated.rhs(i));
            maximal_relative=std::max(maximal_relative,std::abs(defect)/scale);
        }
        result.maximum_normalized_cell_residual=static_cast<double>(maximal_relative);
        long double wallpower=0,volume_power=0,absolute_power=0;
        std::map<std::string,long double> patchpower;
        std::map<std::size_t,double> facepower;
        for(const auto& entry:boundary_faces) {
            const long double q=static_cast<long double>(entry.conductance)*
                (staged[entry.owner]+wall_tangential[entry.face]-entry.driving);
            wallpower+=q;
            absolute_power+=std::abs(q);
            patchpower[entry.patch]+=q;
            facepower[entry.face]=static_cast<double>(q);
        }
        for(std::size_t i=0;i<nc;++i){
            const long double emission=static_cast<long double>(absorption[i])*
                (driving[i]-staged[i])*geometry.cell_volumes[i];
            volume_power+=emission;
            absolute_power+=std::abs(emission);
        }
        result.total_boundary_outward_W=static_cast<double>(wallpower);
        result.volume_emission_minus_absorption_W=static_cast<double>(volume_power);
        result.relative_global_balance=static_cast<double>(std::abs(wallpower-volume_power)/
            std::max(1.0L,absolute_power));
        if (!std::isfinite(result.relative_global_balance) ||
            !std::isfinite(result.maximum_normalized_cell_residual))
            throw std::runtime_error("P1 Marshak: physical audit overflow");
        if(result.maximum_normalized_cell_residual<=controls.physical_tolerance &&
           result.relative_global_balance<=controls.balance_tolerance){
            std::vector<double> staged_source(nc);
            for(std::size_t i=0;i<nc;++i){
                staged_source[i]=absorption[i]*(driving[i]-staged[i]);
                if(!std::isfinite(staged_source[i]))
                    throw std::runtime_error("P1 Marshak: invalid staged energy source");
            }
            for(std::size_t i=0;i<nc;++i){
                irradiation(i)=staged[i];
                radiation_source(i)=staged_source[i];
            }
            for(const auto& [name,power]:patchpower)
                result.patch_outward_W[name]=static_cast<double>(power);
            result.face_outward_W=std::move(facepower);
            result.converged=true;
            return result;
        }
    }
    return result;
}
}
