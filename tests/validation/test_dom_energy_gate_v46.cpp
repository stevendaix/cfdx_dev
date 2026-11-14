// V46 research-only physical acceptance for extreme scattering on an oblique
// 3D finite-volume grid. No thresholds or production DOM settings are changed.
#include "cfdx/physics/dom_face_moment_v45.h"
#include "cfdx/physics/dom_quadrature_v39.h"
#include "thermal_3d_mesh_v36.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
using Cell=DomMomentSchurSystemV44::CellField;
void require(bool ok,const char* why) { if(!ok)throw std::runtime_error(why); }
Mesh grid() {
    auto mesh=cfdx::testing::thermal_v36::make_mesh(3,3,3,.018,0.,.03);
    const double ca=std::cos(.41),sa=std::sin(.41);
    const double cb=std::cos(.33),sb=std::sin(.33);
    for(std::size_t p=0;p<mesh.points().size();++p) {
        const double x=mesh.points().x_data()[p];
        const double y=mesh.points().y_data()[p];
        const double z=mesh.points().z_data()[p];
        const double xr=ca*x-sa*y,yr=sa*x+ca*y;
        mesh.points().set(p,cb*xr+sb*z,yr,-sb*xr+cb*z);
    }
    return mesh;
}
DomWallBoundaryConditions wall_data(const Mesh& mesh) {
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& name=mesh.boundary().patch(p).name;
        const double T=name=="west"?540.:name=="east"?285.:365.;
        const double eps=name=="west"?.35:name=="east"?.65:
            name=="south"?.90:name=="north"?1.:name=="bottom"?0.:1.;
        walls[name]={eps,T};
    }
    return walls;
}
void scenario(bool preconditioned,double scattering) {
    const auto mesh=grid(); const auto geo=build_fv_geometry(mesh);
    const auto dirs=make_dom_tensor_gauss_legendre_v39(4,12);
    const auto walls=wall_data(mesh);
    Cell temperature(mesh.n_cells(),"gas","K",1);
    for(std::size_t c=0;c<mesh.n_cells();++c)
        temperature(c)=370.+30.*geo.cell_centres[c].x+10.*geo.cell_centres[c].y;
    RadiationTransportControls ctl;
    ctl.absorption=.1;ctl.scattering=scattering;
    ctl.tolerance=scattering>100.?1.e-10:1.e-9;ctl.max_iterations=500;
    ctl.linear_tolerance=1.e-12;ctl.linear_max_iterations=2000;
    // Original V45 inner-solve target and Schur: independent energy fails.
    DomMomentSchurSystemV44 original(mesh,geo,temperature,dirs,walls,ctl);
    const auto baseline=original.solve(false,12,160);
    require(baseline.converged,"V46: baseline Schur did not algebraically converge");
    auto baseline_x=baseline.scalar_mean_intensity;
    baseline_x.insert(baseline_x.end(),baseline.reflected_face_mean_intensity.begin(),
                      baseline.reflected_face_mean_intensity.end());
    std::vector<Cell> baseline_I;
    (void)original.project(baseline_x,true,&baseline_I);
    const auto baseline_report=audit_dom_face_moments_v45(
        mesh,geo,temperature,dirs,walls,ctl,original,baseline_x,baseline_I);
    if(scattering>100.)
        require(baseline_report.global_energy_relative>ctl.tolerance,
                "V46: V45 baseline negative control unexpectedly passed");
    else
        require(baseline_report.global_energy_relative<=ctl.tolerance,
                "V46: moderate-scattering baseline did not pass the independent gate");
    if(scattering>100. && preconditioned) {
        // A stricter outer Schur residual alone does not cure the insufficient
        // inner angular FV solve: the old inner accuracy must fail closed.
        const auto loose=original.solve(false,12,100,1.,{},
            [&](const std::vector<double>& x,const std::vector<Cell>& I) {
                const auto r=audit_dom_face_moments_v45(
                    mesh,geo,temperature,dirs,walls,ctl,original,x,I);
                return dom_face_moment_qualified_v45(r,ctl.tolerance);
            });
        require(!loose.converged && loose.scalar_mean_intensity.empty(),
                "V46: unchanged inner accuracy falsely promoted an energy failure");
        std::cout<<"V46_OLD_INNER rejected=1 schur="
                 <<loose.true_schur_relative<<" steps="<<loose.krylov_iterations<<'\n';
    }
    // Inner angular FV accuracy is tightened independently of the physical
    // acceptance criterion, NOT loosened or changed in production.
    auto strict_controls=ctl;
    strict_controls.linear_tolerance=1.e-14;
    DomMomentSchurSystemV44 sys(mesh,geo,temperature,dirs,walls,strict_controls);
    DomLocalAngularJacobiV45 local(sys,geo,dirs,strict_controls);
    std::size_t attempts=0;
    double latest_energy=0.;
    const auto face_gate=[&](const std::vector<double>& x,
                             const std::vector<Cell>& I) {
        const auto r=audit_dom_face_moments_v45(
            mesh,geo,temperature,dirs,walls,strict_controls,sys,x,I);
        latest_energy=r.global_energy_relative;
        ++attempts;
        return dom_face_moment_qualified_v45(r,strict_controls.tolerance);
    };
    const auto strict=sys.solve(false,12,160,1.,
        preconditioned ? std::function<std::vector<double>(const std::vector<double>&)>(
            [&local](const std::vector<double>& r){return local.apply(r);}) :
            std::function<std::vector<double>(const std::vector<double>&)>(),
        face_gate);
    std::cout<<"V46_STRESS scattering="<<scattering<<" precon="<<preconditioned
        <<" original_energy="<<baseline_report.global_energy_relative
        <<" final_energy="<<latest_energy
        <<" steps="<<strict.krylov_iterations
        <<" restarts="<<strict.physical_refinements_v46
        <<" attempts="<<attempts
        <<" schur="<<strict.true_schur_relative
        <<" v42_energy="<<strict.global_energy_relative
        <<" accepted="<<strict.converged
        <<" strict_s="<<strict.seconds
        <<" setup_s="<<(preconditioned?local.setup_seconds():0.)<<'\n';
    require(strict.converged,"V46: strict physical residual gate did not converge");
    require(attempts>=1,"V46: face audit callback never executed");
    if(scattering>100.)
        require(strict.physical_refinements_v46>0,"V46: no physical refinement after rejection");
    else
        require(strict.physical_refinements_v46==0,"V46: needless refinement on moderate case");
    require(latest_energy<=ctl.tolerance,"V46: independent V45 energy threshold violated");
    require(strict.local_transport_relative<=ctl.tolerance &&
            strict.global_energy_relative<=ctl.tolerance,
            "V46: existing V42 acceptance thresholds violated");
    require(strict.scalar_mean_intensity.size()==mesh.n_cells(),
            "V46: accepted moment size mismatch");
    // Deliberately impossible physical predicate: no state can be published.
    if(preconditioned && scattering>100.) {
        std::size_t false_gate_calls=0;
        const auto reject=sys.solve(false,12,160,1.,
            [&local](const std::vector<double>& r){return local.apply(r);},
            [&false_gate_calls](const std::vector<double>&,
                                const std::vector<Cell>&) {
                ++false_gate_calls;
                return false;
            });
        require(false_gate_calls>0 && !reject.converged &&
                reject.scalar_mean_intensity.empty() &&
                reject.reflected_face_mean_intensity.empty(),
                "V46: rejected physical state was published");
        std::cout<<"V46_NEGATIVE rejected_steps="<<reject.krylov_iterations
            <<" physical_refinements="<<reject.physical_refinements_v46
            <<" callback_calls="<<false_gate_calls<<'\n';
    }
}
}
int main() {
    try {
        scenario(false,199.9);
        scenario(true,199.9);
        scenario(true,14.9);
        std::cout<<"V46_DOM_ENERGY_GATE PASS\n";
    } catch(const std::exception& e) {
        std::cerr<<"V46_DOM_ENERGY_GATE FAIL: "<<e.what()<<'\n';
        return 1;
    }
}
