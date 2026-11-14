// V47 research qualification: the continuously refreshed low-order FV face
// closure must give the EXACT high-order transport Schur defect, including
// gray-wall reflection, on rotated and strongly stretched 3D meshes.
#include "cfdx/physics/dom_holo_face_closure_v47.h"
#include "cfdx/physics/dom_quadrature_v39.h"
#include "thermal_3d_mesh_v36.h"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
using Cell=DomMomentSchurSystemV44::CellField;
void check(bool ok,const char* msg) {if(!ok)throw std::runtime_error(msg);}
Mesh grid(std::size_t nx,std::size_t ny,std::size_t nz,double z) {
    auto mesh=cfdx::testing::thermal_v36::make_mesh(nx,ny,nz,.018,0.,z);
    const double ca=std::cos(.41),sa=std::sin(.41);
    const double cb=std::cos(.33),sb=std::sin(.33);
    for(std::size_t p=0;p<mesh.points().size();++p) {
        const double x=mesh.points().x_data()[p],y=mesh.points().y_data()[p];
        const double zz=mesh.points().z_data()[p];
        const double xr=ca*x-sa*y, yr=sa*x+ca*y;
        mesh.points().set(p,cb*xr+sb*zz,yr,-sb*xr+cb*zz);
    }
    return mesh;
}
DomWallBoundaryConditions make_walls(const Mesh& mesh) {
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& name=mesh.boundary().patch(p).name;
        const double T=name=="west"?540.:name=="east"?285.:365.;
        const double eps=name=="west"?.35:name=="east"?.65:
            name=="south"?.90:name=="north"?1.:
            name=="bottom"?0.:1.;
        walls[name]={eps,T};
    }
    return walls;
}
void case_study(std::size_t nx,std::size_t ny,std::size_t nz,double z,double scattering) {
    auto mesh=grid(nx,ny,nz,z);
    const auto geo=build_fv_geometry(mesh);
    auto directions=make_dom_tensor_gauss_legendre_v39(4,8); // 32 ordinates
    const auto walls=make_walls(mesh);
    Cell temp(mesh.n_cells(),"gas","K",1);
    for(std::size_t c=0;c<mesh.n_cells();++c)
        temp(c)=370.+30.*geo.cell_centres[c].x+10.*geo.cell_centres[c].y;
    RadiationTransportControls controls;
    controls.absorption=.1;controls.scattering=scattering;
    controls.tolerance=scattering>100.?1.e-10:1.e-9;
    controls.max_iterations=500;
    controls.linear_tolerance=1.e-14;
    controls.linear_max_iterations=2000;
    DomMomentSchurSystemV44 transport(mesh,geo,temp,directions,walls,controls);
    DomHoloFaceClosureV47 holo(transport,mesh,geo,temp,directions,walls,controls);
    std::vector<double> x(transport.dimension());
    for(std::size_t c=0;c<mesh.n_cells();++c)
        x[c]=blackbody_intensity(temp(c))*(1.+.06*std::sin(double(c)*.37));
    for(const auto& w:transport.wall_faces())
        if(w.reflect_index<transport.n_reflecting_faces())
            x[mesh.n_cells()+w.reflect_index]=w.ib;
    const auto start=std::chrono::steady_clock::now();
    const auto s=holo.evaluate(x);
    check(s.max_discrete_identity_relative<2e-10,
          "V47: refreshed current/moment identity not consistent");
    check(s.max_eddington_trace_error<2e-12,
          "V47: Eddington trace must equal one");
    double face_identity=0.,internal_total=0.;
    for(std::size_t f=0;f<mesh.n_faces();++f) {
        face_identity=std::max(face_identity,std::abs(
            s.face_diffusion[f]+s.face_closure[f]-s.face_current[f]));
        if(mesh.ownership().neighbour(f)>=0)internal_total+=s.face_current[f]-s.face_current[f];
    }
    check(face_identity<1.e-8,"V47: face decomposition inconsistent");
    check(internal_total==0.,"V47: nonconservative interior face accounting");
    std::vector<double> v(x.size());
    for(std::size_t i=0;i<v.size();++i)v[i]=std::sin(.13+double(i)*.37);
    double derivative_identity=0.;
    const auto derivative=holo.linearized_low_order_defect(v,&derivative_identity);
    check(derivative_identity<5e-10,"V47: LO derivative / Schur derivative mismatch");
    const auto exact=transport.apply_exact_schur(v);
    double err=0.,denom=1.;
    for(std::size_t i=0;i<v.size();++i) {
        err=std::max(err,std::abs(derivative[i]-exact[i]));
        denom=std::max(denom,std::abs(exact[i]));
    }
    check(err/denom<5e-10,"V47: directional Schur oracle mismatch");
    // Verify that face closure is NOT frozen when the moment changes.
    const double h=.002;
    auto next=x;
    for(std::size_t i=0;i<next.size();++i)next[i]+=h*v[i];
    const auto s2=holo.evaluate(next);
    double derivative_error=0.,scale=1.,changed=0.;
    for(std::size_t i=0;i<v.size();++i) {
        derivative_error=std::max(derivative_error,std::abs(
            (s2.schur_defect[i]-s.schur_defect[i])/h-derivative[i]));
        scale=std::max(scale,std::abs(derivative[i]));
    }
    for(std::size_t f=0;f<mesh.n_faces();++f)
        changed=std::max(changed,std::abs(s2.face_closure[f]-s.face_closure[f]));
    check(derivative_error/scale<5e-7,"V47: closure not refreshed with moments");
    check(changed>1e-10,"V47: closure incorrectly kept static");
    double second_changed=0.;
    for(std::size_t k=0;k<s.eddington_symmetric.size();++k)
        second_changed=std::max(second_changed,
            std::abs(s2.eddington_symmetric[k]-s.eddington_symmetric[k]));
    check(second_changed>1e-13,
          "V47: Eddington tensor not refreshed with HO intensities");
    // Nonnegative directional intensity guarantees a positive semidefinite
    // discrete angular second moment; test sample unit vector projections.
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        const double* E=&s.eddington_symmetric[6*c];
        const double a=.4,b=-.6,cc=.7;
        const double quadratic=E[0]*a*a+E[1]*b*b+E[2]*cc*cc+
            2.*(E[3]*a*b+E[4]*a*cc+E[5]*b*cc);
        check(quadratic>=-1e-12,"V47: negative second-moment quadratic form");
    }
    // The HOLO operator cannot silently use different wall data/geometry.
    bool mismatch_rejected=false;
    try {
        auto other_walls=walls;
        (void)DomHoloFaceClosureV47(
            transport,mesh,geo,temp,directions,other_walls,controls);
    } catch(const std::invalid_argument&) {mismatch_rejected=true;}
    check(mismatch_rejected,"V47: mismatched wall object accepted");
    // Contradictory input is rejected, not clamped or accepted.
    auto invalid=x;invalid[0]=std::numeric_limits<double>::quiet_NaN();
    bool rejected=false;
    try {(void)holo.evaluate(invalid);}catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"V47: NaN HO moment was accepted");
    bool step_contract=false;
    double step_after=std::numeric_limits<double>::infinity();
    unsigned trials=0;
    if(nx==3 && scattering<100.) {
        const auto step=holo.attempt_corrected_low_order_step(x,.995,12);
        trials=step.line_search_trials;
        step_contract=step.contracted;
        step_after=step.after;
        if(step.contracted) {
            check(!step.accepted_moments.empty(),"V47: contracted state not published");
            const auto check_again=holo.evaluate(step.accepted_moments);
            check(check_again.max_fixed_point_relative < s.max_fixed_point_relative,
                  "V47: accepted correction did not reduce true HO defect");
        } else check(step.accepted_moments.empty(),"V47: failed step published state");
        const auto impossible=holo.attempt_corrected_low_order_step(x,1.e-15,3);
        check(!impossible.contracted && impossible.accepted_moments.empty(),
              "V47: rejected line-search published fields");
    }
    if(scattering<100. && nx<=8) {
        DomLocalAngularJacobiV45 local(transport,geo,directions,controls);
        // A short budget must fail closed without exposing partially updated fields.
        if(nx==3) {
            const auto short_run=holo.iterate_refreshed_closure(x,1,&local);
            check(!short_run.converged && short_run.accepted_moments.empty(),
                  "V47: underbudget HOLO published an unqualified state");
        }
        const auto adapted=holo.iterate_refreshed_closure(x,90,&local);
        if(nx==3) {
            for(const double mix : {0.15,0.50}) {
                const auto mixed=holo.iterate_refreshed_closure(x,90,&local,mix);
                std::cout<<"V47_HOLO_MIX mix="<<mix
                    <<" converged="<<mixed.converged
                    <<" accepted="<<mixed.accepted_steps
                    <<" rejected="<<mixed.rejected_steps
                    <<" last="<<mixed.final_defect
                    <<" cpu_s="<<mixed.seconds<<'\n';
                if(mixed.converged) {
                    const auto confirmed=holo.evaluate(mixed.accepted_moments);
                    check(dom_face_moment_qualified_v45(confirmed.independent,
                           controls.tolerance),"V47: mixed low order falsely qualified");
                } else check(mixed.accepted_moments.empty(),
                    "V47: mixed low-order failure published state");
            }
        }
        std::cout<<"V47_HOLO_LOCAL cells="<<mesh.n_cells()
            <<" converged="<<adapted.converged
            <<" accepted="<<adapted.accepted_steps
            <<" rejected="<<adapted.rejected_steps
            <<" first="<<adapted.initial_defect
            <<" last="<<adapted.final_defect
            <<" energy="<<adapted.physical_energy
            <<" cpu_s="<<adapted.seconds<<'\n';
        if(adapted.converged) {
            const auto check_again=holo.evaluate(adapted.accepted_moments);
            check(dom_face_moment_qualified_v45(check_again.independent,controls.tolerance),
                  "V47: local HOLO certified state fails physical FV audit");
            if(nx<=8) {
                // Independent V44 matrix-free GMRES on the SAME FV equations.
                const auto ref=transport.solve(false,12,160,1.,
                    [&local](const std::vector<double>& v){return local.apply(v);},
                    [&](const std::vector<double>& moment,
                        const std::vector<Cell>& I) {
                        const auto a=audit_dom_face_moments_v45(mesh,geo,temp,
                            directions,walls,controls,transport,moment,I);
                        return dom_face_moment_qualified_v45(a,controls.tolerance);
                    });
                check(ref.converged,"V47: independent V44 Schur reference failed");
                auto ref_m=ref.scalar_mean_intensity;
                ref_m.insert(ref_m.end(),ref.reflected_face_mean_intensity.begin(),
                             ref.reflected_face_mean_intensity.end());
                double max_delta=0.,max_scale=1.;
                for(std::size_t i=0;i<ref_m.size();++i) {
                    max_delta=std::max(max_delta,std::abs(
                        ref_m[i]-adapted.accepted_moments[i]));
                    max_scale=std::max(max_scale,std::abs(ref_m[i]));
                }
                std::cout<<"V47_HOLO_V44_REFERENCE cells="<<mesh.n_cells()
                    <<" discrepancy="<<max_delta/max_scale
                    <<" gmres_steps="<<ref.krylov_iterations
                    <<" gmres_cpu_s="<<ref.seconds<<'\n';
                check(max_delta/max_scale<2e-7,
                      "V47: refreshed HOLO state differs from exact V44 Schur");
            }
        } else check(adapted.accepted_moments.empty(),
                   "V47: failed local HOLO run published moments");
    }
    if(nx==3 && scattering>100.) {
        DomLocalAngularJacobiV45 local(transport,geo,directions,controls);
        const auto very=holo.iterate_refreshed_closure(x,250,&local);
        std::cout<<"V47_HOLO_EXTREME converged="<<very.converged
            <<" accepted="<<very.accepted_steps<<" rejected="<<very.rejected_steps
            <<" last="<<very.final_defect<<" energy="<<very.physical_energy
            <<" cpu_s="<<very.seconds<<'\n';
        check(very.converged || very.accepted_moments.empty(),
              "V47: extreme nonqualified HOLO published moments");
    }
    if(nx==3 && scattering<100.) {
        const auto run=holo.iterate_refreshed_closure(x,100);
        std::cout<<"V47_HOLO_RUN converged="<<run.converged
            <<" accepted="<<run.accepted_steps<<" rejected="<<run.rejected_steps
            <<" first="<<run.initial_defect<<" last="<<run.final_defect
            <<" energy="<<run.physical_energy<<'\n';
        if(run.converged) {
            check(run.accepted_moments.size()==x.size(),
                  "V47: qualified HOLO run has wrong dimension");
            const auto checked=holo.evaluate(run.accepted_moments);
            check(dom_face_moment_qualified_v45(checked.independent,controls.tolerance),
                  "V47: HOLO certified state fails physical FV audit");
        } else check(run.accepted_moments.empty(),
                   "V47: failed HOLO run published moments");
    }
    const double seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count();
    std::cout<<std::setprecision(13)
        <<"V47_FIXTURE cells="<<mesh.n_cells()<<" faces="<<mesh.n_faces()
        <<" nreflect="<<transport.n_reflecting_faces()
        <<" z="<<z<<" ks="<<scattering
        <<" fixed_point="<<s.max_fixed_point_relative
        <<" identity="<<s.max_discrete_identity_relative
        <<" derivative="<<derivative_identity
        <<" eddington_trace="<<s.max_eddington_trace_error
        <<" fd="<<derivative_error/scale<<" closure_change="<<changed
        <<" lo_contracted="<<step_contract<<" after="<<step_after
        <<" trials="<<trials<<" cpu_s="<<seconds<<'\n';
}
}
int main(int argc,char** argv) {
    try {
        const bool small=argc==2 && std::string(argv[1])=="--sanitizer-small";
        if(argc!=1 && !small)
            throw std::invalid_argument("V47: only --sanitizer-small supported");
        case_study(3,3,3,.03,14.9);
        if(small) { std::cout<<"V47_HOLO_SANITIZER_SMALL PASS\n";return 0; }
        case_study(4,4,6,.012,14.9);
        case_study(5,5,8,.006,14.9);
        case_study(8,8,10,.003,14.9);
        // Extreme-scattering closure identity; NOT a full HOLO convergence claim.
        case_study(3,3,3,.03,199.9);
        std::cout<<"V47_HOLO_FACE_CLOSURE PASS\n";
    } catch(const std::exception& e) {
        std::cerr<<"V47_HOLO_FACE_CLOSURE FAIL: "<<e.what()<<'\n';
        return 1;
    }
}
