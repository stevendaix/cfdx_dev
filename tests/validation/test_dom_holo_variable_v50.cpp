// V50: actual cellwise HO matrices, refreshed HO/LO RHS and independent
// gray-wall radiation balance. Experimental: no production integration.
#include "cfdx/physics/dom_holo_variable_v50.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
void req(bool b,const std::string& what){if(!b)throw std::runtime_error(what);}
Mesh mesh_3d(std::size_t nx,std::size_t ny,std::size_t nz,double depth) {
    auto mesh=cfdx::testing::thermal_v36::make_mesh(nx,ny,nz,.018,0.,depth);
    for(std::size_t i=0;i<mesh.points().size();++i) {
        const double x=mesh.points().x_data()[i];
        const double y=mesh.points().y_data()[i];
        const double z=mesh.points().z_data()[i];
        const double a=.41,b=.33;
        const double xr=std::cos(a)*x-std::sin(a)*y;
        const double yr=std::sin(a)*x+std::cos(a)*y;
        mesh.points().set(i,std::cos(b)*xr+std::sin(b)*z,
                           yr,-std::sin(b)*xr+std::cos(b)*z);
    }
    return mesh;
}
DomWallBoundaryConditions boundary(const Mesh& mesh,double base,bool nonuniform) {
    DomWallBoundaryConditions result;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& name=mesh.boundary().patch(p).name;
        result[name]={name=="bottom"?0.:name=="west"?.45:.85,
                      base+(nonuniform?(name=="east"?20.:name=="west"?-12.:0.):0.)};
    }
    return result;
}
void optical_projection_and_holo(bool short_mode) {
    auto mesh=mesh_3d(3,3,3,.03);
    auto geo=build_fv_geometry(mesh);
    const auto directions=make_dom_tensor_gauss_legendre_v39(4,16);
    auto walls=boundary(mesh,420.,false);
    using Field=DomMomentSchurSystemV44::CellField;
    Field uniform(mesh.n_cells(),"temperature","K",1);
    uniform.fill(420.);
    const auto n=mesh.n_cells();
    std::vector<double> ka(n),ks(n);
    for(std::size_t c=0;c<n;++c) {
        ka[c]=1.2+.2*std::abs(geo.cell_centres[c].x);
        ks[c]=.35+.1*std::abs(geo.cell_centres[c].z);
    }
    RadiationTransportControls settings;
    settings.absorption=1.2;settings.scattering=.35;
    settings.linear_tolerance=1e-14;
    settings.linear_max_iterations=2000;
    settings.tolerance=1e-9;
    DomMomentSchurSystemV44 uniform_ho(mesh,geo,uniform,directions,walls,settings,ka,ks);
    std::vector<double> x(uniform_ho.dimension(),blackbody_intensity(420.));
    const auto same=evaluate_dom_variable_v50(uniform_ho,x);
    req(same.true_relative<1e-9,"V50 cellwise HO isothermal invariant failed");
    req(same.global_energy_relative<1e-10,"V50 independent heterogeneous energy imbalance");
    // High-albedo, nonuniform isotropic equilibrium: spatially varying
    // scattering is assembled in each actual high-order transport matrix.
    auto thick_ka=ka,thick_ks=ks;
    for(std::size_t c=0;c<n;++c) {
        thick_ka[c]=.1+.01*std::abs(geo.cell_centres[c].z);
        thick_ks[c]=14.9+.2*std::abs(geo.cell_centres[c].x);
    }
    DomMomentSchurSystemV44 thick(mesh,geo,uniform,directions,walls,settings,
                                  thick_ka,thick_ks);
    const auto thick_state=evaluate_dom_variable_v50(thick,x);
    req(thick_state.true_relative<1e-9,
        "V50 heterogeneous high-scattering isothermal invariant failed");
    std::cout<<"V50_HO_HETEROGENEOUS_SCATTER cells="<<n
             <<" residual="<<thick_state.true_relative<<'\n';

    std::cout<<"V50_HO_HETEROGENEOUS_ISOTHERMAL cells="<<n
             <<" residual="<<same.true_relative
             <<" cell="<<same.physical_cell_relative<<" global="<<same.global_energy_relative<<'\n';
    // Independent control: equal cell coefficients must generate exactly the
    // same angular map as the historical uniform-coefficient V44 matrices.
    DomMomentSchurSystemV44 legacy(mesh,geo,uniform,directions,walls,settings);
    DomMomentSchurSystemV44 identical(mesh,geo,uniform,directions,walls,settings,
        std::vector<double>(n,settings.absorption),
        std::vector<double>(n,settings.scattering));
    const auto ref=legacy.project(x,true);
    const auto same_op=identical.project(x,true);
    double difference=0.;
    for(std::size_t i=0;i<ref.size();++i)
        difference=std::max(difference,std::abs(ref[i]-same_op[i]));
    req(difference<1e-11,"V50 cellwise assembly changed uniform HO transport");
    bool threw=false;
    auto negative=ks;negative[0]=-1.;
    try{DomMomentSchurSystemV44 bad(mesh,geo,uniform,directions,walls,settings,ka,negative);}
    catch(const std::invalid_argument&){threw=true;}
    req(threw,"V50 negative scattering passed validation");
    threw=false;
    try{DomMomentSchurSystemV44 bad(mesh,geo,uniform,directions,walls,settings,ka,{});}
    catch(const std::invalid_argument&){threw=true;}
    req(threw,"V50 unmatched cell coefficient fields passed");
    threw=false;
    try{(void)uniform_ho.solve(false);}
    catch(const std::invalid_argument&){threw=true;}
    req(threw,"V50 variable HO silently used legacy uniform physical oracle");
    // Actual inhomogeneous temperature, scattering, absorption and gray walls.
    walls=boundary(mesh,420.,true);
    Field temperature(n,"temperature","K",1);
    for(std::size_t c=0;c<n;++c)temperature(c)=420.+5.*geo.cell_centres[c].x;
    DomMomentSchurSystemV44 ho(mesh,geo,temperature,directions,walls,settings,ka,ks);
    x.assign(ho.dimension(),blackbody_intensity(420.));
    for(const auto& wall:ho.wall_faces())
        if(wall.reflect_index<ho.n_reflecting_faces())
            x[n+wall.reflect_index]=wall.ib;
    DomV50Options o;
    o.max_outer=short_mode?6:65;
    o.physical_tolerance=1e-9;
    o.linear_tolerance=1e-12;
    o.preconditioner=DomSmmV49Preconditioner::JacobiAuxiliary;
    // HO/LO identity is refreshed for the current real angular face currents;
    // perturbing the state must also perturb the computed face correction.
    const auto snap=evaluate_dom_variable_v50(ho,x);
    std::vector<DomSmmSecondMomentV48::Tensor> tensor=snap.eddington;
    std::vector<double> kt(n);
    for(std::size_t c=0;c<n;++c)kt[c]=ka[c]+ks[c];
    DomSmmGrayWallV49 low(mesh,geo,directions,walls,tensor,kt,ka);
    std::vector<double> low_guess(x.begin(),x.begin()+static_cast<std::ptrdiff_t>(n));
    const auto identity=low.consistent_source(low_guess,snap.gas_ib,snap.face_current);
    req(identity.max_discrete_identity_relative<1e-11,
        "V50 failed HO/LO discrete energy identity");
    auto changed=x;changed[0]*=1.015;
    const auto changed_snap=evaluate_dom_variable_v50(ho,changed);
    DomSmmGrayWallV49 refreshed(mesh,geo,directions,walls,
        changed_snap.eddington,kt,ka);
    std::vector<double> changed_guess(changed.begin(),changed.begin()+static_cast<std::ptrdiff_t>(n));
    const auto new_source=refreshed.consistent_source(changed_guess,
                           changed_snap.gas_ib,changed_snap.face_current);
    req(new_source.max_discrete_identity_relative<1e-11,
        "V50 failed refreshed HO/LO discrete energy identity");
    double source_change=0.;
    for(std::size_t f=0;f<mesh.n_faces();++f)
        source_change=std::max(source_change,
            std::abs(identity.face_ho_minus_lo_watts[f]-
                     new_source.face_ho_minus_lo_watts[f]));
    req(source_change>1e-10,"V50 stale HO/LO face correction");
    std::cout<<"V50_HO_LO_IDENTITY before="<<identity.max_discrete_identity_relative
             <<" after="<<new_source.max_discrete_identity_relative
             <<" correction_change="<<source_change<<'\n';
    const auto solution=solve_dom_holo_v50(ho,x,o);
    if(!short_mode) {
        req(solution.converged,"V50 refreshed HOLO did not close: "+solution.diagnostic);
        req(solution.final_true<=o.physical_tolerance,"V50 false physical acceptance");
        req(solution.accepted_steps>0,"V50 HOLO did not update its state");
        req(solution.lo_solves>0,"V50 HOLO did not use LO correction");
        req(solution.lo_accepted_steps>0,
            "V50 HOLO never accepted a refreshed low-order step");
        req(solution.final_moments.size()==ho.dimension(),"V50 published incomplete state");
        for(std::size_t i=1;i<solution.accepted_defect_history.size();++i)
            req(solution.accepted_defect_history[i]<solution.accepted_defect_history[i-1],
                "V50 nonmonotone physical acceptance");
    }
    // Impossible physical gate must fail closed (no fields published), even
    // when some intermediate accepted iterates have contracted.
    DomV50Options impossible=o;
    impossible.physical_tolerance=1e-30;
    impossible.max_outer=2;
    const auto rejected=solve_dom_holo_v50(ho,x,impossible);
    req(!rejected.converged && rejected.final_moments.empty(),
        "V50 impossible physical gate published a solution");
    std::cout<<"V50_HOLO_VARIABLE cells="<<n
             <<" success="<<solution.converged
             <<" initial="<<solution.initial_true
             <<" final="<<solution.final_true
             <<" iterations="<<solution.accepted_steps
             <<" lo="<<solution.lo_solves
             <<" lo_accepted="<<solution.lo_accepted_steps
             <<" fallbacks="<<solution.picard_fallbacks
             <<" rejected="<<solution.rejected_trials
             <<" negative_gate="<<rejected.converged<<'\n';
}
void amg_stress(std::size_t nx,std::size_t ny,std::size_t nz,double depth) {
    auto mesh=mesh_3d(nx,ny,nz,depth);
    const auto geo=build_fv_geometry(mesh);
    const auto dirs=make_dom_tensor_gauss_legendre_v39(4,16);
    std::vector<DomSmmSecondMomentV48::Tensor> E(mesh.n_cells(),
        {1./3.,1./3.,1./3.,.014,.006,.007});
    std::vector<double> ka(mesh.n_cells()),kt(mesh.n_cells());
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        ka[c]=.3+.02*std::abs(geo.cell_centres[c].x);
        kt[c]=3.+.1*std::abs(geo.cell_centres[c].x);
    }
    DomSmmGrayWallV49 lo(mesh,geo,dirs,boundary(mesh,420.,true),E,kt,ka);
    const auto aux=lo.auxiliary_matrix();
    req(aux.n_rows()==mesh.n_cells() && aux.nnz()>mesh.n_cells(),
        "V50 large-grid auxiliary matrix incomplete");
    const std::vector<double> probe(mesh.n_cells(),blackbody_intensity(420.));
    const auto action=lo.apply(probe);
    for(double v:action)req(std::isfinite(v),"V50 large-grid true action nonfinite");
    // Structural large-grid AMG certificate: genuine setup plus one V-cycle
    // on the actual matrix-free second-moment operator. This does NOT certify
    // convergence of the full 640-cell radiation solve.
    const auto nc=mesh.n_cells();
    cfdx::core::LinearOperator op;
    op.size=nc;
    op.apply=[&lo,nc](const cfdx::core::Vector& v,cfdx::core::Vector& y){
        std::vector<double> a(nc);
        for(std::size_t c=0;c<nc;++c)a[c]=v(c);
        const auto b=lo.apply(a);
        if(y.size()!=nc)y.resize(nc);
        for(std::size_t c=0;c<nc;++c)y(c)=b[c];
    };
    cfdx::core::FunctionalLinearOperator actual(nc,op.apply);
    cfdx::core::MatrixFreeVcyclePreconditioner amg(actual,.55,2,2,.25,12,
        cfdx::core::AMGInterpolationPolicy::SmoothedAggregation,false,false,
        cfdx::core::AMGSmootherPolicy::DampedJacobi);
    const auto t0=std::chrono::steady_clock::now();
    req(amg.setup(aux),"V50 AMG large-grid setup rejected: "+amg.last_setup_error());
    cfdx::core::Vector r(nc,0.),z(nc,0.);
    for(std::size_t c=0;c<nc;++c)r(c)=std::sin(.03*static_cast<double>(c));
    req(amg.apply(r,z),"V50 AMG large-grid V-cycle application failed");
    double size=0.;
    for(std::size_t c=0;c<nc;++c) {
        req(std::isfinite(z(c)),"V50 AMG large-grid nonfinite preconditioned vector");
        size=std::max(size,std::abs(z(c)));
    }
    req(size>0.,"V50 AMG large-grid produced a zero correction");
    const double seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-t0).count();
    std::cout<<"V50_AMG_STRESS cells="<<mesh.n_cells()
             <<" graph_nnz="<<aux.nnz()
             <<" levels="<<amg.hierarchy_level_sizes().size()
             <<" setup_plus_one_cycle_seconds="<<seconds
             <<" status=SETUP_AND_ONE_CYCLE_ONLY_NOT_SOLVER_CERTIFIED\n";
}
}
int main(int argc,char** argv) {
    try {
        const std::string mode=argc==1?"--quick":argc==2?argv[1]:"invalid";
        if(mode=="--quick") optical_projection_and_holo(false);
        else if(mode=="--sanitizer-small") optical_projection_and_holo(true);
        else if(mode=="--stress-640")amg_stress(8,8,10,.003);
        else if(mode=="--stress-10000")amg_stress(20,20,25,.003);
        else throw std::invalid_argument("V50: mode --quick/--sanitizer-small/--stress-640/--stress-10000");
        std::cout<<"V50_ALL_PASS\n";
    } catch(const std::exception& e) {
        std::cerr<<"V50_FAIL: "<<e.what()<<'\n';return 1;
    }
}
