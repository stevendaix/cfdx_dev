// V43 research-only FV-compatible DSA candidate: compares TRUE transport
// fixed-point residuals under the same black-wall optical problem. Never
// changes the production DOM selection or acceptance default.
#include "cfdx/physics/dom_dsa_study_v43.h"
#include "cfdx/physics/dom_quadrature_v39.h"
#include "thermal_3d_mesh_v36.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
void check(bool yes,const char* why) {if(!yes)throw std::runtime_error(why);}
void run() {
    auto mesh=cfdx::testing::thermal_v36::make_mesh(3,3,3,.023,0.,1.0);
    // Rotate the entire physical enclosure: walls, not just cells, are oblique.
    const double ca=std::cos(.41),sa=std::sin(.41);
    const double cb=std::cos(.33),sb=std::sin(.33);
    for(std::size_t p=0;p<mesh.points().size();++p) {
        const double x=mesh.points().x_data()[p];
        const double y=mesh.points().y_data()[p];
        const double z=mesh.points().z_data()[p];
        const double xr=ca*x-sa*y,yr=sa*x+ca*y;
        mesh.points().set(p,cb*xr+sb*z,yr,-sb*xr+cb*z);
    }
    const auto geo=build_fv_geometry(mesh);
    auto dirs=make_dom_tensor_gauss_legendre_v39(4,16);
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& name=mesh.boundary().patch(p).name;
        walls[name]={1.0,name=="west"?530.0:name=="east"?280.0:370.0};
    }
    Field<double,Location::CELL> T(mesh.n_cells(),"T","K",1);
    for(std::size_t c=0;c<mesh.n_cells();++c)
        T(c)=360.+30.*geo.cell_centres[c].x;
    RadiationTransportControls controls;
    controls.absorption=.10;controls.scattering=14.9;
    controls.tolerance=1e-10;controls.max_iterations=400;
    controls.linear_tolerance=1e-12;
    // The diffusive correction operator must be symmetric and strictly SPD.
    const auto B=build_dom_dsa_diffusion_matrix_v43(
        mesh,geo,controls.absorption,controls.scattering);
    std::vector<double> test(mesh.n_cells(),0.0);
    for(std::size_t c=0;c<mesh.n_cells();++c) test[c]=std::sin(double(c)+1.);
    const auto Bt=B.matvec(test);
    long double quadratic=0.;
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        quadratic+=test[c]*Bt[c];
        for(std::size_t j=0;j<mesh.n_cells();++j)
            check(std::abs(B(c,j)-B(j,c))<1e-14,
                "low-order diffusion matrix not symmetric");
    }
    check(quadratic>0.,"low-order diffusion matrix not positive definite");
    const auto baseline=study_dom_dsa_v43(mesh,geo,T,dirs,walls,controls,false);
    const auto candidate=study_dom_dsa_v43(mesh,geo,T,dirs,walls,controls,true);
    std::cout<<"V43_DSA baseline_converged="<<baseline.converged
             <<" baseline_iters="<<baseline.iterations
             <<" baseline_maps="<<baseline.transport_maps
             <<" baseline_seconds="<<baseline.elapsed_s
             <<" baseline_fixed_point="<<baseline.fixed_point_relative
             <<" baseline_transport="<<baseline.local_transport_relative
             <<" baseline_energy="<<baseline.global_energy_relative
             <<" accelerated_converged="<<candidate.converged
             <<" accelerated_iters="<<candidate.iterations
             <<" accelerated_maps="<<candidate.transport_maps
             <<" accelerated_seconds="<<candidate.elapsed_s
             <<" accepted="<<candidate.accepted_accelerations
             <<" rejected="<<candidate.rejected_accelerations
             <<" fixed_point="<<candidate.fixed_point_relative
             <<" transport="<<candidate.local_transport_relative
             <<" energy="<<candidate.global_energy_relative<<'\n';
    check(baseline.converged && candidate.converged,
        "baseline or prototype DSA did not reach physical FV gate");
    check(candidate.accepted_accelerations>0,
        "DSA trial never improved the actual transport residual");
    double diff=0.,scale=1.;
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        diff=std::max(diff,std::abs(baseline.scalar_mean_intensity[c]-
                                   candidate.scalar_mean_intensity[c]));
        scale=std::max(scale,std::abs(baseline.scalar_mean_intensity[c]));
    }
    std::cout<<"V43_DSA relative_scalar_difference="<<diff/scale<<'\n';
    check(diff/scale<2e-6,"DSA changed the accepted radiative fixed point");
    // Cross-check the prototype against the original full production 3D DOM,
    // including its independently audited final angular FV state.
    Field<double,Location::CELL> G(mesh.n_cells(),"G","W/m2",1);
    Field<double,Location::CELL> Q(mesh.n_cells(),"q","W/m3",1);
    G.fill(0.);Q.fill(0.);
    const auto production=solve_participating_radiation_diffuse_gray_walls(
        mesh,geo,T,G,Q,dirs,controls,walls,true,true,true);
    check(production.converged && production.accepted_state_verified_v42,
        "production DOM failed matching high-scattering black-wall case");
    double prod_diff=0.,prod_scale=1.;
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        prod_diff=std::max(prod_diff,
            std::abs(G(c)/(4.*M_PI)-candidate.scalar_mean_intensity[c]));
        prod_scale=std::max(prod_scale,std::abs(G(c)/(4.*M_PI)));
    }
    std::cout<<"V43_DSA production_iters="<<production.iterations
             <<" production_relative_difference="<<prod_diff/prod_scale<<'\n';
    check(prod_diff/prod_scale<2e-6,"DSA candidate differs from production DOM");
    // Outside the diffusion-dominated regime acceleration must not alter
    // the fixed point or bypass the accepted-state FV residual gate.
    auto thin=controls;
    thin.absorption=.4;thin.scattering=.1;thin.max_iterations=80;
    const auto thin_baseline=study_dom_dsa_v43(mesh,geo,T,dirs,walls,thin,false);
    const auto thin_dsa=study_dom_dsa_v43(mesh,geo,T,dirs,walls,thin,true);
    check(thin_baseline.converged && thin_dsa.converged,
        "thin-transport control problem did not converge");
    double thin_diff=0.,thin_scale=1.;
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        thin_diff=std::max(thin_diff,std::abs(
            thin_baseline.scalar_mean_intensity[c]-thin_dsa.scalar_mean_intensity[c]));
        thin_scale=std::max(thin_scale,std::abs(thin_baseline.scalar_mean_intensity[c]));
    }
    check(thin_diff/thin_scale<2e-6,
        "DSA candidate changed optically thinner accepted solution");
    std::cout<<"V43_DSA thinner_baseline_iters="<<thin_baseline.iterations
             <<" thinner_candidate_iters="<<thin_dsa.iterations
             <<" thinner_relative_difference="<<thin_diff/thin_scale<<'\n';
    // Fail closed beyond the prototype scope.
    walls.begin()->second.emissivity=.9;
    bool rejected=false;
    try{(void)study_dom_dsa_v43(mesh,geo,T,dirs,walls,controls,true);}
    catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"prototype must explicitly reject gray diffuse reflection");
}
} // namespace
int main(){try{run();std::cout<<"V43_DSA_STUDY PASS\n";return 0;}
catch(const std::exception& e){std::cerr<<"V43_DSA_STUDY FAIL: "<<e.what()<<'\n';return 1;}}
