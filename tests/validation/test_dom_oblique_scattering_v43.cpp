// V43 physical quadrature and nonlinear DOM checks on genuinely OBLIQUE
// exterior walls, not only on skew interior cells. All kernels are production.
#include "cfdx/physics/dom_quadrature_v39.h"
#include "cfdx/physics/radiation_solver.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
void check(bool yes,const char* why) { if (!yes) throw std::runtime_error(why); }
using CellField = Field<double,Location::CELL>;
Mesh rotated_oblique_mesh() {
    auto mesh=cfdx::testing::thermal_v36::make_mesh(3,3,3,.026,0.,1.);
    const double ca=std::cos(.41),sa=std::sin(.41),cb=std::cos(.33),sb=std::sin(.33);
    for (std::size_t p=0;p<mesh.points().size();++p) {
        const double x=mesh.points().x_data()[p],y=mesh.points().y_data()[p];
        const double z=mesh.points().z_data()[p];
        const double xr=ca*x-sa*y, yr=sa*x+ca*y;
        mesh.points().set(p,cb*xr+sb*z,yr,-sb*xr+cb*z);
    }
    return mesh;
}
DomWallBoundaryConditions walls_for(const Mesh& mesh,bool isothermal) {
    DomWallBoundaryConditions walls;
    for (std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& name=mesh.boundary().patch(p).name;
        const double T=isothermal ? 370.0 :
            name=="west" ? 530.0 : name=="east" ? 310.0 :
            name=="north" ? 420.0 : 360.0;
        walls[name]={0.65,T};
    }
    return walls;
}
struct Case { std::vector<double> G; double moment_error=0.; double flux_norm=0.;
    double energy_error=0.; std::size_t iters=0; };
Case solve(const Mesh& mesh,const FvGeometry& geo,std::size_t nmu,std::size_t nphi) {
    const auto directions=make_dom_tensor_gauss_legendre_v39(nmu,nphi);
    const auto moments=audit_dom_wall_quadrature(mesh,geo,directions);
    const auto walls=walls_for(mesh,false);
    CellField temp(mesh.n_cells(),"T","K",1),G(mesh.n_cells(),"G","W/m2",1);
    CellField source(mesh.n_cells(),"Q","W/m3",1);
    for (std::size_t c=0;c<mesh.n_cells();++c)
        temp(c)=345.+70.*geo.cell_centres[c].x+30.*geo.cell_centres[c].y;
    G.fill(4.*blackbody_emissive_power(370.));source.fill(0.);
    RadiationTransportControls ctrl;
    ctrl.absorption=.20;ctrl.scattering=5.0; // albedo = 0.96153846
    ctrl.max_iterations=300;ctrl.tolerance=1e-9;ctrl.linear_tolerance=1e-12;
    const auto r=solve_participating_radiation_diffuse_gray_walls(
        mesh,geo,temp,G,source,directions,ctrl,walls,true,true,true);
    check(r.converged && r.accepted_state_verified_v42,
        "strict DOM scattering failed on rotated 3D walls");
    check(r.max_cell_transport_residual_v42<ctrl.tolerance,
        "strict DOM cell transport residual not satisfied");
    check(r.relative_energy_closure_v42<ctrl.tolerance,
        "strict DOM energy closure not satisfied");
    double net_wall=0.;for (auto& kv:r.wall_outward_W) net_wall+=kv.second;
    Case out;
    out.G.reserve(mesh.n_cells());
    for (std::size_t c=0;c<mesh.n_cells();++c) out.G.push_back(G(c));
    out.moment_error=moments.maximum_relative_pi_error;
    out.energy_error=r.relative_energy_closure_v42;
    out.flux_norm=net_wall;
    out.iters=r.iterations;
    std::cout<<std::setprecision(12)<<"V43_OBLIQUE angles="<<directions.size()
        <<" max_hemi_pi="<<out.moment_error<<" iterations="<<out.iters
        <<" energy="<<out.energy_error<<" net_wall_W="<<out.flux_norm<<'\n';
    return out;
}
double difference(const Case& a,const Case& fine) {
    check(a.G.size()==fine.G.size(),"inconsistent DOM refinement grids");
    long double sum=0.,scale=0.;
    for(std::size_t i=0;i<a.G.size();++i) {
        sum+=static_cast<long double>(a.G[i]-fine.G[i])*(a.G[i]-fine.G[i]);
        scale+=static_cast<long double>(fine.G[i])*fine.G[i];
    }
    return std::sqrt(double(sum/scale));
}
void analytic_oracles(const Mesh& mesh,const FvGeometry& geo) {
    const auto dirs=make_dom_tensor_gauss_legendre_v39(8,32);
    const auto walls=walls_for(mesh,true);
    CellField T(mesh.n_cells(),"T","K",1);T.fill(370.0);
    std::vector<CellField> I;
    const double Ib=blackbody_intensity(370.);
    for (std::size_t m=0;m<dirs.size();++m) {
        I.emplace_back(mesh.n_cells(),"I","W/m2/sr",1);I.back().fill(Ib);
    }
    RadiationTransportControls controls;controls.absorption=.2;controls.scattering=50.;
    const auto r=audit_dom_accepted_state_v42(mesh,geo,T,I,dirs,controls,walls);
    check(r.max_cell_relative<1e-12 && r.relative_global_energy<1e-11,
        "isothermal oblique-wall exact solution must satisfy FV scattering operator");
    std::size_t oblique=0;
    double maximum_pi=0.;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p)
        for (std::size_t face:mesh.boundary().patch(p).face_ids) {
            const auto sf=geo.face_area_vectors[face];
            const double area=sf.mag();
            const std::array<double,3> n={sf.x/area,sf.y/area,sf.z/area};
            if (std::abs(n[0])>.10 && std::abs(n[1])>.10 && std::abs(n[2])>.10)
                ++oblique;
            const auto h=dom_hemisphere_moments(dirs,n);
            maximum_pi=std::max(maximum_pi,h.relative_pi_error);
            // Analytically: integrate (n.s) [1 + a (n.s)] over the outgoing hemisphere.
            // Integral = pi + a*(2*pi/3); unlike isothermal equilibrium,
            // this exposes inaccurate physical wall flux even when the
            // discrete wall treatment has perfect zero energy residual.
            const double a=.25;
            long double numerical=0.;
            for(const auto& d:dirs) {
                const double mu=d.dx*n[0]+d.dy*n[1]+d.dz*n[2];
                if (mu>0.) numerical+=d.weight*mu*(1.+a*mu);
            }
            const double expected=M_PI*(1.+2.*a/3.);
            check(std::abs(double(numerical)/expected-1.)<.035,
                "oblique analytic angular-flux moment outside 3.5% criterion");
        }
    check(oblique>0,"the new test does not actually include oblique exterior walls");
    check(maximum_pi<.035,"refined quadrature exceeds physical oblique wall moment gate");
    std::cout<<"V43_EXACT_ISOTHERMAL max_fv="<<r.max_cell_relative
             <<" oblique_boundary_faces="<<oblique<<" pi_error="<<maximum_pi<<'\n';
}
} // namespace
int main() {
    try {
        auto mesh=rotated_oblique_mesh();auto geo=build_fv_geometry(mesh);
        analytic_oracles(mesh,geo);
        const auto q64=solve(mesh,geo,4,16);
        const auto q128=solve(mesh,geo,8,16);
        const auto q256=solve(mesh,geo,8,32);
        const auto e64=difference(q64,q256),e128=difference(q128,q256);
        std::cout<<"V43_REFINEMENT relative_L2_64="<<e64
                 <<" relative_L2_128="<<e128<<'\n';
        check(q256.moment_error<q64.moment_error,
            "angular refinement did not reduce oblique physical wall quadrature error");
        check(e128<e64,"radiation angular refinement did not contract against the fine ordinate reference");
        std::cout<<"V43_OBLIQUE_SCATTERING PASS\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"V43_OBLIQUE_SCATTERING FAIL: "<<e.what()<<'\n';
        return 1;
    }
}
