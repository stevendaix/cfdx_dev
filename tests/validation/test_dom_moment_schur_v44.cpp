// V44: exact finite-volume angular Schur moment identity (including gray
// diffuse reflection), preconditioner tests, and matching production DOM.
// Does not change production defaults. Rejects unverified physical states.
#include "cfdx/physics/dom_moment_schur_v44.h"
#include "cfdx/physics/dom_quadrature_v39.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
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
void check(bool yes,const std::string& why){if(!yes)throw std::runtime_error(why);}
using CellField=Field<double,Location::CELL>;
Mesh make_rotated_anisotropic_mesh(){
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
DomWallBoundaryConditions make_walls(const Mesh& mesh,bool gray) {
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& name=mesh.boundary().patch(p).name;
        const double temp=name=="west"?540.:name=="east"?285.:365.;
        const double epsilon=!gray?1.:name=="west"?.35:name=="east"?.65:
            name=="south"?.90:name=="north"?1.:name=="bottom"?0.0:1.;
        walls[name]={epsilon,temp};
    }
    return walls;
}
void compare_assembly(DomMomentSchurSystemV44& system,
                      const Mesh& mesh,const FvGeometry& geo,
                      const CellField& T,
                      const std::vector<DiscreteDirection>& dirs,
                      const RadiationTransportControls& controls) {
    const auto n=system.dimension(),nc=mesh.n_cells();
    std::vector<double> x(n,0.);
    for(std::size_t i=0;i<n;++i) x[i]=1400.+30.*std::sin(1.13*double(i));
    std::vector<CellField> intensities;
    const auto combined=system.project(x,true,&intensities);
    const auto b=system.constant_rhs();
    const auto hom=system.project(x,false);
    const auto Hx=system.apply_exact_schur(x);
    double linearity=0.,schur_identity=0.,scale=1.;
    for(std::size_t i=0;i<n;++i) {
        linearity=std::max(linearity,std::abs(combined[i]-(b[i]+hom[i])));
        schur_identity=std::max(schur_identity,std::abs(
            (x[i]-combined[i])-(Hx[i]-b[i])));
        scale=std::max(scale,std::abs(combined[i]));
    }
    std::cout<<std::setprecision(13)<<"V44_SCHUR dimension="<<n
        <<" reflecting_faces="<<system.n_reflecting_faces()
        <<" affine_identity="<<linearity/scale
        <<" transport_schur_identity="<<schur_identity/scale<<'\n';
    check(linearity/scale<2e-11,"fixed wall/volume angular elimination not affine");
    check(schur_identity/scale<2e-11,"exact Schur identity does not match transport");
    // Independently reconstruct a complete angular FV matrix and RHS with
    // face values instead of modifying the preassembled zero-inlet RHS.
    const auto& d=dirs.front();
    Field<double,Location::FACE> flux(mesh.n_faces(),"flux","m2",1);
    ScalarBoundaryFaceValues face_values;
    ScalarBoundaryConditions bcs;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& patch=mesh.boundary().patch(p);
        bcs[patch.name]={ScalarBoundaryType::ZERO_GRADIENT,0.,0.};
        auto& vals=face_values.values[patch.name];
        vals.assign(mesh.n_faces(),std::numeric_limits<double>::quiet_NaN());
        for(const auto f:patch.face_ids) {
            const auto& info=system.wall_faces();
            const auto it=std::find_if(info.begin(),info.end(),
                [f](const auto& wall){return wall.face==f;});
            check(it!=info.end(),"missing wall face in independent assembly");
            const double reflected=it->reflect_index<system.n_reflecting_faces()
                ? (1.-it->emissivity)*(it->hplus/it->hminus)*
                    x[nc+it->reflect_index] : 0.;
            vals[f]=it->emissivity*it->ib+reflected;
        }
    }
    for(std::size_t f=0;f<mesh.n_faces();++f) {
        const auto sf=geo.face_area_vectors[f];
        flux(f)=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
    }
    CellField source(nc,"su","W/m3/sr",1),sp(nc,"sp","1/m",1),old(nc,"old","",1);
    for(std::size_t c=0;c<nc;++c) {
        source(c)=controls.absorption*blackbody_intensity(T(c))+
            controls.scattering*x[c];
        sp(c)=-(controls.absorption+controls.scattering);
        old(c)=0.;
    }
    auto eq=assemble_scalar_equation(mesh,geo,flux,0.,source,sp,bcs,true,
        &face_values,nullptr,nullptr,nullptr,ConvectionScheme::UPWIND,&old);
    Vector I(nc,0.);
    ScalarSolveControls lc;
    lc.max_iterations=controls.linear_max_iterations;
    lc.tolerance=controls.linear_tolerance;
    const auto sr=solve_scalar_equation(eq,I,lc);
    check(sr.status==SolverStatus::CONVERGED,"independent FV directional system failed");
    double independent=0.,ref=1.;
    for(std::size_t c=0;c<nc;++c){
        independent=std::max(independent,std::abs(I(c)-intensities.front()(c)));
        ref=std::max(ref,std::abs(I(c)));
    }
    std::cout<<"V44_SCHUR independent_assembled_angle="<<independent/ref<<'\n';
    check(independent/ref<2e-11,"Schur angular RHS diverges from reassembled FV");
    // Signed directions must be admitted for linear Schur actions, while NaN
    // must never silently propagate to an accepted radiative state.
    auto bad=x;bad[0]=std::numeric_limits<double>::quiet_NaN();
    bool rejected=false;
    try{(void)system.apply_exact_schur(bad);}
    catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"nonfinite Schur state not rejected");
}
void verify_case(bool gray) {
    auto mesh=make_rotated_anisotropic_mesh();
    const auto geo=build_fv_geometry(mesh);
    const auto dirs=make_dom_tensor_gauss_legendre_v39(4,16);
    auto walls=make_walls(mesh,gray);
    CellField T(mesh.n_cells(),"T","K",1);
    for(std::size_t c=0;c<mesh.n_cells();++c)
        T(c)=370.+30.*geo.cell_centres[c].x+10.*geo.cell_centres[c].y;
    RadiationTransportControls controls;
    controls.absorption=.1;controls.scattering=14.9;
    controls.max_iterations=500;controls.tolerance=1e-9;
    controls.linear_max_iterations=2000;controls.linear_tolerance=1e-12;
    DomMomentSchurSystemV44 system(mesh,geo,T,dirs,walls,controls);
    compare_assembly(system,mesh,geo,T,dirs,controls);
    const auto B=system.diffusion_matrix();
    std::vector<double> test(mesh.n_cells());
    for(std::size_t c=0;c<mesh.n_cells();++c) test[c]=std::sin(double(c)+.25);
    const auto Bt=B.matvec(test);
    double quadratic=0.;
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        quadratic+=test[c]*Bt[c];
        for(std::size_t j=0;j<mesh.n_cells();++j)
            check(std::abs(B(c,j)-B(j,c))<5e-13,"diffusion matrix not symmetric");
    }
    check(quadratic>0.,"gray-wall diffusion preconditioner not SPD");
    const auto before_T=T(0);
    const auto accelerated=system.solve(true,12,100,.02);
    std::cout<<"V44_SCHUR gray="<<gray
        <<" converged="<<accelerated.converged
        <<" krylov="<<accelerated.krylov_iterations
        <<" transport_maps="<<accelerated.transport_maps
        <<" diffusion_solves="<<accelerated.diffusion_solves
        <<" schur="<<accelerated.true_schur_relative
        <<" local="<<accelerated.local_transport_relative
        <<" energy="<<accelerated.global_energy_relative
        <<" seconds="<<accelerated.seconds<<'\n';
    check(accelerated.converged,"preconditioned exact gray moment Schur failed");
    check(T(0)==before_T,"experimental Schur mutated the temperature");
    const auto unpreconditioned=system.solve(false,12,130);
    const auto pure_diffusion=system.solve(true,12,130,1.0);
    check(pure_diffusion.converged,
        "unblended diffusion preconditioning must not bypass FV gates");
    std::cout<<"V44_SCHUR gray="<<gray
        <<" pure_diffusion_krylov="<<pure_diffusion.krylov_iterations
        <<" pure_diffusion_maps="<<pure_diffusion.transport_maps
        <<" pure_diffusion_seconds="<<pure_diffusion.seconds<<'\n';
    check(unpreconditioned.converged,"unpreconditioned exact Schur reference failed");
    std::cout<<"V44_SCHUR gray="<<gray
        <<" plain_krylov="<<unpreconditioned.krylov_iterations
        <<" plain_maps="<<unpreconditioned.transport_maps
        <<" plain_seconds="<<unpreconditioned.seconds<<'\n';
    double mismatch=0.,norm=1.;
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        mismatch=std::max(mismatch,std::abs(accelerated.scalar_mean_intensity[c]-
                                     unpreconditioned.scalar_mean_intensity[c]));
        norm=std::max(norm,std::abs(unpreconditioned.scalar_mean_intensity[c]));
    }
    check(mismatch/norm<1e-6,"preconditioning changed accepted FV state");
    CellField G(mesh.n_cells(),"G","W/m2",1),q(mesh.n_cells(),"q","W/m3",1);
    G.fill(0.);q.fill(0.);
    const auto prod=solve_participating_radiation_diffuse_gray_walls(
        mesh,geo,T,G,q,dirs,controls,walls,true,true,true);
    check(prod.converged && prod.accepted_state_verified_v42,
        "production DOM did not converge matched gray FV problem");
    double mismatch_prod=0.,norm_prod=1.;
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        mismatch_prod=std::max(mismatch_prod,std::abs(G(c)/(4.*M_PI)-
                                          accelerated.scalar_mean_intensity[c]));
        norm_prod=std::max(norm_prod,std::abs(G(c)/(4.*M_PI)));
    }
    std::cout<<"V44_SCHUR gray="<<gray<<" production_iters="<<prod.iterations
        <<" compared_to_production="<<mismatch_prod/norm_prod<<'\n';
    check(mismatch_prod/norm_prod<2e-6,"exact Schur differs from production DOM");
    const auto invalid=system.solve(true,2,1);
    check(!invalid.converged,"underbudgeted Schur must not claim success");
}
}
int main() {
    try{
        verify_case(false);verify_case(true);
        std::cout<<"V44_DOM_MOMENT_SCHUR PASS\n";
        return 0;
    }catch(const std::exception& e) {
        std::cerr<<"V44_DOM_MOMENT_SCHUR FAIL: "<<e.what()<<'\n';
        return 1;
    }
}
