// V49: autonomous discrete half-range P1/gray-wall LO operator, variable
// optical coefficients, AMG auxiliary preconditioner and fail-closed gates.
// This does NOT certify a complete high-order/low-order DSA or production DOM.
#include "cfdx/physics/dom_smm_gray_wall_v49.h"
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
void require(bool condition,const std::string& reason) {
    if(!condition)throw std::runtime_error(reason);
}
Mesh grid(std::size_t nx,std::size_t ny,std::size_t nz,double depth) {
    auto mesh=cfdx::testing::thermal_v36::make_mesh(nx,ny,nz,.018,0.,depth);
    constexpr double a=.41,b=.33;
    for(std::size_t i=0;i<mesh.points().size();++i) {
        const double x=mesh.points().x_data()[i],y=mesh.points().y_data()[i];
        const double z=mesh.points().z_data()[i];
        const double xr=std::cos(a)*x-std::sin(a)*y;
        const double yr=std::sin(a)*x+std::cos(a)*y;
        mesh.points().set(i,std::cos(b)*xr+std::sin(b)*z,
                           yr,-std::sin(b)*xr+std::cos(b)*z);
    }
    return mesh;
}
DomWallBoundaryConditions walls(const Mesh& mesh,double t,bool nonuniform) {
    DomWallBoundaryConditions w;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& name=mesh.boundary().patch(p).name;
        w[name]={name=="bottom"?0.0:name=="west"?.45:.85,
                 !nonuniform?t:name=="east"?t+95.:name=="west"?t-35.:t};
    }
    return w;
}
using Tensor=DomSmmSecondMomentV48::Tensor;
std::vector<Tensor> tensors(const FvGeometry& geo,bool variable) {
    std::vector<Tensor> E(geo.cell_volumes.size());
    for(std::size_t c=0;c<E.size();++c) {
        const double x=variable?geo.cell_centres[c].x:0.;
        E[c]={1./3.+.009*x,1./3.-.009*x,1./3.,.014,.006,.007};
    }
    return E;
}
double ib(double temperature){return blackbody_intensity(temperature);}
void test_case(std::size_t nx,std::size_t ny,std::size_t nz,double depth,
               bool sanitizer_small) {
    auto mesh=grid(nx,ny,nz,depth);
    const auto geo=build_fv_geometry(mesh);
    const auto directions=make_dom_tensor_gauss_legendre_v39(4,16);
    const auto E=tensors(geo,true);
    const auto nc=mesh.n_cells();
    std::vector<double> ka(nc),kt(nc),bb(nc,ib(420.));
    for(std::size_t c=0;c<nc;++c) {
        ka[c]=.3+.02*std::abs(geo.cell_centres[c].x);
        kt[c]=3.0+.10*std::abs(geo.cell_centres[c].x);
    }
    const std::vector<Tensor> Eiso(nc,{1./3.,1./3.,1./3.,0.,0.,0.});
    DomSmmGrayWallV49 same(mesh,geo,directions,walls(mesh,420.,false),Eiso,kt,ka);
    const auto uniform=same.solve(bb,DomSmmV49Preconditioner::None,1.e-10,180);
    if(!uniform.accepted) std::cout<<"UNIFORM_DIAG linear="<<uniform.true_linear_relative<<" cell="<<uniform.max_cell_balance_relative<<" global="<<uniform.global_energy_relative<<" iter="<<uniform.iterations<<" epsilon="<<uniform.maximum_quadrature_pi_error<<std::endl;
    require(uniform.accepted,"V49: heterogeneous isothermal equilibrium failed: "+uniform.diagnostic);
    double uniform_error=0.;
    for(const auto value:uniform.accepted_J)
        uniform_error=std::max(uniform_error,std::abs(value-bb[0])/std::max(1.,bb[0]));
    require(uniform_error<1.e-9,"V49: isothermal gray equilibrium biased");
    double max_q=0.;
    for(const auto& w:same.walls())
        max_q=std::max(max_q,std::abs(uniform.accepted_wall_flux_watts[w.face]));
    require(max_q<1.e-9,"V49: nonzero wall flux at equilibrium");
    // Actual matrix-free LO solve with heterogeneous E, kt, ka and nonuniform
    // wall emissivities/temperatures. No analytic HO physical claim here.
    DomSmmGrayWallV49 system(mesh,geo,directions,walls(mesh,420.,true),E,kt,ka);
    for(std::size_t c=0;c<nc;++c)
        bb[c]=ib(420.+10.*geo.cell_centres[c].x);
    if(nc>96) {
        // Stress the tensor matrix-free action and the SPD auxiliary graph,
        // without treating an unbounded full solve as a CI regression.
        const auto t0=std::chrono::steady_clock::now();
        const auto y=system.apply(bb);
        const auto A=system.auxiliary_matrix();
        require(A.n_rows()==nc && A.nnz()>nc,"V49: large auxiliary graph missing");
        for(const double value:y)
            require(std::isfinite(value),"V49: large anisotropic action nonfinite");
        const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
        std::cout<<"V49_LARGE_OPERATOR_ONLY cells="<<nc
                 <<" nnz="<<A.nnz()<<" seconds="<<elapsed
                 <<" status=NO_LARGE_GMRES_CERTIFICATE"<<std::endl;
        return;
    }
    const auto direct=system.solve(bb,DomSmmV49Preconditioner::None,1.e-10,260);
    if(!direct.accepted) std::cout<<"DIRECT_DIAG cells="<<nc<<" linear="<<direct.true_linear_relative<<" cell="<<direct.max_cell_balance_relative<<" global="<<direct.global_energy_relative<<" iter="<<direct.iterations<<std::endl;
    require(direct.accepted,"V49: no-preconditioner gray solve failed: "+direct.diagnostic);
    const auto pre=system.solve(bb,DomSmmV49Preconditioner::JacobiAuxiliary,1.e-10,260);
    require(pre.accepted,"V49: Jacobi gray solve failed: "+pre.diagnostic);
    double spread=0.;
    for(std::size_t c=0;c<nc;++c)
        spread=std::max(spread,std::abs(direct.accepted_J[c]-pre.accepted_J[c])/
                               std::max({1.,std::abs(direct.accepted_J[c]),
                                            std::abs(pre.accepted_J[c])}));
    require(spread<5.e-8,"V49: true states differ between preconditioners");
    // Physical current must have owner-only sign at the boundary; a perfect
    // reflector (epsilon=0) has identically zero power without a flux repair.
    for(const auto& w:system.walls())
        if(w.emissivity==0.) {
            require(pre.accepted_wall_flux_watts[w.face]==0.,
                    "V49: perfect gray reflector transfers power");
        }
    // Test that the auxiliary matrix is an actual symmetric positive graph.
    const auto auxiliary=system.auxiliary_matrix();
    for(std::size_t c=0;c<nc;++c) {
        require(auxiliary(c,c)>0.,"V49: nonpositive auxiliary diagonal");
        for(std::size_t j=0;j<nc;++j) if(c!=j && auxiliary(c,j)!=0.) {
            require(std::abs(auxiliary(c,j)-auxiliary(j,c))<1.e-12,
                    "V49: nonsymmetric auxiliary graph");
            require(auxiliary(c,j)<0.,"V49: non-M-matrix auxiliary graph");
        }
    }
    if(nc<=96 && !sanitizer_small) {
        const auto amg=system.solve(bb,DomSmmV49Preconditioner::AmgAuxiliary,1.e-10,260);
        if(!amg.accepted)
            std::cout<<"V49_AMG_NOT_QUALIFIED cells="<<nc<<" why="<<amg.diagnostic<<'\n';
        else {
            double amg_spread=0.;
            for(std::size_t c=0;c<nc;++c)
                amg_spread=std::max(amg_spread,std::abs(amg.accepted_J[c]-direct.accepted_J[c])/
                              std::max(1.,std::abs(direct.accepted_J[c])));
            require(amg_spread<5e-8,"V49: AMG accepted a different physical state");
            std::cout<<"V49_AMG_OK cells="<<nc<<" levels="<<amg.preconditioner_levels
                     <<" iterations="<<amg.iterations<<" seconds="<<amg.seconds<<'\n';
        }
    }
    // A corrupted moment is never returned as accepted; invalid coefficient
    // and quadrature inputs must be rejected before a Krylov solve.
    bool refused=false;
    auto invalid=E;invalid[0]={1.,0.,0.,0.,0.,0.};
    try {DomSmmGrayWallV49 bad(mesh,geo,directions,walls(mesh,420.,true),invalid,kt,ka);}
    catch(const std::invalid_argument&){refused=true;}
    require(refused,"V49: nonpositive E tensor accepted");
    refused=false;
    auto invalid_kt=kt;invalid_kt[0]=0.;
    try {DomSmmGrayWallV49 bad(mesh,geo,directions,walls(mesh,420.,true),E,invalid_kt,ka);}
    catch(const std::invalid_argument&){refused=true;}
    require(refused,"V49: zero optical extinction accepted");
    refused=false;
    try {DomSmmGrayWallV49 bad(mesh,geo,directions,walls(mesh,420.,true),E,kt,ka,0.);}
    catch(const std::invalid_argument&){refused=true;}
    require(refused,"V49: failed to enforce strict hemispheric quadrature gate");
    // Never publish a truncated Krylov iterate as physically accepted.
    const auto rejected=system.solve(bb,DomSmmV49Preconditioner::None,1.e-14,1,3);
    require(!rejected.accepted && rejected.accepted_J.empty(),
            "V49: insufficient Krylov budget published a candidate");
    std::cout<<"V49_GRAY_LO cells="<<nc<<" no_prec_iter="<<direct.iterations
             <<" jacobi_iter="<<pre.iterations
             <<" solution_spread="<<spread
             <<" true_linear="<<pre.true_linear_relative
             <<" local_balance="<<pre.max_cell_balance_relative
             <<" global_balance="<<pre.global_energy_relative
             <<" hemisphere_error="<<system.max_pi_error()
             <<" serial_seconds="<<direct.seconds
             <<" jacobi_seconds="<<pre.seconds<<'\n';
}
void true_dom_gray_source() {
    // Actual UPWIND transport on the same FV mesh: the HO face currents and
    // anisotropic E come from angular intensities, not manufactured closure.
    auto mesh=grid(3,3,3,.03);
    const auto geo=build_fv_geometry(mesh);
    const auto dirs=make_dom_tensor_gauss_legendre_v39(4,16);
    auto bc=walls(mesh,410.,true);
    using Cell=DomMomentSchurSystemV44::CellField;
    Cell temp(mesh.n_cells(),"temperature","K",1);
    std::vector<double> gas_bb(mesh.n_cells());
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        temp(c)=390.+25.*geo.cell_centres[c].x;
        gas_bb[c]=ib(temp(c));
    }
    RadiationTransportControls controls;
    controls.absorption=.1;controls.scattering=14.9;
    controls.tolerance=1e-9;controls.linear_tolerance=1e-14;
    controls.linear_max_iterations=2000;
    DomMomentSchurSystemV44 transport(mesh,geo,temp,dirs,bc,controls);
    DomHoloFaceClosureV47 holo(transport,mesh,geo,temp,dirs,bc,controls);
    std::vector<double> x(transport.dimension());
    for(std::size_t c=0;c<mesh.n_cells();++c)x[c]=gas_bb[c];
    for(const auto& wall:transport.wall_faces())
        if(wall.reflect_index<transport.n_reflecting_faces())
            x[mesh.n_cells()+wall.reflect_index]=wall.ib;
    const auto snapshot=holo.evaluate(x);
    const auto source_from=[&](const DomHoloSnapshotV47& snap){
        std::vector<Tensor> tensor(mesh.n_cells());
        std::vector<double> J(mesh.n_cells());
        for(std::size_t c=0;c<mesh.n_cells();++c) {
            J[c]=snap.projected[c];
            for(std::size_t k=0;k<6;++k)
                tensor[c][k]=snap.eddington_symmetric[6*c+k];
        }
        DomSmmGrayWallV49 lo(mesh,geo,dirs,bc,std::move(tensor),
                  std::vector<double>(mesh.n_cells(),15.),
                  std::vector<double>(mesh.n_cells(),.1));
        return lo.consistent_source(J,gas_bb,snap.face_current);
    };
    const auto c1=source_from(snapshot);
    require(c1.max_discrete_identity_relative<1e-11,
        "V49: exact HO/LO energy residual identity failed on production FV");
    // Force a new scattering mean. The face correction MUST be updated when
    // the angular source (and therefore its anisotropic E) changes.
    x[0]*=1.04;
    const auto updated=holo.evaluate(x);
    const auto c2=source_from(updated);
    require(c2.max_discrete_identity_relative<1e-11,
        "V49: refreshed HO/LO energy residual identity failed");
    double change=0.;
    for(std::size_t f=0;f<mesh.n_faces();++f)
        change=std::max(change,std::abs(c1.face_ho_minus_lo_watts[f]-
                                        c2.face_ho_minus_lo_watts[f]));
    require(change>1e-10,"V49: HO/LO closure did not refresh with transport");
    std::cout<<"V49_TRUE_HO_GRAY cells="<<mesh.n_cells()
             <<" identity_before="<<c1.max_discrete_identity_relative
             <<" identity_after="<<c2.max_discrete_identity_relative
             <<" refreshed_change="<<change
             <<" physical_residual_before="<<c1.max_ho_energy_residual_relative
             <<" status=RESIDUAL_IDENTITY_ONLY"<<'\n';
}
void affine_oracle() {
    auto mesh=grid(3,3,3,.2);
    const auto geo=build_fv_geometry(mesh);
    const auto directions=make_dom_tensor_gauss_legendre_v39(4,16);
    const auto E=tensors(geo,false);
    const std::size_t nc=mesh.n_cells();
    const std::vector<double> kt(nc,2.),ka(nc,.25);
    auto wall_bc=walls(mesh,360.,false);
    // First build to obtain the known Robin coefficients. Then prescribe wall
    // T from the exact affine flux, not from the numerical LO solution.
    for(auto& bc:wall_bc)bc.second.emissivity=.7;
    const DomSmmGrayWallV49 trial(mesh,geo,directions,wall_bc,E,kt,ka);
    const Vec3 grad={.018,-.013,.011};
    std::vector<double> exact(nc);
    for(std::size_t c=0;c<nc;++c)
        exact[c]=11.+grad.dot(geo.cell_centres[c]);
    const auto& e=E[0];
    const Vec3 Eg={e[0]*grad.x+e[3]*grad.y+e[4]*grad.z,
                   e[3]*grad.x+e[1]*grad.y+e[5]*grad.z,
                   e[4]*grad.x+e[5]*grad.y+e[2]*grad.z};
    // Face-resolved CHT wall emissions are essential on nonmatching grids;
    // the manufactured affine field requires a different Ib for each wall
    // face, not merely one averaged temperature per patch.
    std::vector<double> face_ib(mesh.n_faces(),std::numeric_limits<double>::quiet_NaN());
    for(const auto& w:trial.walls()) {
        const auto& sf=geo.face_area_vectors[w.face];
        const double target=-4.*std::acos(-1.)/kt[w.owner]*sf.dot(Eg);
        require(w.conductance>0.,"V49: affine manufactured wall not anchored");
        const double required_wall_ib=exact[w.owner]-target/w.conductance;
        require(required_wall_ib>0. && std::isfinite(required_wall_ib),
                "V49: affine wall emission unphysical");
        face_ib[w.face]=required_wall_ib;
    }
    const DomSmmGrayWallV49 system(mesh,geo,directions,wall_bc,E,kt,ka,.03,face_ib);
    const auto matrix_response=system.apply(exact);
    const auto b=system.rhs(exact);
    double max_defect=0.,max_b=1.;
    for(std::size_t c=0;c<nc;++c) {
        max_defect=std::max(max_defect,std::abs(matrix_response[c]-b[c]));
        max_b=std::max(max_b,std::abs(b[c]));
    }
    require(max_defect/max_b<1e-10,
            "V49: independent affine FV/gray Robin manufactured oracle failed");
    const auto sol=system.solve(exact,DomSmmV49Preconditioner::AmgAuxiliary,1.e-10,150);
    require(sol.accepted,"V49: affine gray-wall MMS rejected: "+sol.diagnostic);
    double sol_error=0.;
    for(std::size_t c=0;c<nc;++c)
        sol_error=std::max(sol_error,std::abs(sol.accepted_J[c]-exact[c])/
                                        std::max(1.,std::abs(exact[c])));
    require(sol_error<1e-8,"V49: manufactured affine temperature not recovered");
    // Invalid face-resolved emission must be rejected before mutating a state.
    bool refused=false;
    auto bad_face=face_ib;bad_face[trial.walls().front().face]=-1.;
    try { DomSmmGrayWallV49 bad(mesh,geo,directions,wall_bc,E,kt,ka,.03,bad_face); }
    catch(const std::invalid_argument&){refused=true;}
    require(refused,"V49: invalid per-face blackbody was silently accepted");
    std::cout<<"V49_AFFINE_ROBIN cells="<<nc<<" relative_defect="<<max_defect/max_b
             <<" relative_solution="<<sol_error
             <<" AMG_levels="<<sol.preconditioner_levels<<'\n';
}
}
int main(int argc,char** argv) {
    try {
        if(argc>2 || (argc==2 && std::string(argv[1])!="--sanitizer-small"))
            throw std::invalid_argument("V49: unknown arguments");
        const bool small=argc==2;
        test_case(3,3,3,.03,small);
        if(!small)true_dom_gray_source();
        if(!small) {
            test_case(4,4,6,.012,false);
            test_case(8,8,10,.003,false);
            affine_oracle();
        }
        std::cout<<"V49_ALL_PASS\n";
    } catch(const std::exception& e) {
        std::cerr<<"V49_FAIL: "<<e.what()<<'\n';return 1;
    }
}
