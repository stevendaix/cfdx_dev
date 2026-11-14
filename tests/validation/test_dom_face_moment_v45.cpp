// V45: exact face-current zeroth moment + transport-sampled local block
// preconditioner on three 3D mesh families. Diagnostic, NOT production DSA.
#include "cfdx/physics/dom_face_moment_v45.h"
#include "cfdx/physics/dom_quadrature_v39.h"
#include "thermal_3d_mesh_v36.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
void check(bool v,const char* why) {if(!v)throw std::runtime_error(why);}
using Cell=DomMomentSchurSystemV44::CellField;
std::string spectrum_prefix;
Mesh make_mesh(double length_z) {
    auto mesh=cfdx::testing::thermal_v36::make_mesh(3,3,3,.018,0.,length_z);
    const double ca=std::cos(.41),sa=std::sin(.41);
    const double cb=std::cos(.33),sb=std::sin(.33);
    for(std::size_t p=0;p<mesh.points().size();++p) {
        const double x=mesh.points().x_data()[p],y=mesh.points().y_data()[p];
        const double z=mesh.points().z_data()[p];
        const double xr=ca*x-sa*y,yr=sa*x+ca*y;
        mesh.points().set(p,cb*xr+sb*z,yr,-sb*xr+cb*z);
    }
    return mesh;
}
DomWallBoundaryConditions walls_for(const Mesh& mesh,bool gray) {
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& name=mesh.boundary().patch(p).name;
        const double T=name=="west"?540.:name=="east"?285.:365.;
        const double eps=!gray?1.:name=="west"?.35:name=="east"?.65:
            name=="south"?.90:name=="north"?1.:name=="bottom"?0.:1.;
        walls[name]={eps,T};
    }
    return walls;
}
bool upwind_graph_is_dag(const Mesh& mesh,const FvGeometry& geo,
                         const std::vector<DiscreteDirection>& directions) {
    const auto nc=mesh.n_cells();
    for(const auto& d:directions) {
        std::vector<std::vector<std::size_t>> edges(nc);
        std::vector<std::size_t> incoming(nc,0);
        for(std::size_t f=0;f<mesh.n_faces();++f) {
            const auto ni=mesh.ownership().neighbour(f);
            if(ni<0)continue;
            const auto& sf=geo.face_area_vectors[f];
            const double mu=d.dx*sf.x+d.dy*sf.y+d.dz*sf.z;
            if(mu==0.)continue;
            const auto o=mesh.ownership().owner(f);
            const auto n=static_cast<std::size_t>(ni);
            const auto from=mu>0.?o:n,to=mu>0.?n:o;
            edges[from].push_back(to);
            ++incoming[to];
        }
        std::queue<std::size_t> pending;
        for(std::size_t c=0;c<nc;++c) if(incoming[c]==0)pending.push(c);
        std::size_t traversed=0;
        while(!pending.empty()) {
            const auto c=pending.front();pending.pop();++traversed;
            for(const auto next:edges[c]) if(--incoming[next]==0)pending.push(next);
        }
        if(traversed!=nc)return false;
    }
    return true;
}
void test_case(double length_z,bool gray,double scattering=14.9) {
    auto mesh=make_mesh(length_z);
    const auto geo=build_fv_geometry(mesh);
    const auto dirs=make_dom_tensor_gauss_legendre_v39(4,12);
    const auto walls=walls_for(mesh,gray);
    check(upwind_graph_is_dag(mesh,geo,dirs),
        "V45: fixture has a cyclic upwind graph; local Schur equality must be requalified");
    Cell temperature(mesh.n_cells(),"gas","K",1);
    for(std::size_t c=0;c<mesh.n_cells();++c)
        temperature(c)=370.+30.*geo.cell_centres[c].x+10.*geo.cell_centres[c].y;
    RadiationTransportControls ctl;
    ctl.absorption=.1;ctl.scattering=scattering;
    ctl.tolerance=scattering>100.?1e-10:1e-9;ctl.max_iterations=500;
    ctl.linear_tolerance=1e-12;ctl.linear_max_iterations=2000;
    DomMomentSchurSystemV44 sys(mesh,geo,temperature,dirs,walls,ctl);
    const auto plain=sys.solve(false,12,140);
    check(plain.converged,"V45: plain exact Schur failed");
    std::vector<double> moment=plain.scalar_mean_intensity;
    moment.insert(moment.end(),plain.reflected_face_mean_intensity.begin(),
        plain.reflected_face_mean_intensity.end());
    std::vector<Cell> intensity;
    (void)sys.project(moment,true,&intensity);
    const auto report=audit_dom_face_moments_v45(
        mesh,geo,temperature,dirs,walls,ctl,sys,moment,intensity);
    std::cout<<std::setprecision(13)<<"V45_FACE z="<<length_z
        <<" gray="<<gray<<" scattering="<<scattering<<" fv="<<report.max_transport_balance_relative
        <<" physical="<<report.max_physical_balance_relative
        <<" low_order="<<report.max_low_order_identity_relative
        <<" energy="<<report.global_energy_relative
        <<" mean="<<report.max_moment_mismatch_relative<<'\n';
    const bool energy_qualified=dom_face_moment_qualified_v45(report,ctl.tolerance);
    if(scattering<=100.) check(energy_qualified,
        "V45: independently evaluated face moment gates failed");
    else {
        std::cout<<"V45_STRESS_"<<(energy_qualified?"QUALIFIED":"UNQUALIFIED")
            <<" global_relative="<<report.global_energy_relative
            <<" gate="<<ctl.tolerance<<'\n';
    }
    check(report.max_transport_balance_relative<ctl.tolerance,"V45: discrete transport current not balanced");
    check(report.max_physical_balance_relative<ctl.tolerance,"V45: physical energy not balanced");
    if(scattering<=100.) check(report.global_energy_relative<ctl.tolerance,"V45: global energy not closed");
    check(report.max_low_order_identity_relative<2e-13,"V45: low-order face closure not identical");
    check(report.max_moment_mismatch_relative<ctl.tolerance,"V45: high order projection differs from Schur");
    double drift=0.;
    for(const double v:report.closure_face_current)
        drift=std::max(drift,std::abs(v));
    check(drift>1e-9,"V45: high/low closure trivially zero");
    // Negative: a moment change with unchanged angular fields must be rejected.
    auto wrong=moment; wrong[0]+=300.;
    std::vector<Cell> current_stale;
    (void)sys.project(wrong,true,&current_stale);
    const auto staggered=audit_dom_face_moments_v45(
        mesh,geo,temperature,dirs,walls,ctl,sys,wrong,current_stale);
    check(staggered.max_transport_balance_relative<1e-9,
        "V45: frozen scattering transport face balance not satisfied");
    check(staggered.max_physical_balance_relative>1e-5,
        "V45: lagged scattering residual incorrectly classified as physical equilibrium");
    check(staggered.max_low_order_identity_relative<2e-13,
        "V45: LO face closure broke during lagged transport");
    bool mismatched=false;
    try {
        auto altered=ctl; altered.scattering+=.01;
        (void)audit_dom_face_moments_v45(
            mesh,geo,temperature,dirs,walls,altered,sys,moment,intensity);
    } catch(const std::invalid_argument&) {mismatched=true;}
    check(mismatched,"V45: mismatched transport parameters accepted");
    const auto bad=audit_dom_face_moments_v45(
        mesh,geo,temperature,dirs,walls,ctl,sys,wrong,intensity);
    check(bad.max_transport_balance_relative>1e-5,
        "V45: independently re-evaluated physical transport did not detect stale moment");
    bool rejected=false;
    try {
        wrong=moment;wrong[0]=std::numeric_limits<double>::quiet_NaN();
        (void)audit_dom_face_moments_v45(
            mesh,geo,temperature,dirs,walls,ctl,sys,wrong,intensity);
    } catch(const std::invalid_argument&) {rejected=true;}
    check(rejected,"V45: NaN moment accepted by current audit");
    bool bad_callback=false;
    try {
        (void)sys.solve(true,2,1,1.,
            [](const std::vector<double>& x){return x;});
    } catch(const std::invalid_argument&) {bad_callback=true;}
    check(bad_callback,"V45: conflicting Schur preconditioner modes accepted");
    bad_callback=false;
    try {
        (void)sys.solve(false,2,1,1.,
            [](const std::vector<double>& x){auto v=x;v.pop_back();return v;});
    } catch(const std::invalid_argument&) {bad_callback=true;}
    check(bad_callback,"V45: invalid preconditioner dimension accepted");
    DomLocalAngularJacobiV45 cheap(sys,geo,dirs,ctl);
    const auto cheap_result=sys.solve(false,12,140,1.,
        [&cheap](const std::vector<double>& r){return cheap.apply(r);});
    std::cout<<"V45_CHEAP_STATUS scattering="<<scattering
        <<" converged="<<cheap_result.converged
        <<" schur="<<cheap_result.true_schur_relative
        <<" local="<<cheap_result.local_transport_relative
        <<" energy="<<cheap_result.global_energy_relative<<'\n';
    if(scattering<=100.)
        check(cheap_result.converged,"V45: local angular Jacobi GMRES did not converge");
    else if(!cheap_result.converged) {
        check(cheap_result.scalar_mean_intensity.empty(),
            "V45: underqualified preconditioner exposed a solution");
        const auto fallback=sys.solve(false,12,140);
        check(fallback.converged,"V45: guarded fallback did not restore full Schur");
    }
    if(cheap_result.converged) check(
          cheap_result.local_transport_relative<ctl.tolerance &&
          cheap_result.global_energy_relative<ctl.tolerance,
          "V45: cheap angular block failed FV physical acceptance");
    double cheap_diff=0.,cheap_norm=1.;
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        if(cheap_result.converged) cheap_diff=std::max(cheap_diff,std::abs(
            cheap_result.scalar_mean_intensity[c]-plain.scalar_mean_intensity[c]));
        cheap_norm=std::max(cheap_norm,std::abs(plain.scalar_mean_intensity[c]));
    }
    if(cheap_result.converged) check(cheap_diff/cheap_norm<1e-6,"V45: cheap angular Jacobi changed physical state");
    std::cout<<"V45_CHEAP z="<<length_z<<" gray="<<gray
        <<" scattering="<<scattering<<" plain_krylov="<<plain.krylov_iterations
        <<" cheap_krylov="<<cheap_result.krylov_iterations
        <<" plain_maps="<<plain.transport_maps
        <<" cheap_maps="<<cheap_result.transport_maps
        <<" setup_s="<<cheap.setup_seconds()
        <<" solve_s="<<cheap_result.seconds
        <<" total_s="<<cheap.setup_seconds()+cheap_result.seconds
        <<" compared="<<cheap_diff/cheap_norm<<'\n';
    DomLocalSchurBlocksV45 blocks(sys);
    rejected=false;
    try {(void)DomLocalSchurBlocksV45(sys,20);}
    catch(const std::invalid_argument&) {rejected=true;}
    check(rejected,"V45: dense diagnostic size cap failed");
    check(blocks.dimension()==sys.dimension() &&
        blocks.probe_maps()==sys.dimension(),"V45: inaccurate exact Schur probing count");
    // Each cell+reflecting-face block is inverted by a stable pivoted LU.
    std::vector<double> rhs(blocks.dimension());
    for(std::size_t k=0;k<rhs.size();++k) rhs[k]=std::sin(double(k)+.12);
    const auto z=blocks.apply(rhs);
    const auto& H=blocks.dense_schur();
    double error=0.;
    for(const auto& b:blocks.blocks())
        for(const auto i:b.ids) {
            double value=0.;
            for(const auto j:b.ids) value+=H[i*blocks.dimension()+j]*z[j];
            error=std::max(error,std::abs(value-rhs[i]));
        }
    check(error<2e-11,"V45: local block LU does not invert sampled Schur block");
    const auto zcheap=cheap.apply(rhs);
    double block_delta=0.,block_scale=1.;
    for(std::size_t i=0;i<z.size();++i) {
        block_delta=std::max(block_delta,std::abs(z[i]-zcheap[i]));
        block_scale=std::max(block_scale,std::abs(z[i]));
    }
    // On these three acyclic upwind meshes the inverse angular diagonal
    // equals 1/A_m,cc. A cyclic transport graph would invalidate this claim.
    check(block_delta/block_scale<3e-11,
        "V45: local angular block differs from probed exact cell block");

    if(!spectrum_prefix.empty()) {
        const auto prefix=spectrum_prefix+(gray?"_gray":"_black")+
            "_z"+std::to_string(length_z)+"_ks"+std::to_string(scattering);
        auto save_matrix=[&](const std::string& suffix,
                             const std::vector<double>& vals) {
            std::ofstream file(prefix+suffix+".csv");
            if(!file) throw std::runtime_error("V45: cannot write spectrum fixture");
            file<<std::setprecision(17);
            for(std::size_t row=0;row<blocks.dimension();++row) {
                for(std::size_t col=0;col<blocks.dimension();++col) {
                    if(col)file<<',';
                    file<<vals[row*blocks.dimension()+col];
                }
                file<<'\n';
            }
        };
        save_matrix("_H",blocks.dense_schur());
        const std::size_t n=blocks.dimension();
        std::vector<double> exact(n*n,0.),cheap_matrix(n*n,0.);
        for(std::size_t col=0;col<n;++col) {
            std::vector<double> hc(n);
            for(std::size_t row=0;row<n;++row) hc[row]=H[row*n+col];
            const auto pe=blocks.apply(hc),pc=cheap.apply(hc);
            for(std::size_t row=0;row<n;++row) {
                exact[row*n+col]=pe[row];cheap_matrix[row*n+col]=pc[row];
            }
        }
        save_matrix("_Psampled_invH",exact);
        save_matrix("_Pcheap_invH",cheap_matrix);
    }
    const double row_bound=blocks.max_row_sum_I_minus_preconditioned_H();
    const auto accelerated=sys.solve(false,12,140,1.,
        [&blocks](const std::vector<double>& r){return blocks.apply(r);});
    std::cout<<"V45_PROBED_STATUS scattering="<<scattering
        <<" converged="<<accelerated.converged
        <<" schur="<<accelerated.true_schur_relative
        <<" local="<<accelerated.local_transport_relative
        <<" energy="<<accelerated.global_energy_relative<<'\n';
    if(scattering<=100.) check(accelerated.converged,
        "V45: sampled blocks did not converge");
    else if(!accelerated.converged) check(accelerated.scalar_mean_intensity.empty(),
        "V45: sampled Schur published unaccepted state");
    if(accelerated.converged) check(
          accelerated.local_transport_relative<ctl.tolerance &&
          accelerated.global_energy_relative<ctl.tolerance,
          "V45: preconditioned state fails independent physical acceptance");
    double mismatch=0.,scale=1.;
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        if(accelerated.converged) mismatch=std::max(mismatch,std::abs(
            accelerated.scalar_mean_intensity[c]-plain.scalar_mean_intensity[c]));
        scale=std::max(scale,std::abs(plain.scalar_mean_intensity[c]));
    }
    if(accelerated.converged)
        check(mismatch/scale<1e-6,"V45: preconditioner changed accepted DOM state");
    std::cout<<"V45_BLOCK z="<<length_z<<" gray="<<gray
        <<" scattering="<<scattering<<" dim="<<blocks.dimension()<<" plain_krylov="<<plain.krylov_iterations
        <<" local_krylov="<<accelerated.krylov_iterations
        <<" probes="<<blocks.probe_maps()
        <<" plain_maps="<<plain.transport_maps
        <<" local_maps="<<accelerated.transport_maps
        <<" row_bound="<<row_bound
        <<" state_difference="<<mismatch/scale
        <<" baseline_s="<<plain.seconds
        <<" setup_s="<<blocks.setup_seconds()
        <<" solve_s="<<accelerated.seconds
        <<" total_s="<<blocks.setup_seconds()+accelerated.seconds<<'\n';
}
}
int main(int argc,char** argv) {
    try {
        if(argc==3 && std::string(argv[1])=="--spectrum") spectrum_prefix=argv[2];
        else if(argc!=1) throw std::invalid_argument("V45: --spectrum PREFIX only");
        test_case(1.,false);       // near-isotropic 3D, black walls
        test_case(.03,true);       // strongly stretched, gray/reflecting walls
        test_case(.006,true);      // extreme stretching, same wall physics
        test_case(.03,true,199.9); // optically thicker, scattering dominant
        std::cout<<"V45_DOM_FACE_MOMENTS PASS\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"V45_DOM_FACE_MOMENTS FAIL: "<<e.what()<<'\n';
        return 1;
    }
}
