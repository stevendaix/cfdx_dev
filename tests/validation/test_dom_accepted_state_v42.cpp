// V42: physically re-evaluated DOM scattering/gray-wall balance, on a skew 3D
// mesh. The oracle uses conservative upwind face fluxes, not a Picard delta.
#include "cfdx/physics/radiation_solver.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
using CellField=Field<double,Location::CELL>;
void require(bool condition,const char* why){
    if(!condition) throw std::runtime_error(why);
}
std::vector<DiscreteDirection> eight_ordinates(){
    std::vector<DiscreteDirection> dirs;
    const double a=1.0/std::sqrt(3.0);
    for(int sx:{-1,1}) for(int sy:{-1,1}) for(int sz:{-1,1})
        dirs.push_back({sx*a,sy*a,sz*a,M_PI/2.0});
    return dirs;
}
DomWallBoundaryConditions wall_data(const Mesh& mesh,bool isothermal){
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
        const auto& name=mesh.boundary().patch(p).name;
        walls[name]={0.58,isothermal?370.0:
            name=="west" ? 480.0 : name=="east" ? 285.0 : 340.0};
    }
    return walls;
}
void isothermal_oracle(){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(3,3,3,0.023,0.,0.035);
    auto geo=build_fv_geometry(mesh);
    const auto dirs=eight_ordinates();
    const auto walls=wall_data(mesh,true);
    CellField T(mesh.n_cells(),"T","K",1);
    T.fill(370.0);
    RadiationTransportControls c;
    c.absorption=0.25;c.scattering=7.0;
    std::vector<CellField> I;
    for(std::size_t m=0;m<dirs.size();++m){
        I.emplace_back(mesh.n_cells(),"I","W/m2/sr",1);
        I.back().fill(blackbody_intensity(370.0));
    }
    const auto exact=audit_dom_accepted_state_v42(mesh,geo,T,I,dirs,c,walls);
    require(exact.max_cell_relative<1e-12,"isothermal 3D scattering FV defect");
    require(exact.relative_global_energy<1e-11,"isothermal closed-wall energy");
    I[0](0)*=1.10;
    std::vector<double> face_defects;
    const auto corrupt=audit_dom_accepted_state_v42(
        mesh,geo,T,I,dirs,c,walls,&face_defects);
    require(face_defects.size()==dirs.size()*mesh.n_cells(),
        "DOM V42 face-defect map is incomplete");
    // Compare independent oriented face reconstruction with the ACTUAL
    // production FV matrix reassembled on the same final I and scattering J.
    Field<double,Location::FACE> F(mesh.n_faces(),"angle","m2",1);
    const auto& d=dirs[0];
    for(std::size_t f=0;f<mesh.n_faces();++f) {
        const auto Sf=geo.face_area_vectors[f];
        F(f)=d.dx*Sf.x+d.dy*Sf.y+d.dz*Sf.z;
    }
    CellField Su(mesh.n_cells(),"Su","",1),Sp(mesh.n_cells(),"Sp","",1);
    for(std::size_t cell=0;cell<mesh.n_cells();++cell) {
        double J=0.;
        for(std::size_t m=0;m<dirs.size();++m)
            J+=dirs[m].weight*I[m](cell)/(4.*M_PI);
        Su(cell)=c.absorption*blackbody_intensity(T(cell))+c.scattering*J;
        Sp(cell)=-(c.absorption+c.scattering);
    }
    const auto face_values=build_dom_diffuse_gray_face_values(
        mesh,geo,dirs,I,walls,0);
    ScalarBoundaryConditions bc;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p)
        bc[mesh.boundary().patch(p).name]={ScalarBoundaryType::ZERO_GRADIENT,0.,0.};
    const auto eq=assemble_scalar_equation(mesh,geo,F,0.,Su,Sp,bc,
        true,&face_values,nullptr,nullptr,nullptr,ConvectionScheme::UPWIND,&I[0]);
    const std::vector<double> x(I[0].component_data(0),
        I[0].component_data(0)+mesh.n_cells());
    const auto Ax=eq.matrix.matvec(x);
    double assembly_error=0.;
    for(std::size_t cell=0;cell<mesh.n_cells();++cell)
        assembly_error=std::max(assembly_error,
            std::abs(Ax[cell]-eq.rhs(cell)-face_defects[cell]));
    require(assembly_error<5e-11,
        "independent FV face residual disagrees with production matrix");
    require(corrupt.max_cell_relative>1e-4,
        "local scattered intensity corruption escaped independent audit");
    require(corrupt.max_cell_relative>exact.max_cell_relative*1e5,
        "corruption did not increase physical residual");
    std::cout<<"V42_ISOTHERMAL exact="<<exact.max_cell_relative
             <<" corrupt="<<corrupt.max_cell_relative
             <<" matrix_face_difference="<<assembly_error<<'\n';
}
void nonisothermal_scattering_rollback(){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(3,3,3,0.018,0.,0.04);
    const auto geo=build_fv_geometry(mesh);
    const auto dirs=eight_ordinates();
    const auto walls=wall_data(mesh,false);
    const auto quadrature=audit_dom_wall_quadrature(mesh,geo,dirs);
    require(quadrature.maximum_relative_pi_error>0.1,
        "coarse DOM quadrature unexpectedly passed physical wall accuracy");
    CellField T(mesh.n_cells(),"T","K",1);
    CellField G(mesh.n_cells(),"G","W/m2",1);
    CellField q(mesh.n_cells(),"q","W/m3",1);
    for(std::size_t c=0;c<mesh.n_cells();++c) T(c)=355.+0.6*double(c);
    G.fill(128.75);q.fill(-39.25);
    const auto original_G=G,original_q=q;
    RadiationTransportControls ctl;
    ctl.absorption=0.20;ctl.scattering=5.0;
    ctl.max_iterations=1;ctl.tolerance=1e-9;
    ctl.linear_tolerance=1e-12;
    const auto rejected=solve_participating_radiation_diffuse_gray_walls(
        mesh,geo,T,G,q,dirs,ctl,walls,true,true,true);
    require(!rejected.converged,"one Picard sweep incorrectly accepted");
    for(std::size_t c=0;c<mesh.n_cells();++c){
        require(G(c)==original_G(c) && q(c)==original_q(c),
            "strict DOM rejected trial mutated caller fields");
    }
    ctl.max_iterations=500;
    const auto accepted=solve_participating_radiation_diffuse_gray_walls(
        mesh,geo,T,G,q,dirs,ctl,walls,true,true,true);
    require(accepted.converged && accepted.accepted_state_verified_v42,
        "high-albedo DOM did not achieve physical FV convergence");
    require(accepted.max_cell_transport_residual_v42<ctl.tolerance,
        "local accepted-state FV residual above tolerance");
    require(accepted.relative_energy_closure_v42<ctl.tolerance,
        "accepted-state net radiative energy not conserved");
    const auto accepted_G=G,accepted_q=q;
    G=original_G;q=original_q;
    const auto repeated=solve_participating_radiation_diffuse_gray_walls(
        mesh,geo,T,G,q,dirs,ctl,walls,true,true,true);
    require(repeated.converged && repeated.iterations==accepted.iterations,
        "V42 strict CPU serial solution is not repeatable");
    for(std::size_t c=0;c<mesh.n_cells();++c)
        require(G(c)==accepted_G(c) && q(c)==accepted_q(c),
            "V42 strict CPU serial field is not reproducible");
    G=original_G;q=original_q;
    const auto legacy=solve_participating_radiation_diffuse_gray_walls(
        mesh,geo,T,G,q,dirs,ctl,walls,true,true,false);
    require(legacy.converged,"legacy DOM A/B did not converge");
    double max_G=0.0,max_q=0.0;
    for(std::size_t c=0;c<mesh.n_cells();++c){
        max_G=std::max(max_G,std::abs(G(c)-accepted_G(c)));
        max_q=std::max(max_q,std::abs(q(c)-accepted_q(c)));
        require(std::abs(G(c)-accepted_G(c))<2e-5,
            "strict/legacy final radiation mismatch");
        require(std::abs(q(c)-accepted_q(c))<1e-5,
            "strict/legacy final source mismatch");
    }
    // Independent global energy accounting includes the sum over actual wall
    // patches; wall power in the final accepted result and the cell sources
    // must close with identical signs.
    long double wall=0.0L,gas=0.0L;
    for (const auto& [name,power]:accepted.wall_outward_W){
        (void)name;wall+=power;
    }
    for (std::size_t c=0;c<mesh.n_cells();++c)
        gas+=static_cast<long double>(accepted_q(c))*geo.cell_volumes[c];
    require(std::abs(wall-gas)<1e-8L*std::max(1.0L,std::abs(wall)),
        "accepted DOM wall/gas power mismatch");
    std::cout<<std::setprecision(12)
             <<"V42_SCATTER albedo="<<(ctl.scattering/(ctl.scattering+ctl.absorption))
             <<" legacy_iters="<<legacy.iterations
             <<" strict_iters="<<accepted.iterations
             <<" max_fv="<<accepted.max_cell_transport_residual_v42
             <<" energy_relative="<<accepted.relative_energy_closure_v42
             <<" wall_W="<<double(wall)<<" gas_W="<<double(gas)
             <<" max_G_diff="<<max_G<<" max_q_diff="<<max_q
             <<" max_wall_pi_error="<<quadrature.maximum_relative_pi_error<<'\n';
    // Rejected invalid data and linear failure must not leak mutations.
    G=original_G;q=original_q;
    auto invalid=ctl;
    invalid.linear_max_iterations=1; invalid.linear_tolerance=1e-30;
    bool threw=false;
    try { (void)solve_participating_radiation_diffuse_gray_walls(
        mesh,geo,T,G,q,dirs,invalid,walls,true,true,true); }
    catch(const std::exception&){threw=true;}
    require(threw,"expected linear solve failure was silently accepted");
    for(std::size_t c=0;c<mesh.n_cells();++c)
        require(G(c)==original_G(c) && q(c)==original_q(c),
            "failed strict DOM solve leaked caller state");
}
}
int main(){
    try{
        isothermal_oracle();
        nonisothermal_scattering_rollback();
        std::cout<<"V42_DOM_ACCEPTED_STATE PASS\n";
        return 0;
    }catch(const std::exception& e){
        std::cerr<<"V42_DOM_ACCEPTED_STATE FAIL: "<<e.what()<<'\n';
        return 1;
    }
}
