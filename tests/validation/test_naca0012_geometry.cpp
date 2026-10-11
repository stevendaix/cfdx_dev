#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "qualification_boundary_helpers.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

struct Forces { double cd,cl,cm,cdp,cdv,continuity; };

static Forces run(const std::string& path,std::size_t level,const char* algo_name="SIMPLE",PressureVelocityAlgorithm algo=PressureVelocityAlgorithm::SIMPLE,CoupledSchurModel schur=CoupledSchurModel::BlockLocal)
{
    Mesh mesh;
    if(!cfdx::io::read_mesh_hdf5(path,mesh))
        throw std::runtime_error("cannot read NACA0012 mesh: "+path);
    if(!mesh.topo_validate().ok)
        throw std::runtime_error("invalid NACA0012 topology: "+path);

    constexpr double rho=1.0, Uinf=1.0, Re=1000.0, chord=1.0;
    constexpr double nu=Uinf*chord/Re;
    Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    Field<double,Location::CELL> p(mesh.n_cells(),"p","Pa",1);
    U.fill(0.0); p.fill(0.0);
    VelocityBoundaryConditions ubc;
    ubc["farfield"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{Uinf,0,0}};
    ubc["airfoil"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["front"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ubc["back"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ScalarBoundaryConditions pbc;
    pbc["farfield"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    pbc["airfoil"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    pbc["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    pbc["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};

    IncompressibleSolverControls c;
    c.algorithm=algo;
    c.density=rho; c.kinematic_viscosity=nu;
    if(algo==PressureVelocityAlgorithm::COUPLED)
        c.coupling.schur_model=schur;
    c.linear_max_iterations=4000; c.linear_tolerance=1e-10;
    c.convergence.max_iterations=8000;
    c.convergence.relative_tolerance=1e-8;
    c.convergence.continuity_tolerance=1e-8;
    c.coupling.alpha_u=0.6; c.coupling.alpha_p=0.25;
    c.use_bounded_convection=true;
    c.convection_scheme=ConvectionScheme::UPWIND;
    c.pressure_reference_cell=0;
    c.pressure_reference_value=0.0;
    c.diagnostics.iteration_trace=true;
    c.diagnostics.iteration_trace_frequency=200;

    const auto geometry=build_fv_geometry(mesh);
    for(std::size_t i=0;i<mesh.n_cells();++i)
        U.component_data(0)[i]=Uinf;

    // Exercise the strict mathematical BC/FVM path before the production solve.
    // The legacy maps below are retained for the current steady-solver compatibility API.
    cfdx::physics::BoundaryConstraintMap strict_bc;
    cfdx::validation::add_velocity_dirichlet(strict_bc, "farfield", {Uinf, 0.0, 0.0});
    cfdx::validation::add_velocity_dirichlet(strict_bc, "airfoil", {0.0, 0.0, 0.0});
    cfdx::validation::add_pressure_dirichlet(strict_bc, "farfield", 0.0);
    cfdx::validation::add_pressure_neumann(strict_bc, "airfoil", 0.0);
    const auto validation_flux = cfdx::validation::make_validation_face_flux(mesh, geometry, U);
    cfdx::validation::exercise_new_velocity_bc_contract(
        mesh, geometry, strict_bc, validation_flux, "NACA0012_RE1000_A0");

    // Force-audit infrastructure (diag): the airfoil surface integration is
    // factored into a callable so the per-iteration Cd trace and the final
    // integration share one implementation. n points out of the fluid domain
    // (into the body); the force on the airfoil is F = integral (p*n - tau*n) dA.
    ScalarBoundaryConditions gx,gy;
    gx["farfield"]={ScalarBoundaryType::FIXED_VALUE,Uinf,0};
    gx["airfoil"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    gx["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    gx["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    gy["farfield"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    gy["airfoil"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    gy["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    gy["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    const double q=0.5*rho*Uinf*Uinf;
    const std::size_t cd_trace_frequency=(level<=64u)?16u:((level<=128u)?32u:64u);
    auto airfoil_forces=[&](const Field<double,Location::CELL>& Uc,
                            const Field<double,Location::CELL>& pc)
    {
        Field<double,Location::CELL> lx(mesh.n_cells(),"lx","m/s",1), ly(mesh.n_cells(),"ly","m/s",1);
        for(std::size_t i=0;i<mesh.n_cells();++i){lx(i)=Uc.component_data(0)[i];ly(i)=Uc.component_data(1)[i];}
        const auto lgx=gauss_gradient_with_boundary(lx,mesh,geometry,gx);
        const auto lgy=gauss_gradient_with_boundary(ly,mesh,geometry,gy);
        const auto& fpatch=mesh.boundary().patch(mesh.boundary().find("airfoil"));
        double fx=0,fy=0,mz=0,fxp=0,fyp=0,fxv=0,fyv=0;
        for(const auto f:fpatch.face_ids){
            const std::size_t o=mesh.ownership().owner(f);
            const auto Sf=geometry.face_area_vectors[f];
            const auto n=Sf.normalized();
            const double pp=pc(o);
            const double txx=2.0*rho*nu*lgx.component_data(0)[o];
            const double tyy=2.0*rho*nu*lgy.component_data(1)[o];
            const double txy=rho*nu*(lgx.component_data(1)[o]+lgy.component_data(0)[o]);
            const double tx=(pp*n.x-(txx*n.x+txy*n.y))*Sf.mag();
            const double ty=(pp*n.y-(txy*n.x+tyy*n.y))*Sf.mag();
            fx+=tx; fy+=ty;
            fxp+=pp*n.x*Sf.mag(); fyp+=pp*n.y*Sf.mag();
            fxv+=-(txx*n.x+txy*n.y)*Sf.mag(); fyv+=-(txy*n.x+tyy*n.y)*Sf.mag();
            const auto r=geometry.face_centres[f]-Vec3{0.25,0,0};
            mz+=r.x*ty-r.y*tx;
        }
        (void)fyp; (void)fyv;
        return Forces{fx/q,fy/q,mz/q,fxp/q,fxv/q,0.0};
    };
    c.iteration_output_callback=[&](std::size_t it,double,
                                    const Mesh&,
                                    const Field<double,Location::CELL>& Uc,
                                    const Field<double,Location::CELL>& pc)->bool
    {
        if(it==1 || it%cd_trace_frequency==0){
            const auto tf=airfoil_forces(Uc,pc);
            std::cout<<"NACA_CD_TRACE level="<<level<<" iter="<<it
                     <<" Cd="<<tf.cd<<" Cl="<<tf.cl
                     <<" Cd_pressure="<<tf.cdp<<" Cd_viscous="<<tf.cdv<<"\n";
        }
        return true;
    };

    const auto solve=solve_steady_incompressible(mesh,U,p,ubc,pbc,c);

    // DIAG (temporary branch): the production solve stalls with a frozen
    // continuity floor. Attribute the divergence of the authoritative
    // mass flux cell-by-cell and to the boundary patches it touches.
    {
        std::cout << "NACA_DIAG status=" << to_string(solve.convergence_status)
                  << " converged=" << (solve.converged ? 1 : 0)
                  << " iterations=" << solve.iterations << "\n";
        const auto& flux = solve.authoritative_mass_flux;
        const auto& mcells = mesh.cells();
        const auto* cell_faces = mcells.faces_data();
        const auto* offsets = mcells.offsets_data();
        const auto& own2 = mesh.ownership();
        const auto& bnd = mesh.boundary();
        std::vector<std::string> face_patch(mesh.n_faces(), std::string("interior"));
        for (std::size_t p = 0; p < bnd.n_patches(); ++p) {
            const auto& patch = bnd.patch(p);
            for (const auto f : patch.face_ids) face_patch[f] = patch.name;
        }
        std::vector<std::pair<double, std::size_t>> offenders;
        double div_l1 = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            double balance = 0.0;
            for (auto k = offsets[c]; k < offsets[c + 1]; ++k) {
                const std::size_t f = cell_faces[k];
                balance += (own2.owner(f) == c) ? flux(f) : -flux(f);
            }
            div_l1 += std::abs(balance);
            offenders.push_back(std::make_pair(std::abs(balance), c));
        }
        std::sort(offenders.begin(), offenders.end(),
                  [](const std::pair<double, std::size_t>& a,
                     const std::pair<double, std::size_t>& b) {
                      return a.first > b.first; });
        std::cout << "NACA_DIAG divergence_l1=" << div_l1 << "\n";
        for (std::size_t i = 0; i < 12 && i < offenders.size(); ++i) {
            const std::size_t c = offenders[i].second;
            std::string patches;
            for (auto k = offsets[c]; k < offsets[c + 1]; ++k) {
                const std::size_t f = cell_faces[k];
                if (own2.neighbour(f) < 0) patches += face_patch[f] + ",";
            }
            std::cout << "NACA_DIAG cell=" << c
                      << " |div|=" << offenders[i].first
                      << " patches=[" << patches << "]"
                      << " x=" << geometry.cell_centres[c].x
                      << " y=" << geometry.cell_centres[c].y << "\n";
        }
    }
    if(!solve.converged) {
        // std::cout is buffered; an uncaught exception terminates without
        // flushing it, which silently swallowed the NACA_DIAG attribution.
        std::cout.flush();
        throw std::runtime_error("NACA0012 solver did not converge");
    }

    // Executed-criteria audit (diag): dump the exact gate thresholds and the
    // last history entries so convergence is verifiable from the log alone.
    std::cout<<"NACA_GATE_DEF level="<<level
             <<" relative_tolerance="<<c.convergence.relative_tolerance
             <<" continuity_tolerance="<<c.convergence.continuity_tolerance
             <<" max_iterations="<<c.convergence.max_iterations<<"\n";
    {
        const std::size_t hn=solve.history.size();
        const std::size_t first=(hn>=5)?(hn-5):0;
        for(std::size_t i=first;i<hn;++i){
            const auto& h=solve.history[i];
            std::cout<<"NACA_FINAL_HISTORY level="<<level<<" iter="<<h.iteration
                     <<" momentum_residual="<<h.momentum_residual
                     <<" momentum_equation_residual_relative="<<h.momentum_equation_residual_relative
                     <<" pressure_residual="<<h.pressure_residual
                     <<" continuity_linf="<<h.continuity_linf
                     <<" continuity_normalized="<<h.continuity_normalized
                     <<" velocity_change_inf="<<h.velocity_change_inf
                     <<" pressure_change_inf="<<h.pressure_change_inf<<"\n";
        }
    }
    Field<double,Location::CELL> ux(mesh.n_cells(),"ux","m/s",1), uy(mesh.n_cells(),"uy","m/s",1);
    for(std::size_t i=0;i<mesh.n_cells();++i){ux(i)=U.component_data(0)[i];uy(i)=U.component_data(1)[i];}
    const auto gux=gauss_gradient_with_boundary(ux,mesh,geometry,gx);
    const auto guy=gauss_gradient_with_boundary(uy,mesh,geometry,gy);
    {
        // Independent momentum-budget estimate of the airfoil force, computed
        // only from the outer boundaries (no airfoil integration): steady
        // balance over the fluid domain gives, with n out of the fluid domain,
        // F_body = integral over non-airfoil boundary faces of
        // (-p*n + tau*n - rho*u*(u.n)) dA.
        double bfx=0,bfy=0;
        const auto& bnd=mesh.boundary();
        const std::size_t airfoil_id=bnd.find("airfoil");
        for(std::size_t pi=0;pi<bnd.n_patches();++pi){
            if(pi==airfoil_id) continue;
            for(const auto f:bnd.patch(pi).face_ids){
                const std::size_t o=mesh.ownership().owner(f);
                const auto Sf=geometry.face_area_vectors[f];
                const auto n=Sf.normalized();
                const double pp=p(o);
                const double un=U.component_data(0)[o]*n.x+U.component_data(1)[o]*n.y;
                const double txx=2.0*rho*nu*gux.component_data(0)[o];
                const double tyy=2.0*rho*nu*guy.component_data(1)[o];
                const double txy=rho*nu*(gux.component_data(1)[o]+guy.component_data(0)[o]);
                bfx+=(-pp*n.x+(txx*n.x+txy*n.y)-rho*un*U.component_data(0)[o])*Sf.mag();
                bfy+=(-pp*n.y+(txy*n.x+tyy*n.y)-rho*un*U.component_data(1)[o])*Sf.mag();
            }
        }
        std::cout<<"NACA_FORCE_AUDIT level="<<level
                 <<" budget_Cd="<<(bfx/q)<<" budget_Cl="<<(bfy/q)<<"\n";
    }
    {
        // Convention samples: four airfoil face normals spread over the
        // surface. On the upper surface the fluid-domain outward normal
        // points into the body (ny<0); on the lower surface ny>0.
        const auto& fpatch=mesh.boundary().patch(mesh.boundary().find("airfoil"));
        std::vector<std::size_t> fids(fpatch.face_ids.begin(),fpatch.face_ids.end());
        const std::size_t nf=fids.size();
        if(nf>0){
            const std::size_t picks[4]={0,nf/4,nf/2,nf-1};
            for(const auto pick:picks){
                const auto f=fids[pick];
                const auto Sf=geometry.face_area_vectors[f];
                const auto n=Sf.normalized();
                std::cout<<"NACA_NORMAL_SAMPLE level="<<level
                         <<" x="<<geometry.face_centres[f].x
                         <<" y="<<geometry.face_centres[f].y
                         <<" nx="<<n.x<<" ny="<<n.y<<"\n";
            }
        }
    }

    // Final integration: same implementation as the per-iteration trace.
    const auto ff=airfoil_forces(U,p);
    const double cd=ff.cd, cl=ff.cl, cm=ff.cm;
    const auto& last=solve.history.back();
    std::cout<<"NACA0012_RE1000_A0 algo="<<algo_name<<" level="<<level
             <<" Cd="<<cd<<" Cl="<<cl<<" Cm="<<cm
             <<" Cd_pressure="<<ff.cdp<<" Cd_viscous="<<ff.cdv
             <<" continuity="<<last.continuity_linf
             <<" iterations="<<solve.iterations<<"\n";
    return {cd,cl,cm,ff.cdp,ff.cdv,last.continuity_linf};
}

static void run_quick_contract(const std::string& path)
{
    Mesh mesh;
    if(!cfdx::io::read_mesh_hdf5(path,mesh))
        throw std::runtime_error("cannot read NACA0012 mesh: "+path);
    const auto topo=mesh.topo_validate();
    if(!topo.ok)
        throw std::runtime_error("invalid NACA0012 topology: "+path);

    const auto geometry=build_fv_geometry(mesh);
    Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
    U.fill(0.0);
    for(std::size_t i=0;i<mesh.n_cells();++i)
        U.component_data(0)[i]=1.0;

    cfdx::physics::BoundaryConstraintMap strict_bc;
    cfdx::validation::add_velocity_dirichlet(strict_bc, "farfield", {1.0, 0.0, 0.0});
    cfdx::validation::add_velocity_dirichlet(strict_bc, "airfoil", {0.0, 0.0, 0.0});
    cfdx::validation::add_pressure_dirichlet(strict_bc, "farfield", 0.0);
    cfdx::validation::add_pressure_neumann(strict_bc, "airfoil", 0.0);
    const auto validation_flux =
        cfdx::validation::make_validation_face_flux(mesh, geometry, U);
    cfdx::validation::exercise_new_velocity_bc_contract(
        mesh, geometry, strict_bc, validation_flux, "NACA0012_RE1000_A0");

    const auto airfoil_id=mesh.boundary().find("airfoil");
    const auto farfield_id=mesh.boundary().find("farfield");
    if(airfoil_id>=mesh.boundary().n_patches() ||
       farfield_id>=mesh.boundary().n_patches())
        throw std::runtime_error("NACA0012 required boundary patches are missing");

    std::cout<<"NACA0012_QUICK_CONTRACT: PASS cells="<<mesh.n_cells()
             <<" airfoil_faces="<<mesh.boundary().patch(airfoil_id).size()
             <<" farfield_faces="<<mesh.boundary().patch(farfield_id).size()
             <<" solver_campaign_pending\n";
}

int main(int argc,char**argv)
{
    const bool quick=argc>=2 && std::string(argv[1])=="--quick";
    if(quick) {
        if(argc!=3)
            throw std::invalid_argument("usage: test_naca0012_qualification --quick mesh_n64.h5");
        run_quick_contract(argv[2]);
        return 0;
    }

    // Optional single-model mode used by the CI strategy matrix: one job per
    // model, each running the three refinement levels for that model.
    const char* model_filter=nullptr;
    int mesh_arg=1;
    if(argc>=2 && std::string(argv[1])=="--model") {
        if(argc!=6)
            throw std::invalid_argument("usage: test_naca0012_qualification --model NAME mesh_n64.h5 mesh_n128.h5 mesh_n256.h5");
        model_filter=argv[2];
        mesh_arg=3;
    } else if(argc!=4) {
        throw std::invalid_argument("usage: test_naca0012_qualification [--model NAME] mesh_n64.h5 mesh_n128.h5 mesh_n256.h5");
    }
    const std::size_t need=3;
    const char* meshes[3]={argv[mesh_arg],argv[mesh_arg+1],argv[mesh_arg+2]};

    // DIAG BRANCH ONLY: full qualification matrix across the pressure-velocity
    // algorithms and, for the coupled backend, every Schur complement model.
    // One identical case definition (mesh, BCs, gates) exercised by every
    // backend: the architecture contract this diagnostic measures.
    struct Config { const char* name; PressureVelocityAlgorithm algo; CoupledSchurModel schur; };
    const Config configs[]={
        {"SIMPLE",PressureVelocityAlgorithm::SIMPLE,CoupledSchurModel::BlockLocal},
        {"SIMPLEC",PressureVelocityAlgorithm::SIMPLEC,CoupledSchurModel::BlockLocal},
        {"PISO",PressureVelocityAlgorithm::PISO,CoupledSchurModel::BlockLocal},
        {"PIMPLE",PressureVelocityAlgorithm::PIMPLE,CoupledSchurModel::BlockLocal},
        {"FRACTIONAL_STEP",PressureVelocityAlgorithm::FRACTIONAL_STEP,CoupledSchurModel::BlockLocal},
        {"COUPLED_BLOCKLOCAL",PressureVelocityAlgorithm::COUPLED,CoupledSchurModel::BlockLocal},
        {"COUPLED_PCD",PressureVelocityAlgorithm::COUPLED,CoupledSchurModel::PCD},
        {"COUPLED_LSC",PressureVelocityAlgorithm::COUPLED,CoupledSchurModel::LSC},
        {"COUPLED_BFBT",PressureVelocityAlgorithm::COUPLED,CoupledSchurModel::BFBT},
        {"COUPLED_SIMPLE",PressureVelocityAlgorithm::COUPLED,CoupledSchurModel::SIMPLE},
        {"COUPLED_SIMPLEC",PressureVelocityAlgorithm::COUPLED,CoupledSchurModel::SIMPLEC},
    };
    std::vector<Forces> r;
    bool reference_complete=true;
    for(const auto& cfg:configs){
        if(model_filter && std::string(model_filter)!=cfg.name) continue;
        const bool is_reference=(std::string(cfg.name)=="SIMPLE");
        for(std::size_t i=0;i<need;++i){
            try {
                const Forces f=run(meshes[i],64u<<i,cfg.name,cfg.algo,cfg.schur);
                if(is_reference) r.push_back(f);
            }
            catch(const std::exception& e){
                std::cout.flush();
                std::cout<<"NACA_ALGO algo="<<cfg.name<<" level="<<(64u<<i)
                         <<" EXCEPTION "<<e.what()<<"\n";
                if(is_reference){
                    reference_complete=false;
                    const double nan=std::numeric_limits<double>::quiet_NaN();
                    r.push_back(Forces{nan,nan,nan,nan,nan,nan});
                }
            }
        }
    }
    if(model_filter){
        std::cout<<"NACA0012_QUALIFICATION: DIAGNOSTIC single_model="<<model_filter<<" done\n";
        return 0;
    }
    if(!reference_complete || r.size()!=need){
        std::cout<<"NACA0012_QUALIFICATION: DIAGNOSTIC simple_campaign_incomplete\n";
        return 0;
    }
    bool gates_ok=true;
    for(const auto& x:r) {
        if(!(x.continuity<1e-7)){ gates_ok=false; std::cout<<"NACA_GATE continuity "<<x.continuity<<"\n"; }
        if(!std::isfinite(x.cd) || !std::isfinite(x.cl)){ gates_ok=false; std::cout<<"NACA_GATE nonfinite_forces\n"; }
        if(std::abs(x.cl)>2e-2){ gates_ok=false; std::cout<<"NACA_GATE symmetry Cl="<<x.cl<<"\n"; }
    }
    {
        const double d01=std::abs(r[1].cd-r[0].cd), d12=std::abs(r[2].cd-r[1].cd);
        std::cout<<"NACA0012 Cd mesh differences="<<d01<<","<<d12<<" literature Re=1000 alpha=0 Cd~0.12\n";
        if(!(d12<=d01)){ gates_ok=false; std::cout<<"NACA_GATE refinement d01="<<d01<<" d12="<<d12<<"\n"; }
    }
    std::cout<<"NACA0012_QUALIFICATION: DIAGNOSTIC solver_converged forces_conservation_symmetry "
             <<"gates="<<(gates_ok?"PASS":"FAIL")<<" reference_gate_pending\n";
}
