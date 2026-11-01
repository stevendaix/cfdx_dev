#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "qualification_boundary_helpers.h"
#include "naca_momentum_trace.h"
#include "naca_relaxation_controls.h"
#include "cfdx/physics/wall_face_traction_reconstruction.h"
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

struct Forces { double cd,cl,cm,cdp,cdv,continuity; };

static Forces run(const std::string& path,std::size_t level)
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
    const bool split_farfield = mesh.boundary().find("farfield_inlet") <
                                mesh.boundary().n_patches() &&
                                mesh.boundary().find("farfield_outlet") <
                                mesh.boundary().n_patches();
    VelocityBoundaryConditions ubc;
    if (split_farfield) {
        ubc["farfield_inlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,
                                  {Uinf,0,0}};
        ubc["farfield_outlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,
                                   {0,0,0}};
    } else {
        ubc["farfield"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{Uinf,0,0}};
    }
    ubc["airfoil"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
    ubc["front"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ubc["back"]={VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
    ScalarBoundaryConditions pbc;
    if (split_farfield) {
        pbc["farfield_inlet"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
        pbc["farfield_outlet"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    } else {
        pbc["farfield"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    }
    pbc["airfoil"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    pbc["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    pbc["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};

    IncompressibleSolverControls c;
    c.algorithm=PressureVelocityAlgorithm::SIMPLE;
    c.density=rho; c.kinematic_viscosity=nu;
    c.linear_max_iterations=4000; c.linear_tolerance=1e-10;
    c.convergence.max_iterations=8000;
    c.convergence.relative_tolerance=1e-8;
    c.convergence.continuity_tolerance=1e-8;
    c.coupling.alpha_u=0.6; c.coupling.alpha_p=0.25;
    cfdx::validation::apply_naca_relaxation_ab(c);
    // A/B only: the default NACA qualification is unchanged.
    if (const char* value = std::getenv("CFDX_NACA_SKEW_GRADIENT")) {
        if (std::string(value) != "1")
            throw std::invalid_argument("CFDX_NACA_SKEW_GRADIENT must be 1");
        c.pressure_skew_corrected_gradient = true;
    }
    // Experimental A/B only: the normal production and CI path is unchanged.
    if (const char* line_search = std::getenv("CFDX_NACA_PRESSURE_LINE_SEARCH")) {
        if (std::string(line_search) != "1")
            throw std::invalid_argument("CFDX_NACA_PRESSURE_LINE_SEARCH must be 1");
        c.pressure_continuity_line_search.enabled = true;
    }

    // A/B only: TVD is an existing bounded CFDX convection scheme. Keep the
    // unchanged UPWIND reference for all default NACA qualification gates.
    c.use_bounded_convection=true;
    c.convection_scheme=ConvectionScheme::UPWIND;
    if (const char* tvd = std::getenv("CFDX_NACA_CONVECTION_TVD")) {
        if (std::string(tvd) != "1")
            throw std::invalid_argument("CFDX_NACA_CONVECTION_TVD must be 1");
        c.convection_scheme=ConvectionScheme::TVD;
    }
    if (const char* affine = std::getenv("CFDX_NACA_AFFINE_FLUX_TRIALS")) {
        if (std::string(affine) != "1")
            throw std::invalid_argument("CFDX_NACA_AFFINE_FLUX_TRIALS must be 1");
        // The fast path requires the globally enabled line-search.
        if (!c.pressure_continuity_line_search.enabled)
            throw std::invalid_argument("affine flux trials require line-search");
        c.pressure_continuity_line_search.affine_flux_trials = true;
    }
    c.pressure_reference_cell=0;
    c.pressure_reference_value=0.0;
    c.diagnostics.iteration_trace=true;
    c.diagnostics.iteration_trace_frequency=200;

    const auto geometry=build_fv_geometry(mesh);
    cfdx::validation::NacaMomentumTrace naca_trace(mesh);
    if (const char* forensic = std::getenv("CFDX_NACA_FORENSICS");
        forensic != nullptr && std::string(forensic) == "1") {
        c.momentum_assembly_forensic_callback =
            [&naca_trace](std::size_t iteration, const auto& pressure,
                          const auto& grad_pressure, const auto& velocity,
                          const auto& mass_flux, const auto& mu_eff) {
                naca_trace(iteration, pressure, grad_pressure, velocity,
                           mass_flux, mu_eff);
            };
        c.momentum_equation_forensic_callback =
            [&naca_trace](std::size_t iteration, std::size_t component,
                          const auto& equation, const auto& geometry,
                          const auto& grad_pressure, const auto& body_force) {
                naca_trace.equation(iteration, component, equation, geometry,
                                    grad_pressure, body_force);
            };
    }
    for(std::size_t i=0;i<mesh.n_cells();++i)
        U.component_data(0)[i]=Uinf;

    // Exercise the strict mathematical BC/FVM path before the production solve.
    // The legacy maps below are retained for the current steady-solver compatibility API.
    cfdx::physics::BoundaryConstraintMap strict_bc;
    if (split_farfield)
        cfdx::validation::add_velocity_dirichlet(
            strict_bc, "farfield_inlet", {Uinf, 0.0, 0.0});
    else
        cfdx::validation::add_velocity_dirichlet(
            strict_bc, "farfield", {Uinf, 0.0, 0.0});
    if (split_farfield) {
        cfdx::validation::add_velocity_flux_dependent(
            strict_bc, "farfield_outlet", {Uinf,0.0,0.0}, {0.0,0.0,0.0});
        cfdx::validation::add_pressure_neumann(strict_bc, "farfield_inlet", 0.0);
    }
    cfdx::validation::add_velocity_dirichlet(strict_bc, "airfoil", {0.0, 0.0, 0.0});
    cfdx::validation::add_pressure_dirichlet(
        strict_bc, split_farfield ? "farfield_outlet" : "farfield", 0.0);
    cfdx::validation::add_pressure_neumann(strict_bc, "airfoil", 0.0);
    const auto validation_flux = cfdx::validation::make_validation_face_flux(mesh, geometry, U);
    cfdx::validation::exercise_new_velocity_bc_contract(
        mesh, geometry, strict_bc, validation_flux, "NACA0012_RE1000_A0");

    const auto solve=solve_steady_incompressible(mesh,U,p,ubc,pbc,c);
    if(!solve.converged) {
        std::ostringstream message;
        message << "NACA0012 solver did not converge: status="
                << cfdx::core::to_string(solve.convergence_status)
                << " reason=" << solve.convergence_reason
                << " iterations=" << solve.iterations;
        if(!solve.history.empty()) {
            const auto& last=solve.history.back();
            message << " momentum_rhs_relative="
                    << last.momentum_equation_residual_rhs_relative
                    << " momentum_reference_relative="
                    << last.momentum_equation_residual_relative
                    << " pressure_residual=" << last.pressure_residual
                    << " continuity=" << last.continuity_normalized
                    << " dU=" << last.velocity_change_inf
                    << " dP=" << last.pressure_change_inf
                    << " alpha_u=" << last.effective_alpha_u
                    << " alpha_p=" << last.effective_alpha_p;
        }
        throw std::runtime_error(message.str());
    }

    ScalarBoundaryConditions gx,gy;
    if (split_farfield) {
        gx["farfield_inlet"]={ScalarBoundaryType::FIXED_VALUE,Uinf,0};
        gx["farfield_outlet"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    } else gx["farfield"]={ScalarBoundaryType::FIXED_VALUE,Uinf,0};
    gx["airfoil"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    gx["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    gx["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    if (split_farfield) {
        gy["farfield_inlet"]={ScalarBoundaryType::FIXED_VALUE,0,0};
        gy["farfield_outlet"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    } else gy["farfield"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    gy["airfoil"]={ScalarBoundaryType::FIXED_VALUE,0,0};
    gy["front"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    gy["back"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
    Field<double,Location::CELL> ux(mesh.n_cells(),"ux","m/s",1), uy(mesh.n_cells(),"uy","m/s",1);
    for(std::size_t i=0;i<mesh.n_cells();++i){ux(i)=U.component_data(0)[i];uy(i)=U.component_data(1)[i];}
    const auto gux=gauss_gradient_with_boundary(ux,mesh,geometry,gx);
    const auto guy=gauss_gradient_with_boundary(uy,mesh,geometry,gy);
    // Optional face-traction A/B; no modification to the historical force
    // gate.  This evaluation remains experimental until MMS + grid studies.
    bool face_traction = false;
    if (const char* mode = std::getenv("CFDX_NACA_WALL_TRACTION_FACE")) {
        if (std::string(mode) != "1")
            throw std::invalid_argument("CFDX_NACA_WALL_TRACTION_FACE must be 1");
        face_traction = true;
    }
    Field<double, Location::CELL> wall_pressure_gradient;
    if (face_traction)
        wall_pressure_gradient = gauss_gradient_with_boundary(
            p, mesh, geometry, pbc, c.pressure_skew_corrected_gradient);

    const auto& patch=mesh.boundary().patch(mesh.boundary().find("airfoil"));
    std::ofstream face_profiles;
    if (const char* output = std::getenv("CFDX_NACA_FACE_PROFILES_CSV")) {
        if (!face_traction || std::string(output).empty())
            throw std::invalid_argument("NACA Cp/Cf profiles require face traction and a path");
        face_profiles.open(output);
        if (!face_profiles)
            throw std::runtime_error("cannot open NACA face profiles CSV");
        face_profiles << std::setprecision(17);
        face_profiles << "face,x_over_c,y_over_c,Cp,Cf_tangent,Cd_pressure_face,Cd_viscous_face\n";
    }
    double Fx=0,Fy=0,Mz=0,Fxp=0,Fyp=0,Fxv=0,Fyv=0;
    const double q=0.5*rho*Uinf*Uinf;
    for(const auto f:patch.face_ids){
        const std::size_t o=mesh.ownership().owner(f);
        const auto Sf=geometry.face_area_vectors[f];
        const auto n=Sf.normalized();
        double pp=p(o);
        double txx=2.0*rho*nu*gux.component_data(0)[o];
        double tyy=2.0*rho*nu*guy.component_data(1)[o];
        double txy=rho*nu*(gux.component_data(1)[o]+guy.component_data(0)[o]);
        if (face_traction) {
            const auto wall = reconstruct_wall_face_traction(
                geometry.face_centres[f] - geometry.cell_centres[o], n,
                p(o),
                {wall_pressure_gradient.component_data(0)[o],
                 wall_pressure_gradient.component_data(1)[o],
                 wall_pressure_gradient.component_data(2)[o]},
                {ux(o), uy(o), 0.0}, {0.0,0.0,0.0},
                {{{gux.component_data(0)[o],gux.component_data(1)[o],
                   gux.component_data(2)[o]},
                  {guy.component_data(0)[o],guy.component_data(1)[o],
                   guy.component_data(2)[o]},
                  {0.0,0.0,0.0}}}, rho*nu);
            pp = wall.pressure;
            // Existing total-force integration uses the traction components
            // of the symmetric stress tensor contracted with n.
            // Reconstruct stress components by evaluating the wall gradient;
            // the returned traction is the physically authoritative vector.
            const double face_txv = wall.viscous_on_fluid_normal.x;
            const double face_tyv = wall.viscous_on_fluid_normal.y;
            const double tx = (pp*n.x-face_txv)*Sf.mag();
            const double ty = (pp*n.y-face_tyv)*Sf.mag();
            Fx += tx; Fy += ty;
            Fxp += pp*n.x*Sf.mag(); Fyp += pp*n.y*Sf.mag();
            Fxv -= face_txv*Sf.mag(); Fyv -= face_tyv*Sf.mag();
            const auto r = geometry.face_centres[f]-Vec3{0.25,0,0};
            Mz += r.x*ty-r.y*tx;
            if (face_profiles) {
                Vec3 tangent{-n.y,n.x,0.0};
                if (tangent.x < 0.0) tangent *= -1.0;
                const double cf_tangent =
                    -(face_txv*tangent.x + face_tyv*tangent.y) / q;
                face_profiles << f << ',' << geometry.face_centres[f].x/chord
                              << ',' << geometry.face_centres[f].y/chord
                              << ',' << pp/q << ',' << cf_tangent
                              << ',' << (pp*n.x*Sf.mag())/q
                              << ',' << (-face_txv*Sf.mag())/q << '\n';
            }
            continue;
        }
        // Face normal convention: n = Sf/|Sf| points out of the fluid domain
        // (into the body) on wall faces, because build_fv_geometry()
        // enforces Sf.(Cf - C_owner) > 0 on boundary faces. With this
        // fluid-outward n, the force exerted by the fluid on the body is
        //   F = int (p n - tau n) dS
        // i.e. pressure pushes along +p n and viscous traction acts along
        // -tau n (equivalent to the body->fluid normal contract of
        // src/cfdx/physics/forces.h).
        const double tx=(pp*n.x-txx*n.x-txy*n.y)*Sf.mag();
        const double ty=(pp*n.y-txy*n.x-tyy*n.y)*Sf.mag();
        Fx+=tx; Fy+=ty;
        Fxp+=pp*n.x*Sf.mag(); Fyp+=pp*n.y*Sf.mag();
        Fxv+=-(txx*n.x+txy*n.y)*Sf.mag();
        Fyv+=-(txy*n.x+tyy*n.y)*Sf.mag();
        const auto r=geometry.face_centres[f]-Vec3{0.25,0,0};
        Mz+=r.x*ty-r.y*tx;
    }
    if (face_profiles && !face_profiles.good())
        throw std::runtime_error("failed to write NACA face profiles CSV");
    const double cd=Fx/q, cl=Fy/q, cm=Mz/q;
    const auto& last=solve.history.back();
    std::cout<<"NACA0012_RE1000_A0 level="<<level
             <<" Cd="<<cd<<" Cl="<<cl<<" Cm="<<cm
             <<" Cd_pressure="<<Fxp/q<<" Cd_viscous="<<Fxv/q
             <<" continuity="<<last.continuity_linf
             <<" iterations="<<solve.iterations<<"\n";
    return {cd,cl,cm,Fxp/q,Fxv/q,last.continuity_linf};
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
    const auto split_in = mesh.boundary().find("farfield_inlet");
    const auto split_out = mesh.boundary().find("farfield_outlet");
    const bool split_farfield = split_in < mesh.boundary().n_patches() &&
                                split_out < mesh.boundary().n_patches();
    cfdx::validation::add_velocity_dirichlet(
        strict_bc, split_farfield ? "farfield_inlet" : "farfield", {1.0, 0.0, 0.0});
    if (split_farfield) {
        cfdx::validation::add_velocity_flux_dependent(
            strict_bc, "farfield_outlet", {1.0,0.0,0.0}, {0.0,0.0,0.0});
        cfdx::validation::add_pressure_neumann(strict_bc, "farfield_inlet", 0.0);
    }
    cfdx::validation::add_velocity_dirichlet(strict_bc, "airfoil", {0.0, 0.0, 0.0});
    cfdx::validation::add_pressure_dirichlet(
        strict_bc, split_farfield ? "farfield_outlet" : "farfield", 0.0);
    cfdx::validation::add_pressure_neumann(strict_bc, "airfoil", 0.0);
    const auto validation_flux =
        cfdx::validation::make_validation_face_flux(mesh, geometry, U);
    cfdx::validation::exercise_new_velocity_bc_contract(
        mesh, geometry, strict_bc, validation_flux, "NACA0012_RE1000_A0");

    const auto airfoil_id=mesh.boundary().find("airfoil");
    const auto farfield_id=mesh.boundary().find("farfield");
    if(airfoil_id>=mesh.boundary().n_patches() ||
       (!split_farfield && farfield_id>=mesh.boundary().n_patches()))
        throw std::runtime_error("NACA0012 required boundary patches are missing");

    std::cout<<"NACA0012_QUICK_CONTRACT: PASS cells="<<mesh.n_cells()
             <<" airfoil_faces="<<mesh.boundary().patch(airfoil_id).size()
             <<" farfield_faces="<<(split_farfield ?
                 mesh.boundary().patch(split_in).size() +
                 mesh.boundary().patch(split_out).size() :
                 mesh.boundary().patch(farfield_id).size())
             <<" solver_campaign_pending\n";
}

int main(int argc,char**argv)
{
    try {
    const bool quick=argc>=2 && std::string(argv[1])=="--quick";
    if(quick) {
        if(argc!=3)
            throw std::invalid_argument("usage: test_naca0012_qualification --quick mesh_n64.h5");
        run_quick_contract(argv[2]);
        return 0;
    }

    // Isolated replay for one refinement level. This does not qualify
    // mesh convergence; it preserves the per-level physical acceptance gates.
    if(argc==3 && std::string(argv[1])=="--single") {
        const Forces x=run(argv[2],0);
        if(!(x.continuity<1e-7))
            throw std::runtime_error("NACA0012 continuity gate failed");
        if(!std::isfinite(x.cd) || !std::isfinite(x.cl))
            throw std::runtime_error("NACA0012 force integration produced non-finite QoI");
        if(std::abs(x.cl)>2e-2)
            throw std::runtime_error("NACA0012 alpha=0 symmetry gate failed");
        if(!(x.cdp>0.0))
            throw std::runtime_error("NACA0012 pressure drag sign gate failed");
        std::cout<<"NACA0012_SINGLE: PASS (one level only; not a mesh qualification)\n";
        return 0;
    }

    const std::size_t need=3;
    if(argc!=4)
        throw std::invalid_argument("usage: test_naca0012_qualification [--single mesh.h5 | mesh_n64.h5 mesh_n128.h5 mesh_n256.h5]");
    std::vector<Forces> r;
    for(std::size_t i=0;i<need;++i) r.push_back(run(argv[i+1],64u<<i));
    for(const auto& x:r) {
        if(!(x.continuity<1e-7)) throw std::runtime_error("NACA0012 continuity gate failed");
        if(!std::isfinite(x.cd) || !std::isfinite(x.cl)) throw std::runtime_error("NACA0012 force integration produced non-finite QoI");
        if(std::abs(x.cl)>2e-2) throw std::runtime_error("NACA0012 alpha=0 symmetry gate failed");
        // Regression gate for the wall-normal orientation of the pressure
        // force integration: at alpha=0 the airfoil pressure drag is
        // strictly positive. This catches the historical sign inversion
        // where Fxp accumulated -p n with n pointing into the body.
        if(!(x.cdp>0.0)) throw std::runtime_error("NACA0012 pressure drag sign gate failed (Cd_pressure must be > 0 at alpha=0)");
    }
    if(!quick) {
        const double d01=std::abs(r[1].cd-r[0].cd), d12=std::abs(r[2].cd-r[1].cd);
        std::cout<<"NACA0012 Cd mesh differences="<<d01<<","<<d12<<" literature Re=1000 alpha=0 Cd~0.12\n";
        if(!(d12<=d01)) throw std::runtime_error("NACA0012 drag is not stabilizing under refinement");
    }
    std::cout<<"NACA0012_QUALIFICATION: DIAGNOSTIC solver_converged forces_conservation_symmetry_reference_gate_pending\n";
    } catch (const std::exception& error) {
        std::cerr << "NACA0012_QUALIFICATION: FAIL: " << error.what() << '\n';
        return 1;
    }
}
