#include "cfdx/physics/steady_incompressible_solver.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace cfdx::core;
using namespace cfdx::physics;
static void check(bool ok,const char* label) {if(!ok)throw std::runtime_error(label);}
static bool rejects(const std::function<void()>& call) {
    try {call();} catch(const std::exception&) {return true;}return false;
}
int main() {
    try {
        Mesh mesh;
        mesh.points().resize(8);
        const double xyz[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                                 {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
        for (std::size_t i=0;i<8;++i)
            mesh.points().set(i,xyz[i][0],xyz[i][1],xyz[i][2]);
        mesh.faces().push_face({0,3,2,1}); mesh.faces().push_face({4,5,6,7});
        mesh.faces().push_face({0,1,5,4}); mesh.faces().push_face({3,7,6,2});
        mesh.faces().push_face({0,4,7,3}); mesh.faces().push_face({1,2,6,5});
        mesh.ownership().resize(6);
        for (std::size_t f=0;f<6;++f) {
            mesh.ownership().set_owner(f,0);
            mesh.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
        }
        mesh.cells().push_cell({0,1,2,3,4,5});
        Patch inlet; inlet.name="inlet"; inlet.type=PatchType::INLET;
        inlet.face_ids={4}; mesh.boundary().add_patch(inlet);
        Patch outlet; outlet.name="outlet"; outlet.type=PatchType::OUTLET;
        outlet.face_ids={5}; mesh.boundary().add_patch(outlet);
        Patch aux; aux.name="aux_outlet"; aux.type=PatchType::OUTLET;
        aux.face_ids={3}; mesh.boundary().add_patch(aux);
        Patch wall; wall.name="wall"; wall.type=PatchType::WALL;
        wall.face_ids={0,1,2}; mesh.boundary().add_patch(wall);
        const auto geom=build_fv_geometry(mesh);
        Field<double,Location::FACE> old_phi(6,"phi","kg/s",1);
        old_phi.fill(0);
        VelocityBoundaryConditions ubc{
            {"inlet",{VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}}},
            {"outlet",{VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}}},
            {"aux_outlet",{VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}}},
            {"wall",{VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}}}};
        const std::map<std::string,Vec3> mixed{{"outlet",{0.3,-0.2,0.1}}};
        old_phi(5)=+0.25;
        auto out=resolve_inlet_outlet_velocity_faces(mesh,geom,old_phi,ubc,mixed);
        check(out[0].has(5),"outlet must be resolved");
        check(out[0].conditions[5].type == ScalarBoundaryFaceCondition::Type::FIXED_GRADIENT,
              "positive outward phi must select zero-gradient");
        check(out[0].conditions[4].type == ScalarBoundaryFaceCondition::Type::FIXED_VALUE,
              "inlet must retain fixed value");
        old_phi(5)=-0.25;
        auto in=resolve_inlet_outlet_velocity_faces(mesh,geom,old_phi,ubc,mixed);
        check(in[0].conditions[5].type == ScalarBoundaryFaceCondition::Type::FIXED_VALUE,
              "negative outward phi must select fixed inflow");
        check(std::abs(in[0].conditions[5].value-0.3)<1e-15 &&
              std::abs(in[1].conditions[5].value+0.2)<1e-15,
              "all velocity components must have matching inflow values");
        Field<double,Location::CELL> U(1,"U","m/s",3); U.fill(0.0);
        const auto out_flux=make_mass_flux(mesh,geom,U,1.0,ubc,{},false,&out);
        const auto in_flux=make_mass_flux(mesh,geom,U,1.0,ubc,{},false,&in);
        check(std::abs(out_flux(5))<1e-15,"zeroGradient outflow uses owner U");
        check(std::abs(in_flux(5)-0.3)<1e-15,"mixed inflow must override H/A face velocity");
        Field<double,Location::CELL> p(1,"p","Pa",1);p.fill(0.0);
        std::array<std::vector<double>,3> inv{{{1.0},{1.0},{1.0}}};
        ScalarBoundaryConditions pbc{{"outlet",{ScalarBoundaryType::FIXED_VALUE,0,0}},
                                     {"aux_outlet",{ScalarBoundaryType::ZERO_GRADIENT,0,0}},
                                     {"inlet",{ScalarBoundaryType::ZERO_GRADIENT,0,0}},
                                     {"wall",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}};
        auto rc=make_rhie_chow_mass_flux(mesh,geom,U,p,inv,1.0,ubc,pbc,
                  RhieChowVelocityContract::HBYA_PREDICTOR,false,false,&in);
        check(std::abs(rc(5)-in_flux(5))<1e-15,
              "zero pressure difference: RC and H/A must use identical mixed face BC");
        // Independent FV transport identity, including an outlet face that
        // reverses its direction. The assembled linear system and the physical
        // face-by-face flux must use the SAME per-face inlet/outlet mode.
        // This is deliberately not a NACA convergence or GCI certificate.
        Field<double,Location::CELL> scalar(1,"scalar","-",1);
        Field<double,Location::CELL> su(1,"Su","-",1);
        Field<double,Location::CELL> sp(1,"Sp","-",1);
        scalar(0)=0.7; su.fill(0.0); sp.fill(0.0);
        const ScalarBoundaryConditions scalar_bc{
            {"inlet",{ScalarBoundaryType::FIXED_VALUE,1.0,0.0}},
            {"outlet",{ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0}},
            {"aux_outlet",{ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0}},
            {"wall",{ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0}}};
        for (const double outlet_flux : {0.25, -0.25}) {
            old_phi.fill(0.0);
            old_phi(3)=1.25; // a separate free outlet keeps A_PP invertible
            old_phi(4)=-1.0;
            old_phi(5)=outlet_flux;
            const auto modes=resolve_inlet_outlet_velocity_faces(mesh,geom,old_phi,ubc,mixed);
            const auto matrix=assemble_scalar_equation(
                mesh,geom,old_phi,0.0,su,sp,scalar_bc,false,nullptr,nullptr,
                nullptr,nullptr,ConvectionScheme::UPWIND,&scalar,&modes[0]);
            const auto fluxes=reconstruct_scalar_transport_flux(
                mesh,geom,old_phi,scalar,0.0,scalar_bc,false,
                ConvectionScheme::UPWIND,&modes[0]);
            double face_residual=0.0;
            for (std::size_t f=0;f<mesh.n_faces();++f)
                face_residual+=fluxes(f);
            const double matrix_residual=matrix.matrix(0,0)*scalar(0)-matrix.rhs(0);
            const double expected=-1.0 + 1.25*0.7 + outlet_flux*(outlet_flux<0.0?0.3:0.7);
            check(std::abs(face_residual-matrix_residual)<1e-12,
                  "inletOutlet mismatch: assembled momentum / independent face flux");
            check(std::abs(face_residual-expected)<1e-12,
                  "inletOutlet wrong inflow/outflow scalar convective flux");
            // The bounded (advective) convective operator intentionally adds
            // the local -u*div(phi) term. Do not mistake its residual for the
            // conservative finite-volume face divergence when div(phi)!=0.
            const auto bounded=assemble_scalar_equation(
                mesh,geom,old_phi,0.0,su,sp,scalar_bc,true,nullptr,nullptr,
                nullptr,nullptr,ConvectionScheme::UPWIND,&scalar,&modes[0]);
            const double advective_residual=
                bounded.matrix(0,0)*scalar(0)-bounded.rhs(0);
            const double local_divergence=1.25-1.0+outlet_flux;
            check(std::abs(advective_residual-(face_residual-scalar(0)*local_divergence))<1e-12,
                  "bounded convection must differ from conservative budget only by -u div(phi)");
        }
        old_phi(5)=std::numeric_limits<double>::quiet_NaN();
        check(rejects([&]{(void)resolve_inlet_outlet_velocity_faces(mesh,geom,old_phi,ubc,mixed);}),
              "NaN phi must fail");
        old_phi(5)=0.0;
        check(rejects([&]{(void)resolve_inlet_outlet_velocity_faces(mesh,geom,old_phi,ubc,
              {{"unknown",{1,0,0}}});}),"unknown patch must fail");
        check(rejects([&]{(void)resolve_inlet_outlet_velocity_faces(mesh,geom,old_phi,ubc,
              {{"wall",{1,0,0}}});}),"fixed U patch cannot also be mixed");
        // A paired checkpoint contains accepted U,p,phi and residual references,
        // NOT the Anderson/NGMRES sample history. A restart claiming exact
        // multisecant continuity must fail rather than silently reset history.
        IncompressibleSolverControls restart_controls;
        restart_controls.convection_scheme = ConvectionScheme::TVD;
        restart_controls.initial_face_mass_flux_from_paired_checkpoint = true;
        restart_controls.initial_face_mass_flux.assign(mesh.n_faces(), 0.0);
        restart_controls.initial_momentum_residual_reference = 1.0;
        restart_controls.initial_nonlinear_metric_reference = 1.0;
        restart_controls.initial_accepted_iteration = 10;
        validate_incompressible_controls(restart_controls, mesh.n_cells());
        restart_controls.nonlinear_multisecant.method = NonlinearMultisecantMethod::ANDERSON;
        check(rejects([&]{validate_incompressible_controls(restart_controls, mesh.n_cells());}),
              "paired restart must reject missing Anderson nonlinear history");
        restart_controls.nonlinear_multisecant.method = NonlinearMultisecantMethod::NGMRES;
        check(rejects([&]{validate_incompressible_controls(restart_controls, mesh.n_cells());}),
              "paired restart must reject missing NGMRES nonlinear history");
        std::cout<<"INLET_OUTLET_VELOCITY_FACES: PASS\n";
    }catch(const std::exception& e) {
        std::cerr<<"INLET_OUTLET_VELOCITY_FACES: FAIL: "<<e.what()<<'\n'; return 1;
    }
}
