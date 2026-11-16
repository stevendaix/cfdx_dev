#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"
#include <array>
#include <cmath>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {
Mesh skewed_two_cell_mesh() {
    Mesh mesh;
    mesh.points().resize(12);
    const double xyz[12][3] = {
        {0,0,0},{2,0,0},{3,0,0},{0,1,0},{2,1,0},{3,1,0},
        {0,0,1},{2,0,1},{3,0,1},{0,1,1},{2,1,1},{3,1,1}};
    for (std::size_t i=0;i<12;++i) {
        const double shift=(i==1 || i==4 || i==7 || i==10) ? 0.15 :
            ((i==2 || i==5 || i==8 || i==11) ? 0.55 : 0.0);
        mesh.points().set(i,xyz[i][0],xyz[i][1]+shift,xyz[i][2]);
    }
    for (const auto& face : std::vector<std::vector<FaceIndex>>{
        {0,6,9,3},{0,1,7,6},{3,9,10,4},{0,3,4,1},{6,7,10,9},
        {7,10,4,1},{2,5,11,8},{1,2,8,7},{4,10,11,5},
        {1,4,5,2},{7,8,11,10}}) mesh.faces().push_face(face);
    mesh.ownership().resize(11);
    for (std::size_t f=0;f<11;++f) {
        mesh.ownership().set_owner(f,f<=5?0:1);
        mesh.ownership().set_neighbour(f,f==5?1:FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0,1,2,3,4,5});
    mesh.cells().push_cell({5,6,7,8,9,10});
    Patch left,right;
    left.name="left";left.type=PatchType::WALL;left.face_ids={0,1,2,3,4};
    right.name="right";right.type=PatchType::WALL;right.face_ids={6,7,8,9,10};
    mesh.boundary().add_patch(left);
    mesh.boundary().add_patch(right);
    return mesh;
}
}

int main() {
    run_case("full_rhie_chow_compensated_step_matches_fp64_on_skewed_grid", [] {
        const auto mesh=skewed_two_cell_mesh();
        const auto geom=build_fv_geometry(mesh);
        const std::array<std::vector<double>,3> rAU{{{0.2,0.3},{0.4,0.5},{0.6,0.8}}};
        ScalarBoundaryConditions pbcs;
        pbcs["left"]={ScalarBoundaryType::FIXED_VALUE,3.0,0.0};
        pbcs["right"]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        Field<double,Location::CELL> p(2,"p","Pa",1);
        Field<double,Location::CELL> H(2,"H","m/s",3);
        p(0)=0.7;p(1)=1.9;
        H.set(0,3.0,0.1,0.0);H.set(1,2.0,-0.5,0.0);
        Vector high(2,0.0),low(2,0.0);
        high(0)=0.29;high(1)=-0.61;
        low(0)=3.0e-13;low(1)=-1.0e-13;
        for (const bool projected: {false,true}) {
            const auto original=make_rhie_chow_mass_flux(
                mesh,geom,H,p,rAU,1.2,{},pbcs,
                RhieChowVelocityContract::HBYA_PREDICTOR,false,projected);
            const auto split=reconstruct_compensated_rhie_chow_pressure_step(
                mesh,geom,original,high,low,rAU,1.2,pbcs,projected);
            auto updated=p;
            for (std::size_t i=0;i<2;++i)
                updated(i)+=high(i)+low(i);
            const auto direct=make_rhie_chow_mass_flux(
                mesh,geom,H,updated,rAU,1.2,{},pbcs,
                RhieChowVelocityContract::HBYA_PREDICTOR,false,projected);
            for(std::size_t f=0;f<mesh.n_faces();++f)
                EXPECT_NEAR(static_cast<double>(split[f]),direct(f),1e-11);
        }
    });
    run_case("low_component_survives_when_fp64_rounded_pressure_loses_it", [] {
        const auto mesh=skewed_two_cell_mesh();
        const auto geom=build_fv_geometry(mesh);
        const std::array<std::vector<double>,3> rAU{{{0.2,0.3},{0.4,0.5},{0.6,0.8}}};
        Field<double,Location::CELL> H(2,"H","m/s",3);
        Field<double,Location::CELL> p(2,"p","Pa",1);
        H.fill(0.0);p.fill(0.0);
        const auto base=make_rhie_chow_mass_flux(
            mesh,geom,H,p,rAU,1.2,{}, {},
            RhieChowVelocityContract::HBYA_PREDICTOR,false,false);
        Vector high(2,1.0e8),low(2,0.0);
        const auto only_high=reconstruct_compensated_rhie_chow_pressure_step(
            mesh,geom,base,high,low,rAU,1.2,{},false);
        low(1)=1.0e-9;
        const auto true_split=reconstruct_compensated_rhie_chow_pressure_step(
            mesh,geom,base,high,low,rAU,1.2,{},false);
        // FP64 cannot represent this low part next to 1e8, but the internal
        // pressure-face response must preserve it.
        EXPECT_TRUE(1.0e8 + low(1) == 1.0e8);
        EXPECT_TRUE(std::abs(true_split[5]-only_high[5])>1.0e-12L);
        EXPECT_TRUE(std::isfinite(true_split[5]));
    });
    run_case("canonical_high_low_full_flux_matches_direct_fp64_on_skewed_grid", [] {
        const auto mesh=skewed_two_cell_mesh();
        const auto geom=build_fv_geometry(mesh);
        const std::array<std::vector<double>,3> rAU{{{0.2,0.3},{0.4,0.5},{0.6,0.8}}};
        ScalarBoundaryConditions bcs;
        bcs["left"]={ScalarBoundaryType::FIXED_VALUE,3.0,0.0};
        Field<double,Location::CELL> p(2,"p","Pa",1), H(2,"H","m/s",3);
        p(0)=0.7;p(1)=1.9;
        H.set(0,3.0,0.1,0.0);H.set(1,2.0,-0.5,0.0);
        Vector low(2,0.0);
        for (const bool projected: {false,true}) {
            const auto actual=make_rhie_chow_mass_flux_high_low(
                mesh,geom,H,p,low,rAU,1.2,{},bcs,projected);
            const auto rounded=make_rhie_chow_mass_flux(
                mesh,geom,H,p,rAU,1.2,{},bcs,
                RhieChowVelocityContract::HBYA_PREDICTOR,false,projected);
            for (std::size_t f=0;f<mesh.n_faces();++f)
                EXPECT_NEAR(static_cast<double>(actual[f]),rounded(f),1e-11);
            const auto gradient=gauss_gradient_pressure_high_low(
                mesh,geom,p,low,bcs,projected);
            const auto reference=gauss_gradient_with_boundary(
                p,mesh,geom,bcs,false,projected);
            for(std::size_t c=0;c<mesh.n_cells();++c)
                for(std::size_t d=0;d<3;++d)
                    EXPECT_NEAR(gradient.component_data(d)[c],
                                reference.component_data(d)[c],1e-11);
        }
    });
    run_case("physical_pressure_low_is_retained_by_velocity_relaxation", [] {
        Field<double,Location::CELL> p(2,"p","Pa",1);
        p(0)=1.0e8;p(1)=-3.0e8;
        Vector low(2,0.0),dhi(2,0.0),dlo(2,0.0);
        dhi(0)=1.0e-9;dhi(1)=-2.0e-9;
        add_pressure_high_low(p,low,dhi,dlo,1.0);
        EXPECT_TRUE(p(0)==1.0e8);
        EXPECT_TRUE(p(1)==-3.0e8);
        EXPECT_TRUE(std::abs(low(0)-1.0e-9)<1e-12);
        EXPECT_TRUE(std::abs(low(1)+2.0e-9)<1e-12);
        // Low pressure affects the flux operator even though FP64 p did not move.
        const auto mesh=skewed_two_cell_mesh();
        const auto geom=build_fv_geometry(mesh);
        Field<double,Location::CELL> H(2,"H","m/s",3);H.fill(0.0);
        const std::array<std::vector<double>,3> rAU{{{0.2,0.3},{0.4,0.5},{0.6,0.8}}};
        const auto base=make_rhie_chow_mass_flux_high_low(
            mesh,geom,H,p,Vector(2,0.0),rAU,1.0,{},{});
        const auto exact=make_rhie_chow_mass_flux_high_low(
            mesh,geom,H,p,low,rAU,1.0,{},{});
        EXPECT_TRUE(std::abs(exact[5]-base[5])>1e-12L);
    });
    run_case("true_face_residual_orientation_and_signed_conservation", [] {
        const auto mesh=skewed_two_cell_mesh();
        std::vector<long double> flux(mesh.n_faces(),0.0L);
        flux[5]=3.0L;
        Vector rhs(2,3.0);
        const auto result=pressure_true_face_residual(mesh,flux,rhs);
        EXPECT_NEAR(result.first(0),-3.0,1e-15);
        EXPECT_NEAR(result.first(1),3.0,1e-15);
        EXPECT_NEAR(result.second,1.0,1e-15);
    });
    run_case("true_face_refinement_rejects_wrong_preconditioner", [] {
        Vector hi(2,0.0),lo(2,0.0);
        const auto evaluator=[](const Vector& h,const Vector& l) {
            Vector residual(2,0.0);
            residual(0)=1.0-h(0)-l(0);
            residual(1)=2.0-h(1)-l(1);
            return std::make_pair(residual,
                std::hypot(residual(0),residual(1))/std::sqrt(5.0));
        };
        const auto malicious=[](const Vector& rhs, Vector& delta) {
            for(std::size_t i=0;i<delta.size();++i) delta(i)=-rhs(i);
            return true;
        };
        const auto fail=refine_pressure_true_face_high_low(
            hi,lo,1e-11,4,evaluator,malicious);
        EXPECT_TRUE(!fail.converged);
        EXPECT_TRUE(fail.accepted==0);
        EXPECT_NEAR(hi(0),0.0,1e-15);
        const auto good=[](const Vector& rhs,Vector& delta) {
            for(std::size_t i=0;i<delta.size();++i) delta(i)=rhs(i);
            return true;
        };
        const auto pass=refine_pressure_true_face_high_low(
            hi,lo,1e-11,3,evaluator,good);
        EXPECT_TRUE(pass.converged);
        EXPECT_TRUE(pass.accepted==1);
        EXPECT_NEAR(hi(0),1.0,1e-15);
        EXPECT_NEAR(hi(1),2.0,1e-15);
    });
    run_case("split_pressure_low_contributes_to_momentum_gradient", [] {
        const auto mesh=skewed_two_cell_mesh();
        const auto geom=build_fv_geometry(mesh);
        Field<double,Location::CELL> p(2,"p","Pa",1);
        p(0)=1.0e8;p(1)=1.0e8;
        Vector zero(2,0.0),low(2,0.0);
        low(1)=1.0e-9;
        const auto old_gradient=gauss_gradient_pressure_high_low(
            mesh,geom,p,zero,{},false);
        const auto new_gradient=gauss_gradient_pressure_high_low(
            mesh,geom,p,low,{},false);
        bool changed=false;
        for (std::size_t c=0;c<mesh.n_cells();++c)
            for (std::size_t d=0;d<3;++d)
                if (std::abs(new_gradient.component_data(d)[c]-
                             old_gradient.component_data(d)[c])>1.0e-12)
                    changed=true;
        EXPECT_TRUE(changed);
    });
    run_case("pressure_refinement_globalizes_overshooting_direction", [] {
        Vector hi(1,0.0), lo(1,0.0);
        const auto physical=[](const Vector& h,const Vector& l) {
            Vector residual(1,1.0-h(0)-l(0));
            return std::make_pair(residual,std::abs(residual(0)));
        };
        const auto oversized=[](const Vector& rhs,Vector& delta) {
            delta(0)=2.0*rhs(0);
            return true;
        };
        const auto report=refine_pressure_true_face_high_low(
            hi,lo,1.0e-11,3,physical,oversized);
        EXPECT_TRUE(report.converged);
        EXPECT_TRUE(report.attempted==1 && report.accepted==1);
        EXPECT_NEAR(hi(0)+lo(0),1.0,1.0e-15);
    });
    run_case("invalid_production_high_low_contract_is_rejected", [] {
        IncompressibleSolverControls c;
        c.pressure_high_low_production=true;
        bool refused=false;
        try { validate_incompressible_controls(c,2); }
        catch (const std::invalid_argument&) { refused=true; }
        EXPECT_TRUE(refused);
    });
    return run_all();
}
