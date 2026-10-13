#include "cfdx/physics/turbulence_transport.h"
#include "cfdx/physics/turbulence_models.h"
#include "cfdx/physics/turbulence_solver.h"
#include "cfdx/physics/radiation_advanced.h"
#include "cfdx/physics/energy_solver.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "common/test_harness.h"
#include <limits>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh one_cell_mesh() {
    Mesh m;
    m.points().resize(8);
    const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                          {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1}); m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4}); m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3}); m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0;f<6;++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    for (std::size_t f = 0; f < 6; ++f) wall.face_ids.push_back(f);
    m.boundary().add_patch(wall);
    return m;
}

int main() {
    // --- N11 F2: turbulence_nu_t read-level floor reporting ---
    run_case("nu_t_without_report_keeps_historical_clamp", [] {
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::KEPSILON;
        const double clamped=turbulence_nu_t(-1.0,-1.0,0.0,0.0,c);
        const double floored=turbulence_nu_t(c.k_min,c.epsilon_min,0.0,0.0,c);
        EXPECT_NEAR(clamped,floored,1e-24);
    });
    run_case("nu_t_report_counts_k_and_second_floors", [] {
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::KEPSILON;
        TurbulenceBoundDiagnostics d;
        turbulence_nu_t(1e-14,1e-14,0.0,0.0,c,&d);
        EXPECT_TRUE(d.first_floor_activations==1);
        EXPECT_TRUE(d.second_floor_activations==1);
        EXPECT_TRUE(d.first_total_correction>0.0);
        EXPECT_TRUE(d.second_total_correction>0.0);
        EXPECT_TRUE(d.first_max_correction>0.0);
        EXPECT_TRUE(d.second_max_correction>0.0);
        EXPECT_TRUE(d.nonfinite_values()==0);
    });
    run_case("nu_t_report_counts_nonfinite_without_clamping", [] {
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::KEPSILON;
        TurbulenceBoundDiagnostics d;
        const double nut=turbulence_nu_t(
            std::numeric_limits<double>::quiet_NaN(),c.epsilon_min,0.0,0.0,c,&d);
        EXPECT_TRUE(std::isnan(nut));
        EXPECT_TRUE(d.first_nonfinite==1);
        EXPECT_TRUE(d.first_floor_activations==0);
    });
    run_case("nu_t_report_omega_based_second_floor", [] {
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::KOMEGA;
        TurbulenceBoundDiagnostics d;
        turbulence_nu_t(c.k_min,1e-14,0.0,0.0,c,&d);
        EXPECT_TRUE(d.second_floor_activations==1);
        EXPECT_TRUE(d.second_total_correction>0.0);
    });
    run_case("nu_t_positive_inputs_report_nothing", [] {
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::KEPSILON;
        TurbulenceBoundDiagnostics d;
        turbulence_nu_t(1.0,10.0,0.0,0.0,c,&d);
        EXPECT_TRUE(d.first_floor_activations==0);
        EXPECT_TRUE(d.second_floor_activations==0);
        EXPECT_TRUE(d.nonfinite_values()==0);
        EXPECT_TRUE(d.first_total_correction==0.0);
    });

    // --- N11 F2: compute_kepsilon_production_field read-level floors ---
    run_case("production_field_report_counts_floors_and_keeps_values", [] {
        const Mesh m=one_cell_mesh();
        Field<double,Location::CELL> k(m.n_cells(),"k","m2/s2",1);
        Field<double,Location::CELL> eps(m.n_cells(),"epsilon","m2/s3",1);
        Field<double,Location::CELL> strain(m.n_cells(),"S","1/s",1);
        Field<double,Location::CELL> p_ref(m.n_cells(),"Pref","W/m3",1);
        Field<double,Location::CELL> p(m.n_cells(),"P","W/m3",1);
        k(0)=1e-14; eps(0)=1e-14; strain(0)=3.0;
        TurbulenceTransportControls c;
        c.model=TurbulenceModel::KEPSILON;
        compute_kepsilon_production_field(m,strain,k,eps,p_ref,c);
        TurbulenceBoundDiagnostics d;
        compute_kepsilon_production_field(m,strain,k,eps,p,c,&d);
        EXPECT_NEAR(p(0),p_ref(0),1e-18);
        EXPECT_TRUE(d.first_floor_activations==1);
        EXPECT_TRUE(d.second_floor_activations==1);
        EXPECT_TRUE(d.nonfinite_values()==0);
    });

    // --- N11 F3: production limiter diagnostics ---
    run_case("production_limiter_diagnostics_count_activation", [] {
        TurbulenceCorrectionControls c;
        TurbulenceInvariants x;
        TurbulenceProductionDiagnostics d;
        const double limited=corrected_turbulence_production(50.0,x,c,1.0,&d);
        const double reference=corrected_turbulence_production(50.0,x,c);
        EXPECT_NEAR(limited,c.production_limit,1e-15);
        EXPECT_NEAR(limited,reference,1e-15);
        EXPECT_TRUE(d.production_limiter_activations==1);
        EXPECT_NEAR(d.production_limiter_total_reduction,
                    50.0-c.production_limit,1e-15);
        EXPECT_NEAR(d.production_limiter_max_reduction,
                    50.0-c.production_limit,1e-15);
    });
    run_case("production_limiter_diagnostics_silent_below_limit", [] {
        TurbulenceCorrectionControls c;
        TurbulenceInvariants x;
        TurbulenceProductionDiagnostics d;
        const double p=corrected_turbulence_production(5.0,x,c,1.0,&d);
        EXPECT_NEAR(p,5.0,1e-15);
        EXPECT_TRUE(d.production_limiter_activations==0);
        EXPECT_TRUE(d.production_limiter_total_reduction==0.0);
    });
    run_case("production_limiter_disabled_reports_nothing", [] {
        TurbulenceCorrectionControls c;
        c.production_limiter=false;
        TurbulenceInvariants x;
        TurbulenceProductionDiagnostics d;
        const double p=corrected_turbulence_production(50.0,x,c,1.0,&d);
        EXPECT_NEAR(p,50.0,1e-15);
        EXPECT_TRUE(d.production_limiter_activations==0);
    });

    // --- N11 F-R2: Rosseland minimum-temperature floor reporting ---
    run_case("rosseland_minimum_temperature_floor_is_counted", [] {
        Mesh m=one_cell_mesh();
        auto g=build_fv_geometry(m);
        Field<double,Location::FACE> mass_flux(m.n_faces(),"phi","kg/s",1);
        mass_flux.fill(0.0);
        Field<double,Location::CELL> temperature(m.n_cells(),"T","K",1);
        Field<double,Location::CELL> source(m.n_cells(),"S","W/m3",1);
        Field<double,Location::CELL> absorption(m.n_cells(),"kappa","1/m",1);
        temperature(0)=300.0;
        source(0)=-1.0e5;
        absorption(0)=1.0;
        EnergySolverControls energy;
        energy.density=1.0; energy.cp=1000.0; energy.conductivity=1.0; energy.dt=0.0;
        RosselandSolveControls controls;
        controls.max_iterations=3;
        controls.minimum_temperature=1.0;
        ScalarBoundaryConditions bc;
        bc["wall"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
        const auto r=solve_rosseland_energy(m,g,mass_flux,temperature,source,
                                            absorption,energy,controls,bc);
        EXPECT_TRUE(r.minimum_temperature_floor_activations>0);
        EXPECT_TRUE(r.minimum_temperature_total_correction>0.0);
        EXPECT_TRUE(r.minimum_temperature_max_correction>0.0);
        EXPECT_NEAR(temperature(0),controls.minimum_temperature,1e-12);
        // Healthy nonlinear solve: no non-finite candidate was seen by the
        // acceptance clamp, so the dedicated counter must stay at zero.
        EXPECT_TRUE(r.nonfinite_candidates==0);
    });

    // --- N11 F-R2 follow-up (audit of #864): the Rosseland floor must pin
    // the exact clamp contract for non-finite candidates. The historical
    // std::max(minimum_temperature, candidate) converts NaN (and -inf)
    // into the finite floor -- that conversion must be observable through
    // the helper flags instead of remaining indistinguishable from a
    // healthy iteration. ---
    run_case("rosseland_floor_helper_pins_nonfinite_contract", [] {
        const double minimum=1.0;
        // Finite candidate below the floor: measurable correction.
        {
            const auto f=rosseland_apply_minimum_temperature_floor(0.5,minimum);
            EXPECT_TRUE(f.floored && !f.nonfinite && !f.nonfinite_floored);
            EXPECT_TRUE(f.bounded==minimum);
        }
        // Finite candidate exactly at the floor: untouched, not an
        // activation (strictly-below contract, as in #864).
        {
            const auto f=rosseland_apply_minimum_temperature_floor(minimum,minimum);
            EXPECT_TRUE(!f.floored && !f.nonfinite);
            EXPECT_TRUE(f.bounded==minimum);
        }
        // Finite candidate above the floor: untouched.
        {
            const auto f=rosseland_apply_minimum_temperature_floor(350.0,minimum);
            EXPECT_TRUE(!f.floored && !f.nonfinite);
            EXPECT_TRUE(f.bounded==350.0);
        }
        // NaN candidate: std::max(minimum, NaN) returns the floor, so the
        // non-finite candidate silently became a finite temperature. The
        // helper must expose that conversion without counting it as a
        // measurable floor activation.
        {
            const auto f=rosseland_apply_minimum_temperature_floor(
                std::numeric_limits<double>::quiet_NaN(),minimum);
            EXPECT_TRUE(f.nonfinite && f.nonfinite_floored && !f.floored);
            EXPECT_TRUE(f.bounded==minimum);
        }
        // -inf candidate: also converted into the finite floor.
        {
            const auto f=rosseland_apply_minimum_temperature_floor(
                -std::numeric_limits<double>::infinity(),minimum);
            EXPECT_TRUE(f.nonfinite && f.nonfinite_floored && !f.floored);
            EXPECT_TRUE(f.bounded==minimum);
        }
        // +inf candidate: propagates unclamped (documented historical
        // behaviour) -- non-finite but NOT converted to the floor.
        {
            const auto f=rosseland_apply_minimum_temperature_floor(
                std::numeric_limits<double>::infinity(),minimum);
            EXPECT_TRUE(f.nonfinite && !f.nonfinite_floored && !f.floored);
            EXPECT_TRUE(std::isinf(f.bounded) && f.bounded>0.0);
        }
    });

    return run_all();
}
