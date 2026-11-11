#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/sst_row_trace_v33.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
using namespace cfdx::testing;
int main() {
    run_case("v33_trace_parser_validates_before_simulation_mutation", []{
        const auto rows=sst_v33_trace_cells("341,329,425",712);
        EXPECT_TRUE(rows.size()==3 && rows[0]==341 && rows[1]==329 && rows[2]==425);
        EXPECT_TRUE(sst_v33_trace_cells(nullptr,712).empty());
        for (const char* invalid : {"", ",1", "1,", "1,,2", "1,01", "-1", "+1",
                                    "1,x", "712", "1844674407370955161612345",
                                    "0,1,2,3,4,5,6,7,8"})
            EXPECT_THROW(sst_v33_trace_cells(invalid,712),std::invalid_argument);
        EXPECT_THROW(sst_v33_trace_cells("1",4097),std::invalid_argument);
    });
    run_case("v33_branch_monitor_matches_real_sst_sources", []{
        TurbulenceTransportControls c;c.model=TurbulenceModel::SST;
        c.density=1.2;c.k_min=1e-12;c.omega_min=1e-12;
        for (double strain : {1.,10.,100.}) {
            for (double wall_g : {-1.,0.05}) {
                auto b=sst_v33_source_branches(
                    0.1,5.,strain,0.75,wall_g,c.density,c.a1,
                    c.beta_star,c.sst_production_limiter,c.k_min,c.omega_min);
                auto s=sst_cell_sources(0.1,5.,strain,0.7,0.75,0.,c,wall_g);
                EXPECT_NEAR(b.turbulent_viscosity,s.nut,1e-14);
                EXPECT_NEAR(b.limited_production,s.sk,1e-14);
            }
        }
        auto capped=sst_v33_source_branches(1.,2.,100.,1.,-1.,1.,.31,.09,10.,1e-10,1e-10);
        EXPECT_TRUE(capped.production_limited);
        auto free=sst_v33_source_branches(1.,2.,1.,1.,-1.,1.,.31,.09,10.,1e-10,1e-10);
        EXPECT_TRUE(!free.production_limited);
        EXPECT_TRUE(sst_v33_branch_changed(capped,free));
    });
    run_case("v33_branch_switches_and_floors_are_diagnostic_only", []{
        auto a=sst_v33_source_branches(1.,1.,.31,1.,-1.,1.,.31,.09,10.,1e-10,1e-10);
        EXPECT_TRUE(a.close_to_viscosity_switch);
        auto b=sst_v33_source_branches(1.,1.,.311,1.,-1.,1.,.31,.09,10.,1e-10,1e-10);
        EXPECT_TRUE(sst_v33_branch_changed(a,b));
        auto floored=sst_v33_source_branches(1e-11,1e-11,1.,1.,-1.,1.,.31,.09,10.,1e-10,1e-10);
        EXPECT_TRUE(floored.k_on_floor && floored.omega_on_floor);
        EXPECT_THROW(sst_v33_source_branches(std::numeric_limits<double>::quiet_NaN(),1.,1.,1.,-1.,1.,.31,.09,10.,1e-10,1e-10),std::invalid_argument);
    });
    run_case("v33_semi_smooth_k_source_derivative_switches_sign", []{
        TurbulenceTransportControls c;c.model=TurbulenceModel::SST;
        c.density=1.;c.k_min=1e-12;c.omega_min=1e-12;
        const auto check=[](const auto& state,double k,double w,double strain,
                            const auto& c) {
            const double h=1e-5;
            auto fun=[&](double kk,double ww) {
                const auto source=sst_cell_sources(
                    kk,ww,strain,1.,1.,0.,c,-1.);
                return c.density*c.beta_star*kk*ww-source.sk;
            };
            const double dw=(fun(k,w*std::exp(h))-fun(k,w*std::exp(-h)))/(2*h);
            const double dk=(fun(k*std::exp(h),w)-fun(k*std::exp(-h),w))/(2*h);
            EXPECT_NEAR(state.d_rk_dlog_k_frozen_source,dk,1e-6);
            EXPECT_NEAR(state.d_rk_dlog_omega_frozen_source,dw,1e-6);
        };
        const auto capped=sst_v33_source_branches(1.,2.,100.,1.,-1.,
            c.density,c.a1,c.beta_star,c.sst_production_limiter,c.k_min,c.omega_min);
        const auto uncapped=sst_v33_source_branches(1.,2.,0.1,1.,-1.,
            c.density,c.a1,c.beta_star,c.sst_production_limiter,c.k_min,c.omega_min);
        EXPECT_TRUE(capped.frozen_source_derivative_valid);
        EXPECT_TRUE(uncapped.frozen_source_derivative_valid);
        EXPECT_TRUE(capped.d_rk_dlog_omega_frozen_source < 0.);
        EXPECT_TRUE(uncapped.d_rk_dlog_omega_frozen_source > 0.);
        check(capped,1.,2.,100.,c);
        check(uncapped,1.,2.,0.1,c);
        const auto near_switch=sst_v33_source_branches(1.,1.,.31,1.,-1.,
            c.density,c.a1,c.beta_star,c.sst_production_limiter,c.k_min,c.omega_min);
        EXPECT_TRUE(!near_switch.frozen_source_derivative_valid);
    });
    run_case("v33_ptc_attribution_is_not_an_outer_commit", []{
        SSTV30ContinuationState state;
        state.accept_inner(1.,0.8);
        state.record_selected_source(false);
        state.accept_inner(1.,0.7);
        state.record_selected_source(true);
        EXPECT_TRUE(state.accepted_inner==2);
        EXPECT_TRUE(state.selected_own==1 && state.selected_other==1);
        EXPECT_TRUE(state.accepted_outer==0);
        // Outer rejection rolls back the state and all proposed CFL changes.
        const auto original=SSTV30ContinuationState{};
        state=original;
        EXPECT_TRUE(state.selected_own==0 && state.selected_other==0);
        EXPECT_NEAR(state.cfl,0.5,0.);
    });
    return run_all();
}
