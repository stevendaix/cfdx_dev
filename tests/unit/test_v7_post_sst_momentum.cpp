#include "../audit/sst_channel_fixture.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("post_sst_gate_rejects_missing_sst_runtime", [] {
        SSTChannelCase c;
        c.controls.sst_runtime=nullptr;
        c.controls.require_post_sst_momentum_reassembly=true;
        EXPECT_THROW_WITH(
            solve_steady_incompressible(c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls),
            std::invalid_argument,"SST");
    });
    run_case("post_sst_reassembles_actual_production_carrier_operator", [] {
        SSTChannelCase c;
        c.controls.require_post_sst_momentum_reassembly=true;
        c.runtime.max_transport_iterations=30;
        const auto result=solve_steady_incompressible(
            c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls);
        EXPECT_TRUE(result.history.size()==1);
        const auto& h=result.history.back();
        EXPECT_TRUE(std::isfinite(h.post_sst_momentum_rhs_relative));
        EXPECT_TRUE(h.post_sst_momentum_rhs_relative>=0.0);
        EXPECT_TRUE(std::isfinite(h.momentum_equation_residual_rhs_active_component_max));
        std::cout << "post_sst=" << h.post_sst_momentum_rhs_relative
                  << " pre_sst=" << h.momentum_equation_residual_rhs_active_component_max
                  << " sst_converged=" << result.converged << "\n";
    });
    run_case("post_sst_gate_probes_changed_viscosity_in_shear_fixture", [] {
        SSTChannelCase c;
        c.controls.require_post_sst_momentum_reassembly=true;
        c.ubc["wall"].value={0.0,0.0,0.0};
        const auto result=solve_steady_incompressible(
            c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls);
        EXPECT_TRUE(result.history.size()==1);
        const auto& h=result.history.back();
        EXPECT_TRUE(std::isfinite(h.post_sst_momentum_rhs_relative));
        EXPECT_TRUE(h.post_sst_momentum_rhs_relative>=0.0);
        EXPECT_TRUE(std::abs(h.post_sst_momentum_rhs_relative-
            h.momentum_equation_residual_rhs_active_component_max)>1e-8);
        std::cout << "post_sst_shear=" << h.post_sst_momentum_rhs_relative
                  << " pre_sst_shear=" << h.momentum_equation_residual_rhs_active_component_max
                  << " nu_t=" << c.state.nu_t(0) << "\n";
    });
    return run_all();
}
