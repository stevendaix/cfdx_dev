#include "../audit/chain_fixture.h"
#include "cfdx/physics/sst_openfoam_wall_coupling_v8.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>

using namespace cfdx::physics;
using namespace cfdx::testing;
int main() {
    run_case("spalding_viscous_limit_and_monotone_wall_nut", [] {
        const double nu=1e-5, y=1e-4;
        EXPECT_NEAR(spalding_wall_nut_v8(0.,y,nu),0.,0.);
        const double slow=spalding_wall_nut_v8(0.05,y,nu);
        const double fast=spalding_wall_nut_v8(20.,y,nu);
        EXPECT_TRUE(std::isfinite(slow) && slow>=0.);
        EXPECT_TRUE(std::isfinite(fast) && fast>=slow);
    });
    run_case("openfoam_wall_viscous_branch_preserves_g_proxy", [] {
        const double g=0.123;
        const auto r=openfoam_wall_coupling_point_v8(0.,1e-5,1e-4,0.,g);
        EXPECT_NEAR(r.production,g,1e-14);
        EXPECT_NEAR(r.nut_spalding,0.,0.);
        EXPECT_TRUE(std::isfinite(r.omega) && r.omega>0.);
    });
    run_case("openfoam_wall_log_branch_changes_G_and_nut", [] {
        const auto r=openfoam_wall_coupling_point_v8(0.09,1e-5,1e-3,5.,0.);
        EXPECT_TRUE(r.yplus_from_k>11.06);
        EXPECT_TRUE(r.nut_spalding>0.);
        EXPECT_TRUE(std::isfinite(r.production) && r.production>0.);
    });
    run_case("actual_faces_wall_triple_audit_is_diagnostic_only", [] {
        auto mesh=make_1d_chain(3);
        const auto geom=build_fv_geometry(mesh);
        SSTProductionState st(3);
        st.k.fill(.05);st.omega.fill(2.);
        st.nu_t.fill(.005);st.strain_rate.fill(2.);st.wall_distance.fill(.3);
        Field<double,Location::CELL> U(3,"U","m/s",3);
        U.fill(0.);for(std::size_t i=0;i<3;++i)U(i,0)=2.;
        TurbulenceTransportControls c; c.model=TurbulenceModel::SST;
        c.molecular_viscosity=1e-5;
        const auto audit=audit_openfoam_wall_coupling_v8(
            mesh,geom,{"wall"},U,st,{{"wall",{0.,0.,0.}}},c);
        EXPECT_TRUE(audit.faces==12);
        EXPECT_TRUE(audit.max_omega_owner_gap>0.);
        EXPECT_TRUE(audit.max_k_yplus>0.);
        EXPECT_TRUE(!audit.openfoam_equivalent);
    });
    run_case("v8_wall_audit_rejects_invalid_geometry", [] {
        EXPECT_THROW(spalding_wall_nut_v8(1.,0.,1e-5),std::invalid_argument);
        EXPECT_THROW(openfoam_wall_coupling_point_v8(0.1,0.,1e-3,1.,1.),
                     std::invalid_argument);
    });
    return run_all();
}
