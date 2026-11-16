#include "cfdx/physics/anisotropic_fv_oracle_v33.h"
#include "common/test_harness.h"
#include <vector>
using namespace cfdx::physics;
using namespace cfdx::testing;
int main(){
    run_case("v33_QR_anisotropy_1e5_in_2D_extruded", []{
        const std::vector<SSTV33Vec3> offsets{
            {1e-5,1,0},{-2e-5,1,0},{3e-5,-2,0},{-1e-5,-2,0},
            {4e-5,3,0},{-4e-5,-3,0}};
        std::vector<SSTV33WlsSample> samples;
        for(const auto d:offsets) samples.push_back({d,3*d.x-2*d.y,1.});
        const auto g=sst_v33_scaled_qr_gradient(samples,2);
        EXPECT_NEAR(g.gradient.x,3.,1e-8);
        EXPECT_NEAR(g.gradient.y,-2.,1e-11);
        EXPECT_NEAR(g.gradient.z,0.,0.);
        EXPECT_THROW(sst_v33_scaled_qr_gradient(samples,3),std::domain_error);
    });
    run_case("v33_QR_explicit_rank_failure_and_no_silent_zero", []{
        std::vector<SSTV33WlsSample> samples{
            {{1,0,0},1,1},{{2,0,0},2,1},{{3,0,0},3,1}};
        EXPECT_THROW(sst_v33_scaled_qr_gradient(samples,2),std::domain_error);
    });
    run_case("v33_diffusion_affine_exact_and_conservative_on_skew_face", []{
        const SSTV33Vec3 d{2.,1.,0.},area{1.5,.5,0.},grad{.25,.75,0.};
        const auto phi_neighbor=v33_dot(d,grad);
        const auto flux=sst_v33_conservative_diffusion(d,area,grad,0.,phi_neighbor,2.);
        EXPECT_NEAR(flux.owner_flux,-2.*v33_dot(grad,area),1e-13);
        EXPECT_NEAR(flux.owner_flux+flux.neighbour_flux,0.,0.);
        EXPECT_TRUE(flux.two_point_conductance>0.);
    });
    run_case("v33_diffusion_rejects_inverted_and_invalid_face", []{
        EXPECT_THROW(sst_v33_conservative_diffusion({1,0,0},{-1,0,0},
            {1,0,0},0,1,1),std::domain_error);
        EXPECT_THROW(sst_v33_conservative_diffusion({1,0,0},{0,1,0},
            {1,0,0},0,1,1),std::domain_error);
    });
    return run_all();
}
