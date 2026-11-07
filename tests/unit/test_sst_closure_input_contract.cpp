#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/turbulence_transport.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>
#include <stdexcept>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("sst_wls_planar_empty_requires_2d_gradients", [] {
        cfdx::core::Mesh mesh;
        mesh.boundary().add_patch(cfdx::core::Patch{
            "front", cfdx::core::PatchType::EMPTY, {}, {}});
        TurbulenceTransportControls controls;
        controls.model = TurbulenceModel::SST;
        controls.sst_use_wls_boundary_gradients = true;
        controls.sst_gradient_dimension = 3;
        EXPECT_THROW_WITH(require_sst_empty_wls_dimension(mesh, controls),
                          std::invalid_argument, "dimension=2");
        controls.sst_gradient_dimension = 2;
        require_sst_empty_wls_dimension(mesh, controls);
        controls.sst_use_wls_boundary_gradients = false;
        controls.sst_gradient_dimension = 3;
        require_sst_empty_wls_dimension(mesh, controls);
    });
    run_case("healthy_sst_source_and_eddy_viscosity", [] {
        TurbulenceTransportControls c;
        c.model = TurbulenceModel::SST;
        const auto q=sst_cell_sources(2.0,5.0,3.0,0.3,0.8,0.0,c);
        EXPECT_NEAR(q.nut, 0.31*2.0/(3.0*0.8), 1e-14);
        EXPECT_NEAR(turbulence_nu_t(2.0,5.0,3.0,0.1,c,-1.0,0.8),q.nut,1e-14);
    });
    run_case("nonfinite_blending_factors_rejected", [] {
        TurbulenceTransportControls c;
        c.model = TurbulenceModel::SST;
        const double nan=std::numeric_limits<double>::quiet_NaN();
        const double inf=std::numeric_limits<double>::infinity();
        for (const double invalid : {nan,inf,-inf}) {
            EXPECT_THROW_WITH(sst_cell_sources(2.0,5.0,3.0,invalid,0.8,0.0,c),
                              std::invalid_argument,"SST");
            EXPECT_THROW_WITH(sst_cell_sources(2.0,5.0,3.0,0.3,invalid,0.0,c),
                              std::invalid_argument,"SST");
            EXPECT_THROW_WITH(turbulence_nu_t(2.0,5.0,3.0,0.1,c,-1.0,invalid),
                              std::invalid_argument,"SST");
        }
    });
    run_case("nonphysical_sst_blending_and_strain_rejected", [] {
        TurbulenceTransportControls c;
        c.model = TurbulenceModel::SST;
        for (const double invalid_f : {-0.01, 1.01}) {
            EXPECT_THROW_WITH(sst_cell_sources(2.0,5.0,3.0,invalid_f,0.8,0.0,c),
                              std::invalid_argument,"SST");
            EXPECT_THROW_WITH(sst_cell_sources(2.0,5.0,3.0,0.3,invalid_f,0.0,c),
                              std::invalid_argument,"SST");
        }
        EXPECT_THROW_WITH(sst_cell_sources(2.0,5.0,-1.0,0.3,0.8,0.0,c),
                          std::invalid_argument,"SST");
        EXPECT_THROW_WITH(turbulence_nu_t(2.0,5.0,-1.0,0.1,c,-1.0,0.8),
                          std::invalid_argument,"SST");
    });
    run_case("sst_eddy_viscosity_rejects_nonfinite_inputs", [] {
        const double nan=std::numeric_limits<double>::quiet_NaN();
        const double inf=std::numeric_limits<double>::infinity();
        for (const double invalid : {nan, inf, -inf}) {
            EXPECT_THROW_WITH(sst_eddy_viscosity(invalid,5.0,3.0),
                              std::invalid_argument,"SST");
            EXPECT_THROW_WITH(sst_eddy_viscosity(2.0,invalid,3.0),
                              std::invalid_argument,"SST");
            EXPECT_THROW_WITH(sst_eddy_viscosity(2.0,5.0,invalid),
                              std::invalid_argument,"SST");
            EXPECT_THROW_WITH(sst_eddy_viscosity(2.0,5.0,3.0,invalid),
                              std::invalid_argument,"SST");
        }
    });
    run_case("nonfinite_sst_strain_and_cross_gradient_rejected", [] {
        TurbulenceTransportControls c;
        c.model = TurbulenceModel::SST;
        const double nan=std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW_WITH(sst_cell_sources(2.0,5.0,nan,0.3,0.8,0.0,c),
                          std::invalid_argument,"SST");
        EXPECT_THROW_WITH(sst_cell_sources(2.0,5.0,3.0,0.3,0.8,nan,c),
                          std::invalid_argument,"SST");
        EXPECT_THROW_WITH(turbulence_nu_t(2.0,5.0,nan,0.1,c,-1.0,0.8),
                          std::invalid_argument,"SST");
    });
    return run_all();
}
