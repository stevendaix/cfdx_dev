#include "../audit/chain_fixture.h"
#include "cfdx/physics/sst_openfoam_wall_triplet_v15.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
using namespace cfdx::physics;
using namespace cfdx::testing;
int main() {
    run_case("openfoam_wall_triplet_three_fields_one_state", [] {
        auto mesh=make_1d_chain(3);
        const auto geometry=build_fv_geometry(mesh);
        SSTProductionState state(3);
        state.k.fill(0.05);
        state.omega.fill(2.0);
        state.nu_t.fill(0.005);
        state.strain_rate.fill(2.0);
        state.wall_distance.fill(0.3);
        Field<double,Location::CELL> velocity(3,"U","m/s",3);
        velocity.fill(0.0);
        for (std::size_t i=0;i<3;++i) velocity(i,0)=2.0;
        TurbulenceTransportControls controls;
        controls.model=TurbulenceModel::SST;
        controls.molecular_viscosity=1e-5;
        auto prepared=prepare_openfoam_wall_triplet_v15(
            mesh,geometry,{"wall"},velocity,state,{{"wall",{0.0,0.0,0.0}}},controls);
        EXPECT_TRUE(prepared.wall_faces==12);
        EXPECT_TRUE(prepared.wall_cells>=1);
        EXPECT_TRUE(!prepared.openfoam_equivalent && !prepared.applied_to_momentum &&
                    !prepared.applied_to_sst);
        std::size_t nf=0;
        for (std::size_t f=0;f<mesh.n_faces();++f) {
            if (!prepared.has_wall_face[f]) continue;
            ++nf;
            EXPECT_TRUE(std::isfinite(prepared.omega_by_face[f]) &&
                        prepared.omega_by_face[f]>0.0);
            EXPECT_TRUE(std::isfinite(prepared.nut_by_face[f]) &&
                        prepared.nut_by_face[f]>=0.0);
        }
        EXPECT_TRUE(nf==12);
        for (std::size_t i=0;i<mesh.n_cells();++i)
            if (prepared.has_wall_cell[i])
                EXPECT_TRUE(std::isfinite(prepared.production_by_cell[i]) &&
                            prepared.production_by_cell[i]>=0.0);
        EXPECT_THROW(prepare_openfoam_wall_triplet_v15(
            mesh,geometry,{"wall","wall"},velocity,state,
            {{"wall",{0.0,0.0,0.0}}},controls),std::invalid_argument);
        EXPECT_THROW(prepare_openfoam_wall_triplet_v15(
            mesh,geometry,{"wall"},velocity,state,{},controls),std::invalid_argument);
    });
    return run_all();
}
