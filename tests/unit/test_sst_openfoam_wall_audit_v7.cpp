#include "../audit/chain_fixture.h"
#include "cfdx/physics/sst_openfoam_wall_audit_v7.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
#include <limits>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("wall_reference_audits_all_actual_faces_without_mutation", [] {
        auto mesh=make_1d_chain(3);
        const auto geometry=build_fv_geometry(mesh);
        Field<double,Location::CELL> U(3,"U","m/s",3), k(3,"k","m2/s2",1);
        U.fill(0.0); k.fill(0.05);
        for (std::size_t i=0; i<3; ++i) U(i,0)=2.0;
        TurbulenceTransportControls ctl;
        ctl.model=TurbulenceModel::SST;
        ctl.molecular_viscosity=1e-5;
        const auto before=make_sst_wall_omega_face_values(mesh,geometry,{"wall"},ctl);
        const std::map<std::string,std::array<double,3>> wall_u{{"wall",{0,0,0}}};
        const auto audit=audit_openfoam_omega_wall_reference(
            mesh,geometry,{"wall"},before,U,k,wall_u,ctl);
        EXPECT_TRUE(audit.face_count==mesh.boundary().patch(2).face_ids.size());
        EXPECT_TRUE(audit.face_count==12);
        EXPECT_TRUE(audit.viscous_faces+audit.buffer_faces+audit.logarithmic_faces==12);
        EXPECT_TRUE(std::isfinite(audit.worst_relative_omega_gap));
        EXPECT_TRUE(audit.worst_relative_omega_gap>=0 && audit.worst_relative_omega_gap<=1);
        EXPECT_TRUE(audit.maximum_k_based_y_plus>=audit.minimum_k_based_y_plus);
        EXPECT_TRUE(before.values.at("wall")[audit.worst_face]>0);
    });
    run_case("wall_audit_rejects_uncomparable_or_corrupted_face", [] {
        auto mesh=make_1d_chain(1);
        const auto geometry=build_fv_geometry(mesh);
        Field<double,Location::CELL> U(1,"U","m/s",3), k(1,"k","",1);
        U.fill(0.0); k.fill(0.01);
        TurbulenceTransportControls ctl;
        ctl.model=TurbulenceModel::SST;
        ctl.molecular_viscosity=1e-5;
        auto faces=make_sst_wall_omega_face_values(mesh,geometry,{"wall"},ctl);
        const std::map<std::string,std::array<double,3>> wall_u{{"wall",{0,0,0}}};
        EXPECT_THROW(audit_openfoam_omega_wall_reference(
            mesh,geometry,{"wall"},faces,U,k,{},ctl),std::invalid_argument);
        faces.values.at("wall")[mesh.boundary().patch(2).face_ids.front()]=
            std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(audit_openfoam_omega_wall_reference(
            mesh,geometry,{"wall"},faces,U,k,wall_u,ctl),std::domain_error);
    });
    return run_all();
}
