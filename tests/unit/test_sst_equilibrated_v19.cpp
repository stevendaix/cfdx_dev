#include "../audit/chain_fixture.h"
#include "cfdx/physics/sst_wall_operator_v17.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
#include <limits>
using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("v19_equilibrated_defect_on_unchanged_original_operator", [] {
        ScalarEquation equation;
        cfdx::core::SparseMatrix matrix(3,3);
        matrix.push_back(0,0,1.0);
        matrix.push_back(1,1,1.e-6);
        matrix.push_back(2,2,1.0);
        matrix.finalize();
        equation.matrix=matrix;
        equation.rhs=cfdx::core::Vector(3,0.0);
        equation.rhs(0)=1.e7;
        equation.rhs(1)=0.001;
        equation.rhs(2)=10.0;
        equation.diagonal={1.,1.e-6,1.};
        const auto original=equation;
        cfdx::core::Vector solution(3,0.0);
        solution(0)=1.e7;
        solution(1)=1200.;
        solution(2)=10.;
        const auto r=sst_equilibrated_defect_correction_v19(equation,solution,1.e-8);
        EXPECT_TRUE(r.attempted);
        EXPECT_TRUE(r.accepted);
        EXPECT_TRUE(r.after<r.before*1.e-4);
        EXPECT_NEAR(solution(1),1000.,1.e-6);
        EXPECT_NEAR(equation.matrix(1,1),original.matrix(1,1),1.e-18);
        EXPECT_NEAR(equation.rhs(1),original.rhs(1),1.e-18);
        const auto repeat=sst_equilibrated_defect_correction_v19(equation,solution,1.e-8);
        EXPECT_TRUE(!repeat.attempted);
    });
    run_case("v19_rejects_missing_original_diagonal", [] {
        ScalarEquation equation;
        cfdx::core::SparseMatrix matrix(1,1);
        matrix.push_back(0,0,1.0);
        matrix.finalize();
        equation.matrix=matrix;
        equation.rhs=cfdx::core::Vector(1,1.0);
        cfdx::core::Vector solution(1,0.0);
        bool rejected=false;
        try { (void)sst_equilibrated_defect_correction_v19(equation,solution,1e-8); }
        catch (const std::invalid_argument&) { rejected=true; }
        EXPECT_TRUE(rejected);
    });
    run_case("v19_zero_defect_preserves_field", [] {
        ScalarEquation eq;
        cfdx::core::SparseMatrix a(2,2);
        a.push_back(0,0,1.);a.push_back(1,1,1.);a.finalize();
        eq.matrix=a;eq.rhs=cfdx::core::Vector(2,1.);
        eq.diagonal={1.,1.};
        cfdx::core::Vector x(2,1.);
        const auto diag=sst_equilibrated_defect_correction_v19(eq,x,1e-8);
        EXPECT_TRUE(!diag.attempted && !diag.accepted);
        EXPECT_NEAR(x(0),1.,1e-14);
        EXPECT_NEAR(x(1),1.,1e-14);
    });
    run_case("v19_invalid_inputs_are_not_silent", [] {
        ScalarEquation eq;
        cfdx::core::SparseMatrix a(2,2);
        a.push_back(0,0,0.);a.push_back(1,1,1.);a.finalize();
        eq.matrix=a;eq.rhs=cfdx::core::Vector(2,1.);
        eq.diagonal={0.,1.};
        cfdx::core::Vector x(2,0.);
        EXPECT_THROW(sst_equilibrated_defect_correction_v19(eq,x,1e-8),std::invalid_argument);
    });
    run_case("v19_openfoam_corner_face_count_vs_legacy_area_weight", [] {
        const auto mesh=make_1d_chain(3);
        const auto geometry=build_fv_geometry(mesh);
        SSTProductionState state(mesh.n_cells());
        state.k.fill(0.05);state.omega.fill(2.);
        state.nu_t.fill(0.005);state.strain_rate.fill(2.);
        state.wall_distance.fill(0.3);
        cfdx::core::Field<double,cfdx::core::Location::CELL> U(mesh.n_cells(),"U","m/s",3);
        U.fill(0.);
        for (std::size_t i=0;i<mesh.n_cells();++i) U(i,0)=2.0;
        TurbulenceTransportControls ctl;
        ctl.model=TurbulenceModel::SST;
        ctl.molecular_viscosity=1e-5;ctl.density=1.;
        const std::map<std::string,std::array<double,3>> fixed{
            {"inlet",{0.,0.,0.}},{"wall",{0.,0.,0.}}};
        const auto area=prepare_openfoam_wall_triplet_v15(
            mesh,geometry,{"inlet","wall"},U,state,fixed,ctl);
        const auto equal=prepare_openfoam_wall_triplet_v15(
            mesh,geometry,{"inlet","wall"},U,state,fixed,ctl,true);
        EXPECT_TRUE(area.corner_cells>0 && equal.corner_cells>0);
        // Explicitly verify the face-count policy, not just a difference from area.
        long double omega_sum=0.0L;
        std::size_t count=0;
        for (std::size_t f=0;f<mesh.n_faces();++f)
            if (equal.has_wall_face[f] && mesh.ownership().owner(f)==0) {
                omega_sum+=equal.omega_by_face[f];
                ++count;
            }
        EXPECT_TRUE(count>1);
        EXPECT_NEAR(equal.omega_by_cell[0],static_cast<double>(omega_sum/count),1e-10);
        EXPECT_TRUE(std::isfinite(area.omega_by_cell[0]));
        EXPECT_TRUE(std::isfinite(equal.omega_by_cell[0]));
    });
    return run_all();
}
