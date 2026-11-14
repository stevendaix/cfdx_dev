#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/sst_pseudo_transient_v29.h"
#include "cfdx/physics/sst_openfoam_weighted_owner_v29.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {
ScalarEquation two_cell(double b0=2.,double b1=8.) {
    ScalarEquation q;
    q.matrix=cfdx::core::SparseMatrix(2,2);
    q.matrix.push_back(0,0,4.);
    q.matrix.push_back(0,1,-1.);
    q.matrix.push_back(1,0,-2.);
    q.matrix.push_back(1,1,5.);
    q.matrix.finalize();
    q.rhs=cfdx::core::Vector(2,0.);
    q.rhs(0)=b0;q.rhs(1)=b1;
    q.diagonal={4.,5.};
    return q;
}
double coeff(const ScalarEquation& q,int i,int j) {
    for (auto p=q.matrix.row_offsets_data()[i];p<q.matrix.row_offsets_data()[i+1];++p)
        if (q.matrix.columns_data()[p]==static_cast<std::uint32_t>(j))
            return q.matrix.values_data()[p];
    return 0.;
}
}
int main() {
    run_case("v29_pseudo_time_adds_local_mass_and_preserves_stationary_root",[]{
        auto eq=two_cell();
        cfdx::core::Vector old(2,0.);old(0)=1.;old(1)=2.;
        auto a=sst_add_local_pseudo_time_v29(eq,old,0.5);
        EXPECT_TRUE(a.shifted_rows==2 && a.skipped_wall_rows==0);
        EXPECT_NEAR(coeff(eq,0,0),12.,1e-13);
        EXPECT_NEAR(coeff(eq,1,1),15.,1e-13);
        EXPECT_NEAR(coeff(eq,0,1),-1.,1e-13);
        EXPECT_NEAR(eq.rhs(0),10.,1e-13);
        EXPECT_NEAR(eq.rhs(1),28.,1e-13);
        // F(x*)=0 implies (A+M)x*=b+Mx*: the fixed point is unchanged.
        EXPECT_NEAR(coeff(eq,0,0)*old(0)+coeff(eq,0,1)*old(1)-eq.rhs(0),0.,1e-13);
        EXPECT_NEAR(coeff(eq,1,0)*old(0)+coeff(eq,1,1)*old(1)-eq.rhs(1),0.,1e-13);
    });
    run_case("v29_wall_owner_row_is_not_given_pseudo_mass",[]{
        auto eq=two_cell();
        cfdx::core::Vector old(2,1.);
        auto a=sst_add_local_pseudo_time_v29(eq,old,1.,{0,1});
        EXPECT_TRUE(a.shifted_rows==1 && a.skipped_wall_rows==1);
        EXPECT_NEAR(coeff(eq,0,0),8.,1e-13);
        EXPECT_NEAR(coeff(eq,1,1),5.,1e-13);
        EXPECT_NEAR(eq.rhs(1),8.,1e-13);
    });
    run_case("v29_invalid_pseudo_mass_is_transactional",[]{
        auto eq=two_cell();
        const auto original=eq;
        cfdx::core::Vector old(2,1.);
        old(1)=std::numeric_limits<double>::infinity();
        EXPECT_THROW(sst_add_local_pseudo_time_v29(eq,old,0.5),std::invalid_argument);
        EXPECT_NEAR(coeff(eq,0,0),coeff(original,0,0),0.);
        EXPECT_NEAR(eq.rhs(0),original.rhs(0),0.);
        EXPECT_THROW(sst_add_local_pseudo_time_v29(eq,old,0.),std::invalid_argument);
    });
    run_case("v29_cfl_changes_only_after_physical_selection",[]{
        SSTV29PseudoTimeState a;
        EXPECT_NEAR(a.cfl,0.5,0.);
        a.selected(1.,0.8);
        EXPECT_NEAR(a.cfl,1.,1e-13);
        a.selected(1.,0.999);
        EXPECT_NEAR(a.cfl,1.25,1e-13);
        a.rejected();
        EXPECT_NEAR(a.cfl,0.625,1e-13);
        EXPECT_THROW(a.selected(1.,1.1),std::invalid_argument);
        EXPECT_NEAR(a.cfl,0.625,1e-13);
    });
    run_case("v29_openfoam_full_owner_eliminates_offcell_column",[]{
        auto q=two_cell();
        const auto a=sst_apply_weighted_owner_matrix_v29(
            q,{1,0},{3.,0.},{1.,0.});
        EXPECT_TRUE(a.full_rows==1 && a.eliminated_offdiagonals==1);
        EXPECT_NEAR(coeff(q,0,0),4.,0.);
        EXPECT_NEAR(coeff(q,0,1),0.,0.);
        EXPECT_NEAR(q.rhs(0),12.,1e-13);
        EXPECT_NEAR(coeff(q,1,0),0.,0.);
        EXPECT_NEAR(q.rhs(1),14.,1e-13);
        EXPECT_NEAR(coeff(q,1,1),5.,0.);
        // Original interior residual at x=[3, 14/5] is unchanged.
        EXPECT_NEAR(-2.*3.+5.*(14./5.)-8.,0.,1e-13);
    });
    run_case("v29_openfoam_fractional_weight_matches_fractional_setvalue",[]{
        EXPECT_NEAR(sst_v29_openfoam_wall_weight(0.1),0.,0.);
        EXPECT_NEAR(sst_v29_openfoam_wall_weight(.55),.5,1e-14);
        EXPECT_NEAR(sst_v29_openfoam_wall_weight(1.),1.,1e-14);
        auto q=two_cell();
        EXPECT_THROW(sst_apply_weighted_owner_matrix_v29(
            q,{1,0},{10.,0.},{.55,0.},true),std::invalid_argument);
        EXPECT_NEAR(coeff(q,0,0),4.,0.);
        auto a=sst_apply_weighted_owner_matrix_v29(
            q,{1,0},{10.,0.},{.55,0.},false);
        EXPECT_TRUE(a.partial_rows==1 && a.full_rows==0);
        EXPECT_NEAR(coeff(q,0,0),8.,1e-13);
        EXPECT_NEAR(q.rhs(0),42.,1e-12);
        EXPECT_NEAR(coeff(q,1,0),-2.,0.);
    });
    run_case("v29_weighted_matrix_rejects_bad_input_before_mutation",[]{
        auto q=two_cell();
        EXPECT_THROW(sst_apply_weighted_owner_matrix_v29(
            q,{1,0},{1.,0.},{2.,0.}),std::invalid_argument);
        EXPECT_NEAR(coeff(q,0,0),4.,0.);
        EXPECT_NEAR(q.rhs(0),2.,0.);
    });
    run_case("v29_options_off_by_default_and_independent",[]{
        TurbulenceTransportControls c;c.model=TurbulenceModel::SST;
        EXPECT_TRUE(!c.sst_pseudo_transient_v29_experimental);
        EXPECT_TRUE(!c.sst_weighted_wall_matrix_v29_experimental);
        c.sst_pseudo_transient_v29_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_jfnk_v22_experimental=true;
        validate_turbulence_controls(c);
        c.sst_weighted_wall_matrix_v29_experimental=true;
        EXPECT_THROW(validate_turbulence_controls(c),std::invalid_argument);
        c.sst_owner_consistent_projection_v21=true;
        validate_turbulence_controls(c);
    });
    return run_all();
}
