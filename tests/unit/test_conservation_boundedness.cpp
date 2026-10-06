#include "cfdx/core/numerics/conservation.h"
#include "cfdx/core/mesh/mesh.h"
#include "common/test_harness.h"
#include <limits>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

static Mesh make_unit_cube()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i=0; i<8; ++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for (std::size_t f=0; f<6; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    return m;
}

static Mesh make_two_cell_chain()
{
    Mesh m;
    m.points().resize(12);
    const double p[12][3] = {
        {0,0,0},{0.5,0,0},{0.5,1,0},{0,1,0},{0,0,1},{0.5,0,1},
        {1,0,0},{1,1,0},{1,0,1},{1,1,1},{0.5,1,1},{0.5,0,1}
    };
    for (std::size_t i=0; i<12; ++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);

    m.faces().push_face({0,4,5,3});
    m.faces().push_face({1,2,10,11});
    m.faces().push_face({0,1,11,4});
    m.faces().push_face({3,5,10,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,11,10,5});
    m.faces().push_face({6,7,9,8});
    m.faces().push_face({1,6,8,11});
    m.faces().push_face({2,10,9,7});
    m.faces().push_face({1,2,7,6});
    m.faces().push_face({11,8,9,10});

    m.ownership().resize(11);
    for (std::size_t f=0; f<11; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1,0);
    m.ownership().set_neighbour(1,1);
    for (std::size_t f=6; f<11; ++f) m.ownership().set_owner(f,1);

    m.cells().push_cell({0,1,2,3,4,5});
    m.cells().push_cell({1,6,7,8,9,10});
    return m;
}

int main()
{
    run_case("uniform_closed_flux_has_zero_local_and_global_balance", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::FACE> phi(m.n_faces(),"phi","m3/s",1);
        phi.fill(0.0);
        const auto d = audit_face_flux_conservation(m, phi);
        EXPECT_NEAR(d.global_boundary_flux, 0.0, 1e-15);
        EXPECT_NEAR(d.max_local_imbalance, 0.0, 1e-15);
        EXPECT_NEAR(d.l2_local_imbalance, 0.0, 1e-15);
        EXPECT_NEAR(d.global_cell_balance, 0.0, 1e-15);
        EXPECT_TRUE(d.worst_cell == 0);
        EXPECT_TRUE(d.boundary_faces == 6);
    });

    run_case("internal_face_is_accounted_once", [] {
        const Mesh m = make_two_cell_chain();
        Field<double,Location::FACE> phi(m.n_faces(),"phi","m3/s",1);
        phi.fill(0.0);
        phi(1) = 2.5;
        const auto d = audit_face_flux_conservation(m, phi);
        EXPECT_NEAR(d.l1_local_imbalance, 5.0, 1e-15);
        EXPECT_NEAR(d.global_boundary_flux, 0.0, 1e-15);
        EXPECT_NEAR(d.global_cell_balance, 0.0, 1e-15);
        EXPECT_TRUE(d.worst_cell_imbalance > 0.0);
    });

    run_case("nonzero_boundary_flux_matches_global_cell_balance", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::FACE> phi(m.n_faces(),"phi","m3/s",1);
        phi.fill(0.0);
        phi(0) = 3.0;
        const auto d = audit_face_flux_conservation(m, phi);
        EXPECT_NEAR(d.global_boundary_flux, 3.0, 1e-15);
        EXPECT_NEAR(d.global_cell_balance, 3.0, 1e-15);
    });

    run_case("boundedness_detects_overshoot_without_tolerance_inflation", [] {
        const auto d = audit_boundedness({0.0, 0.25, 0.75, 1.0}, 0.0, 1.0);
        EXPECT_TRUE(d.bounded());
        EXPECT_NEAR(d.minimum, 0.0, 1e-15);
        EXPECT_NEAR(d.maximum, 1.0, 1e-15);

        const auto b = audit_boundedness({0.0, 0.25, 1.0001, -0.01}, 0.0, 1.0);
        EXPECT_TRUE(!b.bounded());
        EXPECT_TRUE(b.above_upper == 1);
        EXPECT_TRUE(b.below_lower == 1);
        EXPECT_TRUE(b.maximum_violation > 0.0);
        EXPECT_TRUE(b.minimum_violation > 0.0);
    });

    run_case("cell_field_boundedness_uses_native_field_storage", [] {
        Field<double,Location::CELL> k(3,"k","m2/s2",1);
        k(0)=1e-12; k(1)=0.2; k(2)=2.0;
        EXPECT_TRUE(audit_positive_field(k).bounded());
        k(1)=-1e-8;
        EXPECT_TRUE(!audit_positive_field(k).bounded());
    });

    run_case("positive_fields_reject_negative_and_nonfinite_values", [] {
        EXPECT_TRUE(audit_positive_field({1e-12, 0.1, 10.0}).bounded());
        const auto bad = audit_positive_field({0.1, -1e-8, 1.0});
        EXPECT_TRUE(!bad.bounded());
        EXPECT_TRUE(bad.below_lower == 1);
        const auto nan = audit_positive_field({1.0, std::numeric_limits<double>::quiet_NaN()});
        EXPECT_TRUE(!nan.bounded());
        EXPECT_TRUE(nan.nonfinite == 1);
    });

    run_case("generic_cell_balance_is_reusable_for_any_conserved_quantity", [] {
        const std::vector<double> face_balance{2.0, -1.0, -1.0};
        const std::vector<double> source_integral{-2.0, 0.5, 0.5};
        // Each cell must close independently: 2-2=0, -1+0.5-(-0.5)=0,
        // and -1+0.5-(-0.5)=0.
        const std::vector<double> accumulation{0.0, -0.5, -0.5};
        const auto r = audit_cell_balance(face_balance, source_integral, &accumulation);
        EXPECT_TRUE(r.finite());
        EXPECT_TRUE(r.closed(1e-15));
        EXPECT_NEAR(r.residual, 0.0, 1e-15);
        EXPECT_NEAR(r.l1_cell_residual, 0.0, 1e-15);
        EXPECT_NEAR(r.l2_cell_residual, 0.0, 1e-15);
        EXPECT_NEAR(r.max_cell_residual, 0.0, 1e-15);
    });

    run_case("generic_cell_balance_detects_nonfinite_component", [] {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const std::vector<double> face_balance{1.0, nan};
        const std::vector<double> source_integral{-1.0, 0.0};
        const auto r = audit_cell_balance(face_balance, source_integral);
        EXPECT_TRUE(!r.finite());
        EXPECT_TRUE(r.nonfinite_balance == 1);
    });

    run_case("generic_mesh_balance_reuses_independent_face_reconstruction", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::FACE> flux(m.n_faces(),"flux","kg/s",1);
        Field<double,Location::CELL> source(m.n_cells(),"source","kg/s",1);
        flux.fill(0.0);
        source.fill(0.0);
        flux(0) = 2.0;
        source(0) = -2.0;
        const auto r = audit_cell_balance(m, flux, source);
        EXPECT_TRUE(r.closed(1e-15));
        EXPECT_NEAR(r.normalized_residual, 0.0, 1e-15);
    });

    run_case("integrated_balance_is_independent_of_linear_residual", [] {
        const Mesh m = make_unit_cube();
        Field<double,Location::FACE> flux(m.n_faces(),"flux","u",1);
        Field<double,Location::CELL> source(m.n_cells(),"source","u/V",1);
        flux.fill(0.0); source.fill(0.0);
        flux(0) = 2.0;
        source(0) = -2.0;
        const std::vector<double> volumes{1.0};
        const auto r = audit_integrated_balance(m, flux, source, volumes);
        EXPECT_NEAR(r.residual, 0.0, 1e-15);
        EXPECT_NEAR(r.normalized_residual, 0.0, 1e-15);
        EXPECT_TRUE(r.worst_cell == 0);
    });

    run_case("boundedness_reports_worst_cell_without_clipping", [] {
        Field<double,Location::CELL> phi(4,"phi","1",1);
        phi(0)=0.0; phi(1)=1.0; phi(2)=1.25; phi(3)=-0.2;
        const auto r = audit_boundedness(phi,0.0,1.0);
        EXPECT_TRUE(!r.bounded());
        EXPECT_TRUE(r.above_upper == 1);
        EXPECT_TRUE(r.below_lower == 1);
        EXPECT_TRUE(r.worst_cell == 2 || r.worst_cell == 3);
        EXPECT_TRUE(r.worst_violation > 0.0);
    });

    return run_all();
}
