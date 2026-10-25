// A conductive CHT interface is conformal and impermeable, and both regions
// advance the same physical step. Invalid contracts fail before any mutation.
#include "audit/chain_fixture.h"
#include "cfdx/physics/cht_solver.h"
#include "common/test_harness.h"
#include <limits>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

struct PairFixture {
    Mesh a=make_1d_chain(1), b=make_1d_chain(1);
    FvGeometry ga,gb;
    Field<double,Location::FACE> fa{a.n_faces(),"F","kg/s",1};
    Field<double,Location::FACE> fb{b.n_faces(),"F","kg/s",1};
    Field<double,Location::CELL> ta{1,"T","K",1},tb{1,"T","K",1};
    Field<double,Location::CELL> qa{1,"Q","W/m3",1},qb{1,"Q","W/m3",1};
    EnergySolverControls ea,eb;
    ChtInterfaceControls c;
    PairFixture()
    {
        for (std::size_t i=0;i<b.n_points();++i)
            b.points().set(i,b.points().x(i)+1.0,b.points().y(i),b.points().z(i));
        ga=build_fv_geometry(a); gb=build_fv_geometry(b);
        fa.fill(0.0); fb.fill(0.0); qa.fill(0.0); qb.fill(0.0);
        ta.fill(300.0); tb.fill(300.0);
        ea.density=eb.density=1.0; ea.cp=eb.cp=1.0;
        ea.conductivity=eb.conductivity=1.0; ea.dt=eb.dt=1.0;
        ea.relaxation=eb.relaxation=1.0;
        c.region1_patch="outlet"; c.region2_patch="inlet";
        c.max_iterations=200;
    }
    ChtSolveResult solve()
    {
        return solve_two_region_cht(a,ga,b,gb,fa,fb,ta,tb,qa,qb,ea,eb,c);
    }
    void unchanged() const
    {
        EXPECT_NEAR(ta(0),300.0,0.0); EXPECT_NEAR(tb(0),300.0,0.0);
    }
};

int main()
{
    run_case("cht_accepts_adjacent_opposed_interface", [] {
        PairFixture f;
        const auto pairs=match_cht_interface(f.a,f.ga,f.b,f.gb,f.c);
        EXPECT_TRUE(pairs.size()==1);
        EXPECT_TRUE(f.solve().converged);
        f.unchanged();
    });
    run_case("cht_rejects_coincident_domains_with_parallel_outward_normals", [] {
        PairFixture f;
        f.c.region2_patch="outlet";
        EXPECT_THROW_WITH(match_cht_interface(f.a,f.ga,f.a,f.ga,f.c),
                          std::invalid_argument,"opposite outward normals");
    });
    run_case("cht_rejects_empty_interface", [] {
        PairFixture f;
        f.a.boundary().patch(1).face_ids.clear();
        f.b.boundary().patch(0).face_ids.clear();
        f.ga=build_fv_geometry(f.a); f.gb=build_fv_geometry(f.b);
        EXPECT_THROW_WITH(f.solve(),std::invalid_argument,"non-empty");
        f.unchanged();
    });
    run_case("cht_rejects_nonfinite_interface_controls", [] {
        PairFixture f;
        f.c.matching_tolerance=std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW_WITH(f.solve(),std::invalid_argument,"interface controls");
        f.unchanged();
    });
    run_case("cht_rejects_different_physical_steps", [] {
        for (const double dt : {0.0,0.5,2.0}) {
            PairFixture f;
            f.eb.dt=dt;
            EXPECT_THROW_WITH(f.solve(),std::invalid_argument,"same physical time step");
            f.unchanged();
        }
    });
    run_case("cht_rejects_mass_transport_at_impermeable_interface", [] {
        for (const bool first : {true,false}) {
            PairFixture f;
            if(first) f.fa(1)=1e-12;
            else f.fb(0)=-1e-12;
            EXPECT_THROW_WITH(f.solve(),std::invalid_argument,"impermeable");
            f.unchanged();
        }
    });
    run_case("cht_rejects_malformed_second_region_before_first_solve", [] {
        PairFixture f;
        f.qb(0)=std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW_WITH(f.solve(),std::invalid_argument,"region source");
        f.unchanged();
    });
    run_case("cht_rejects_nonboundary_interface_face", [] {
        PairFixture f;
        f.a.ownership().set_neighbour(1,0);
        EXPECT_THROW_WITH(f.solve(),std::invalid_argument,"boundary faces");
        f.unchanged();
    });
    return run_all();
}
