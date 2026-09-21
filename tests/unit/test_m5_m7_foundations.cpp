#include "cfdx/physics/vof.h"
#include "cfdx/runtime/dynamic_mesh.h"
#include "cfdx/physics/fsi.h"
#include "common/test_harness.h"
#include <cmath>\n#include <limits>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::runtime;
using namespace cfdx::testing;

static Mesh make_cube()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for (std::size_t i=0;i<8;++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1}); m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4}); m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3}); m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for (std::size_t f=0;f<6;++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    return m;
}

int main()
{
    run_case("vof_advection_is_bounded_and_conservative", [] {
        Field<double,Location::CELL> a(8,"alpha","1",1);
        for (std::size_t i=0;i<a.size();++i) a(i)=i<4?1.0:0.0;
        const auto next = advect_volume_fraction_1d(a,0.5,1.0,1.0);
        double old_sum=0.0,new_sum=0.0;
        for (std::size_t i=0;i<a.size();++i) {
            old_sum += a(i); new_sum += next(i);
            EXPECT_TRUE(next(i)>=0.0 && next(i)<=1.0);
        }
        EXPECT_NEAR(new_sum,old_sum,1e-12);
    });

    run_case("vof_contact_angle_and_surface_force", [] {
        EXPECT_NEAR(contact_angle_wall_normal_component(60.0),0.5,1e-12);
        EXPECT_NEAR(continuum_surface_force(0.07,2.0,3.0),0.42,1e-12);
    });

    run_case("dynamic_mesh_translation_preserves_volume", [] {
        auto mesh = make_cube();
        std::vector<Vec3> d(mesh.n_points(), Vec3{0.1,-0.2,0.3});
        const auto result = apply_point_displacement(mesh,d);
        EXPECT_TRUE(result.valid);
        EXPECT_NEAR(result.min_volume,1.0,1e-12);
        EXPECT_NEAR(result.min_volume_ratio,1.0,1e-12);
    });

    run_case("dynamic_mesh_rejects_invalid_displacement", [] {
        auto mesh = make_cube();
        std::vector<Vec3> d(mesh.n_points(), Vec3{0,0,0});
        d[0].x = std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(apply_point_displacement(mesh,d),std::invalid_argument);
    });

    run_case("fsi_work_and_aitken_are_finite", [] {
        const std::vector<double> force{2.0,3.0};
        const std::vector<double> du{0.5,1.0};
        EXPECT_NEAR(interface_work(force,du),4.0,1e-12);
        const std::vector<double> r0{1.0,0.5};
        const std::vector<double> r1{0.5,0.25};
        const double omega=aitken_relaxation(0.5,r1,r0,0.05,1.0);
        EXPECT_TRUE(std::isfinite(omega));
        EXPECT_TRUE(omega>=0.05 && omega<=1.0);
    });

    return run_all();
}
