#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/numerics/divergence.h"
#include "cfdx/core/numerics/flux.h"
#include "cfdx/core/numerics/integrate.h"
#include "common/test_harness.h"

#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

static Mesh make_two_cell_cartesian()
{
    Mesh m;
    m.points().resize(12);
    const double p[12][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i,p[i][0],p[i][1],p[i][2]);

    const std::vector<std::vector<VertexIndex>> faces = {
        {0,3,2,1}, {4,5,6,7}, {0,1,5,4}, {3,7,6,2}, {0,4,7,3},
        {1,2,6,5}, {8,9,11,10}, {5,10,11,6},
        {1,8,10,5}, {2,6,11,9}, {1,2,9,8}
    };
    for (const auto& f : faces) m.faces().push_face(f);

    m.ownership().resize(faces.size());
    for (std::size_t f = 0; f < faces.size(); ++f) {
        const bool shared = f == 5;
        m.ownership().set_owner(f, shared ? 0 : (f < 6 ? 0 : 1));
        m.ownership().set_neighbour(f, shared ? 1 : FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    m.cells().push_cell({5,6,7,8,9,10});
    return m;
}

int main()
{
    run_case("gauss_gradient_respects_internal_face_orientation", [] {
        const Mesh m = make_two_cell_cartesian();
        Field<double,Location::CELL> phi(2,"phi","m",1);
        phi(0) = 1.0;
        phi(1) = 1.0;
        const auto grad = compute_gradient_gauss(phi,m);
        EXPECT_NEAR(grad.component_data(0)[0],0.0,1e-12);
        EXPECT_NEAR(grad.component_data(0)[1],0.0,1e-12);
        EXPECT_NEAR(grad.component_data(1)[0],0.0,1e-12);
        EXPECT_NEAR(grad.component_data(2)[1],0.0,1e-12);
    });

    run_case("divergence_uses_oriented_cell_volumes", [] {
        const Mesh m = make_two_cell_cartesian();
        Field<double,Location::FACE> flux(m.n_faces(),"phi","m3/s",1);
        for (std::size_t f = 0; f < m.n_faces(); ++f) {
            const auto off = m.faces().offsets_data()[f];
            const auto n = m.faces().offsets_data()[f+1] - off;
            const auto fg = compute_face_geometry(
                m.points().x_data(),m.points().y_data(),m.points().z_data(),
                m.faces().vertices_data(),off,n);
            flux(f) = fg.Sf.x;
        }
        const auto div = compute_divergence(flux,m);
        EXPECT_NEAR(div(0),0.0,1e-12);
        EXPECT_NEAR(div(1),0.0,1e-12);
    });

    run_case("flux_divergence_conserves_internal_face_flux", [] {
        const Mesh m = make_two_cell_cartesian();
        Field<double,Location::FACE> U(m.n_faces(),"U","m/s",3);
        U.fill(0.0);
        for (std::size_t f = 0; f < m.n_faces(); ++f)
            U(f,0) = 1.0;

        const auto phi = compute_flux(U,m);
        const auto div = compute_divergence(phi,m);

        // The only internal face is shared by both cells. Its contribution
        // must enter one cell with the opposite sign in the neighbour.
        EXPECT_NEAR(div(0) * 1.0 + div(1) * 1.0, 0.0, 1e-12);
        EXPECT_NEAR(phi(5), 1.0, 1e-12);
    });

    run_case("volume_integral_uses_oriented_cell_volumes", [] {
        const Mesh m = make_two_cell_cartesian();
        Field<double,Location::CELL> phi(2,"phi","1",1);
        phi.fill(2.0);
        EXPECT_NEAR(volume_integrate(phi,m),4.0,1e-12);
    });

    return run_all();
}
