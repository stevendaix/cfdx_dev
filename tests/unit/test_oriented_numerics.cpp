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


static Mesh make_unit_tetrahedron()
{
    Mesh m;
    m.points().resize(4);
    m.points().set(0,0.0,0.0,0.0);
    m.points().set(1,1.0,0.0,0.0);
    m.points().set(2,0.0,1.0,0.0);
    m.points().set(3,0.0,0.0,1.0);

    // The geometry cache canonicalizes every face to owner -> exterior.
    const std::vector<std::vector<VertexIndex>> faces = {
        {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}
    };
    for (const auto& face : faces) m.faces().push_face(face);

    m.ownership().resize(faces.size());
    for (std::size_t f = 0; f < faces.size(); ++f) {
        m.ownership().set_owner(f, 0);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3});
    return m;
}


static Mesh make_skew_tetrahedron()
{
    Mesh m;
    m.points().resize(4);
    m.points().set(0,0.0,0.0,0.0);
    m.points().set(1,1.0,0.0,0.0);
    m.points().set(2,0.0,2.0,0.0);
    m.points().set(3,1.0,1.0,1.0);

    // Same owner->exterior winding pattern as the unit reference tetrahedron.
    // Edge dot products (1,0,0).(1,1,1) = 1 and (0,2,0).(1,1,1) = 2 are
    // nonzero, so this cell is genuinely skewed, unlike the reference cell.
    const std::vector<std::vector<VertexIndex>> faces = {
        {1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}
    };
    for (const auto& face : faces) m.faces().push_face(face);

    m.ownership().resize(faces.size());
    for (std::size_t f = 0; f < faces.size(); ++f) {
        m.ownership().set_owner(f, 0);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3});
    return m;
}

int main()
{
    run_case("reference_tetrahedron_geometry_cache_satisfies_gauss_moment_identity", [] {
        const Mesh m = make_unit_tetrahedron();
        const GeometryCache geometry = make_geometry_cache(m);
        EXPECT_NEAR(geometry.cell_volumes[0], 1.0 / 6.0, 1e-12);

        const auto component = [](const Vec3& value, std::size_t d) {
            return d == 0 ? value.x : (d == 1 ? value.y : value.z);
        };

        // Check the same oriented area closure used by the Lean reference cell.
        Vec3 area_sum;
        for (std::size_t f = 0; f < m.n_faces(); ++f)
            area_sum = area_sum + geometry.face_Sf[f];
        EXPECT_NEAR(area_sum.x, 0.0, 1e-12);
        EXPECT_NEAR(area_sum.y, 0.0, 1e-12);
        EXPECT_NEAR(area_sum.z, 0.0, 1e-12);

        // The first-moment tensor sum_f C_f outer Sf_f must equal V * I.
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                double moment = 0.0;
                for (std::size_t f = 0; f < m.n_faces(); ++f) {
                    moment += component(geometry.face_centres[f], i) *
                              component(geometry.face_Sf[f], j);
                }
                const double expected = (i == j) ? geometry.cell_volumes[0] : 0.0;
                EXPECT_NEAR(moment, expected, 1e-12);
            }
        }
    });

    run_case("skew_tetrahedron_geometry_cache_satisfies_gauss_moment_identity", [] {
        const Mesh m = make_skew_tetrahedron();
        const GeometryCache geometry = make_geometry_cache(m);

        // Edge determinant |det(p1-p0, p2-p0, p3-p0)| = 2, hence
        // V = 2/6 = 1/3, matching the Lean skew_tetra_volume_edge_determinant.
        EXPECT_NEAR(geometry.cell_volumes[0], 1.0 / 3.0, 1e-12);
        EXPECT_TRUE(is_valid(geometry, m));

        const auto component = [](const Vec3& value, std::size_t d) {
            return d == 0 ? value.x : (d == 1 ? value.y : value.z);
        };

        // Outward-oriented face area vectors, matching the Lean fixture
        // skewTetraOrientedArea in CFDX.Geometry.SkewTetrahedron.
        const double expectedSf[4][3] = {
            {1.0, 0.5, -0.5}, {-1.0, 0.0, 1.0},
            {0.0, -0.5, 0.5}, {0.0, 0.0, -1.0}};
        for (std::size_t f = 0; f < 4; ++f) {
            EXPECT_NEAR(component(geometry.face_Sf[f], 0), expectedSf[f][0], 1e-12);
            EXPECT_NEAR(component(geometry.face_Sf[f], 1), expectedSf[f][1], 1e-12);
            EXPECT_NEAR(component(geometry.face_Sf[f], 2), expectedSf[f][2], 1e-12);
        }

        // Exact oriented surface closure on the skewed cell.
        Vec3 area_sum;
        for (std::size_t f = 0; f < m.n_faces(); ++f)
            area_sum = area_sum + geometry.face_Sf[f];
        EXPECT_NEAR(area_sum.x, 0.0, 1e-12);
        EXPECT_NEAR(area_sum.y, 0.0, 1e-12);
        EXPECT_NEAR(area_sum.z, 0.0, 1e-12);

        // The first-moment tensor sum_f C_f outer Sf_f equals V * I on the
        // skewed cell: all nine components, matching the Lean theorem
        // skew_tetra_first_moment. This extends the moment conformance from
        // the axis-aligned reference simplex to an oblique cell.
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                double moment = 0.0;
                for (std::size_t f = 0; f < m.n_faces(); ++f) {
                    moment += component(geometry.face_centres[f], i) *
                              component(geometry.face_Sf[f], j);
                }
                const double expected = (i == j) ? geometry.cell_volumes[0] : 0.0;
                EXPECT_NEAR(moment, expected, 1e-12);
            }
        }
    });

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

    run_case("gauss_gradient_constant_field_zero_for_every_cell_component", [] {
        const Mesh m = make_two_cell_cartesian();
        const GeometryCache geometry = make_geometry_cache(m);
        const auto& ownership = m.ownership();
        const auto* cell_faces = m.cells().faces_data();
        const auto* cell_offsets = m.cells().offsets_data();

        // Verify the geometric closure hypothesis used by the Lean theorem
        // on the actual geometry cache and cell-face CSR for this mesh.
        for (std::size_t c = 0; c < m.n_cells(); ++c) {
            Vec3 area_sum;
            for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
                const std::size_t f = cell_faces[k];
                const Vec3 oriented_area = (ownership.owner(f) == c)
                    ? geometry.face_Sf[f]
                    : geometry.face_Sf[f] * (-1.0);
                area_sum = area_sum + oriented_area;
            }
            EXPECT_NEAR(area_sum.x, 0.0, 1e-12);
            EXPECT_NEAR(area_sum.y, 0.0, 1e-12);
            EXPECT_NEAR(area_sum.z, 0.0, 1e-12);
        }

        Field<double,Location::CELL> phi(m.n_cells(), "phi", "1", 1);
        phi.fill(3.25);
        const auto grad = compute_gradient_gauss(phi, m, geometry);
        for (std::size_t c = 0; c < m.n_cells(); ++c) {
            for (std::size_t d = 0; d < 3; ++d) {
                EXPECT_NEAR(grad(c, d), 0.0, 1e-12);
            }
        }
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
        Field<double,Location::FACE> phi(m.n_faces(),"phi","m3/s",1);
        phi.fill(0.0);
        phi(5) = 1.0;

        const auto div = compute_divergence(phi,m);

        // Only the shared internal face carries flux. Owner and neighbour
        // must receive equal and opposite contributions.
        EXPECT_NEAR(div(0), 1.0, 1e-12);
        EXPECT_NEAR(div(1), -1.0, 1e-12);
        EXPECT_NEAR(div(0) + div(1), 0.0, 1e-12);
    });

    run_case("volume_integral_uses_oriented_cell_volumes", [] {
        const Mesh m = make_two_cell_cartesian();
        Field<double,Location::CELL> phi(2,"phi","1",1);
        phi.fill(2.0);
        EXPECT_NEAR(volume_integrate(phi,m),4.0,1e-12);
    });

    return run_all();
}
