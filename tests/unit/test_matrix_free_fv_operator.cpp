// M0.8/M0.7 — Assembled vs matrix-free finite-volume diffusion equivalence

#include "cfdx/core/linalg/matrix_free_fv_operator.h"
#include "cfdx/core/numerics/matrix_free.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/solvers/scalar_diffusion.h"
#include "common/test_harness.h"
#include <limits>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

static Mesh make_two_cell_mesh() {
    Mesh m;
    m.points().resize(12);
    const double p[12][3] = {
        {0,0,0},{1,0,0},{2,0,0},{0,1,0},{1,1,0},{2,1,0},
        {0,0,1},{1,0,1},{2,0,1},{0,1,1},{1,1,1},{2,1,1}};
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, p[i][0], p[i][1], p[i][2]);

    m.faces().push_face({0,6,9,3});
    m.faces().push_face({0,1,7,6});
    m.faces().push_face({3,9,10,4});
    m.faces().push_face({0,3,4,1});
    m.faces().push_face({6,7,10,9});
    m.faces().push_face({7,10,4,1});
    m.faces().push_face({2,5,11,8});
    m.faces().push_face({1,2,8,7});
    m.faces().push_face({4,10,11,5});
    m.faces().push_face({1,4,5,2});
    m.faces().push_face({7,8,11,10});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 5; ++f) {
        m.ownership().set_owner(f, 0);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(5, 0);
    m.ownership().set_neighbour(5, 1);
    for (std::size_t f = 6; f < 11; ++f) {
        m.ownership().set_owner(f, 1);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    m.cells().push_cell({5,6,7,8,9,10});
    return m;
}

int main() {
    run_case("matrix_free_matches_assembled_diffusion", []() {
        Mesh m = make_two_cell_mesh();

        std::vector<Vec3> centres(m.n_cells());
        std::vector<Vec3> Sf(m.n_faces());
        const GeometryCache geometry = make_geometry_cache(m);
        centres = geometry.cell_centres;
        Sf = geometry.face_Sf;

        constexpr double gamma = 2.5;
        FvDiffusionOperator mf(m, centres, Sf, gamma);
        MatrixFreeFvDiffusionOperator mf_public(m, geometry, gamma);
        EXPECT_TRUE(mf_public.rows() == m.n_cells());
        EXPECT_TRUE(mf_public.cols() == m.n_cells());

        DirichletBoundary boundary;
        boundary.face_values.assign(m.n_faces(), std::numeric_limits<double>::quiet_NaN());
        Vector rhs;
        SparseMatrix A = assemble_cell_diffusion_matrix(
            m, boundary, gamma, {0.0, 0.0}, rhs, geometry);

        Vector x(2);
        x(0) = 3.0;
        x(1) = -1.0;
        Vector y_mf;
        mf.apply(x, y_mf);
        const std::vector<double> y_a = A.matvec(x);

        EXPECT_NEAR(y_mf(0), y_a[0], 1e-12);
        EXPECT_NEAR(y_mf(1), y_a[1], 1e-12);
        EXPECT_NEAR(y_mf(0) + y_mf(1), 0.0, 1e-12);
    });

    return run_all();
}
