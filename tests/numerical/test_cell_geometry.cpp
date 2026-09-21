// M0.2-T02 — Tests for CellGeometry
// Validation Level 1 (invariants) + Level 2 (analytical cube).

#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/mesh/index_types.h"
#include "cfdx/core/mesh/mesh.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    // Cube unité [0,1]^3. Faces orientées vers l'extérieur (CCW vue de l'extérieur).
    //   bas   (z=0) : 0 3 2 1   → Sf = (0,0,-1)
    //   haut  (z=1) : 4 5 6 7   → Sf = (0,0,+1)
    //   avant (y=0) : 0 1 5 4   → Sf = (0,-1,0)
    //   arrière(y=1) : 3 7 6 2   → Sf = (0,+1,0)
    //   gauche (x=0): 0 4 7 3   → Sf = (-1,0,0)
    //   droite(x=1) : 1 2 6 5   → Sf = (+1,0,0)

    run_case("cube_volume", []() {
        std::vector<double> px = {0,1,1,0, 0,1,1,0};
        std::vector<double> py = {0,0,1,1, 0,0,1,1};
        std::vector<double> pz = {0,0,0,0, 1,1,1,1};
        VertexIndex fv[6][4] = {
            {0,3,2,1}, {4,5,6,7}, {0,1,5,4},
            {3,7,6,2}, {0,4,7,3}, {1,2,6,5}
        };
        Vec3 fc[6], fs[6];
        for (int i = 0; i < 6; ++i) {
            auto g = compute_face_geometry(px, py, pz, fv[i], 0, 4);
            fc[i] = g.centre; fs[i] = g.Sf;
        }
        FaceIndex face_ids[6] = {0,1,2,3,4,5};
        auto cell = compute_cell_geometry(fc, fs, face_ids, 6);
        EXPECT_NEAR(cell.volume, 1.0, 1e-12);
    });

    run_case("cube_centre", []() {
        std::vector<double> px = {0,1,1,0, 0,1,1,0};
        std::vector<double> py = {0,0,1,1, 0,0,1,1};
        std::vector<double> pz = {0,0,0,0, 1,1,1,1};
        VertexIndex fv[6][4] = {
            {0,3,2,1}, {4,5,6,7}, {0,1,5,4},
            {3,7,6,2}, {0,4,7,3}, {1,2,6,5}
        };
        Vec3 fc[6], fs[6];
        for (int i = 0; i < 6; ++i) {
            auto g = compute_face_geometry(px, py, pz, fv[i], 0, 4);
            fc[i] = g.centre; fs[i] = g.Sf;
        }
        FaceIndex face_ids[6] = {0,1,2,3,4,5};
        auto cell = compute_cell_geometry(fc, fs, face_ids, 6);
        EXPECT_NEAR(cell.centre.x, 0.5, 1e-12);
        EXPECT_NEAR(cell.centre.y, 0.5, 1e-12);
        EXPECT_NEAR(cell.centre.z, 0.5, 1e-12);
    });

    run_case("cube_closure_invariant", []() {
        std::vector<double> px = {0,1,1,0, 0,1,1,0};
        std::vector<double> py = {0,0,1,1, 0,0,1,1};
        std::vector<double> pz = {0,0,0,0, 1,1,1,1};
        VertexIndex fv[6][4] = {
            {0,3,2,1}, {4,5,6,7}, {0,1,5,4},
            {3,7,6,2}, {0,4,7,3}, {1,2,6,5}
        };
        Vec3 fs[6];
        for (int i = 0; i < 6; ++i) {
            auto g = compute_face_geometry(px, py, pz, fv[i], 0, 4);
            fs[i] = g.Sf;
        }
        Vec3 sum;
        for (int i = 0; i < 6; ++i) sum = sum + fs[i];
        EXPECT_NEAR(sum.x, 0.0, 1e-12);
        EXPECT_NEAR(sum.y, 0.0, 1e-12);
        EXPECT_NEAR(sum.z, 0.0, 1e-12);
    });

    run_case("tetrahedron_volume", []() {
        std::vector<double> px = {0,1,0,0};
        std::vector<double> py = {0,0,1,0};
        std::vector<double> pz = {0,0,0,1};
        VertexIndex verts[4][3] = {
            {0,2,1}, {0,1,3}, {0,3,2}, {1,2,3}
        };
        Vec3 fc[4], fs[4];
        for (int i = 0; i < 4; ++i) {
            auto g = compute_face_geometry(px, py, pz, verts[i], 0, 3);
            fc[i] = g.centre; fs[i] = g.Sf;
        }
        FaceIndex face_ids[4] = {0,1,2,3};
        auto cell = compute_cell_geometry(fc, fs, face_ids, 4);
        EXPECT_NEAR(cell.volume, 1.0 / 6.0, 1e-12);
    });

    run_case("oriented_two_cell_shared_face_has_positive_volumes", []() {
        Mesh m;
        m.points().resize(12);
        const double p[12][3] = {
            {0,0,0},{1,0,0},{1,1,0},{0,1,0},
            {0,0,1},{1,0,1},{1,1,1},{0,1,1},
            {2,0,0},{2,1,0},{2,0,1},{2,1,1}
        };
        for (std::size_t i=0; i<12; ++i)
            m.points().set(i,p[i][0],p[i][1],p[i][2]);

        const std::vector<std::vector<VertexIndex>> faces = {
            {0,3,2,1}, {4,5,6,7}, {0,1,5,4}, {3,7,6,2}, {0,4,7,3},
            {1,2,6,5}, {8,9,11,10}, {5,10,11,6},
            {1,8,10,5}, {2,6,11,9}, {4,7,9,8}
        };
        for (const auto& face : faces) m.faces().push_face(face);
        m.ownership().resize(faces.size());
        for (std::size_t f=0; f<faces.size(); ++f) {
            m.ownership().set_owner(f, f == 5 ? 0 : (f < 6 ? 0 : 1));
            m.ownership().set_neighbour(f, f == 5 ? 1 : FaceOwnership::BOUNDARY);
        }
        m.cells().push_cell({0,1,2,3,4,5});
        m.cells().push_cell({5,6,7,8,9,10});

        Vec3 fc[11], fs[11];
        for (std::size_t f=0; f<11; ++f) {
            const auto off=m.faces().offsets_data()[f];
            const auto n=m.faces().offsets_data()[f+1]-off;
            const auto g=compute_face_geometry(
                m.points().x_data(),m.points().y_data(),m.points().z_data(),
                m.faces().vertices_data(),off,n);
            fc[f]=g.centre; fs[f]=g.Sf;
        }

        const auto g0=compute_cell_geometry_oriented(
            fc,fs,m.cells().faces_data(),6,0,m.ownership());
        const auto g1=compute_cell_geometry_oriented(
            fc,fs,m.cells().faces_data()+6,6,1,m.ownership());
        EXPECT_NEAR(g0.volume,1.0,1e-12);
        EXPECT_NEAR(g1.volume,1.0,1e-12);
        EXPECT_NEAR(g0.centre.x,0.5,1e-12);
        EXPECT_NEAR(g1.centre.x,1.5,1e-12);
    });

    run_case("empty_cell_throws", []() {
        Vec3 fc[1]; Vec3 fs[1];
        FaceIndex faces[1] = {0};
        EXPECT_THROW(compute_cell_geometry(fc, fs, faces, 0), std::runtime_error);
    });

    return run_all();
}
