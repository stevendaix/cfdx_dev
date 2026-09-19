// M0.2-T02 — Tests for CellGeometry
// Validation Level 1 (invariants) + Level 2 (analytical cube).

#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    // Cube unité [0,1]^3. 6 faces, 8 sommets.
    // Sommets :
    //   0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0)
    //   4:(0,0,1) 5:(1,0,1) 6:(1,1,1) 7:(0,1,1)
    // Faces orientées vers l'extérieur (CCW vue de l'extérieur) :
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
        std::uint32_t fv[6][4] = {
            {0,3,2,1}, {4,5,6,7}, {0,1,5,4},
            {3,7,6,2}, {0,4,7,3}, {1,2,6,5}
        };
        Vec3 fc[6], fs[6];
        for (int i = 0; i < 6; ++i) {
            auto g = compute_face_geometry(px, py, pz, fv[i], 0, 4);
            fc[i] = g.centre; fs[i] = g.Sf;
        }
        std::uint32_t face_ids[6] = {0,1,2,3,4,5};
        auto cell = compute_cell_geometry(fc, fs, face_ids, 6);
        EXPECT_NEAR(cell.volume, 1.0, 1e-12);
    });

    run_case("cube_centre", []() {
        std::vector<double> px = {0,1,1,0, 0,1,1,0};
        std::vector<double> py = {0,0,1,1, 0,0,1,1};
        std::vector<double> pz = {0,0,0,0, 1,1,1,1};
        std::uint32_t fv[6][4] = {
            {0,1,2,3}, {4,7,6,5}, {0,4,5,1},
            {3,2,6,7}, {0,3,7,4}, {1,5,6,2}
        };
        Vec3 fc[6], fs[6];
        for (int i = 0; i < 6; ++i) {
            auto g = compute_face_geometry(px, py, pz, fv[i], 0, 4);
            fc[i] = g.centre; fs[i] = g.Sf;
        }
        std::uint32_t face_ids[6] = {0,1,2,3,4,5};
        auto cell = compute_cell_geometry(fc, fs, face_ids, 6);
        EXPECT_NEAR(cell.centre.x, 0.5, 1e-12);
        EXPECT_NEAR(cell.centre.y, 0.5, 1e-12);
        EXPECT_NEAR(cell.centre.z, 0.5, 1e-12);
    });

    run_case("cube_closure_invariant", []() {
        // Pour une cellule fermée : Σ Sf = 0 (§18)
        std::vector<double> px = {0,1,1,0, 0,1,1,0};
        std::vector<double> py = {0,0,1,1, 0,0,1,1};
        std::vector<double> pz = {0,0,0,0, 1,1,1,1};
        std::uint32_t fv[6][4] = {
            {0,1,2,3}, {4,7,6,5}, {0,4,5,1},
            {3,2,6,7}, {0,3,7,4}, {1,5,6,2}
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
        // Tétraèdre : (0,0,0), (1,0,0), (0,1,0), (0,0,1) → volume = 1/6
        std::vector<double> px = {0,1,0,0};
        std::vector<double> py = {0,0,1,0};
        std::vector<double> pz = {0,0,0,1};
        // Faces orientées vers l'extérieur :
        //   base (z=0) : 0 2 1  → normale -Z
        //   face xz   : 0 1 3  → normale -Y
        //   face yz   : 0 3 2  → normale -X
        //   face oblique : 1 2 3 → normale +X+Y+Z (extérieure)
        std::uint32_t verts[4][3] = {
            {0,2,1}, {0,1,3}, {0,3,2}, {1,2,3}
        };
        Vec3 fc[4], fs[4];
        for (int i = 0; i < 4; ++i) {
            auto g = compute_face_geometry(px, py, pz, verts[i], 0, 3);
            fc[i] = g.centre; fs[i] = g.Sf;
        }
        std::uint32_t face_ids[4] = {0,1,2,3};
        auto cell = compute_cell_geometry(fc, fs, face_ids, 4);
        EXPECT_NEAR(cell.volume, 1.0 / 6.0, 1e-12);
    });

    run_case("empty_cell_throws", []() {
        Vec3 fc[1]; Vec3 fs[1];
        std::uint32_t faces[1] = {0};
        EXPECT_THROW(compute_cell_geometry(fc, fs, faces, 0), std::runtime_error);
    });

    return run_all();
}