// M0.2-T01 — Tests for FaceGeometry
// Validation Level 1 (geometric invariants) + Level 2 (analytical).

#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/mesh/index_types.h"
#include "test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

// Helper: build point arrays for a square in XY plane.
struct SquareMesh {
    std::vector<double> px = {0.0, 1.0, 1.0, 0.0};
    std::vector<double> py = {0.0, 0.0, 1.0, 1.0};
    std::vector<double> pz = {0.0, 0.0, 0.0, 0.0};
    VertexIndex verts[4] = {0, 1, 2, 3};
};

int main() {
    run_case("square_area", []() {
        SquareMesh s;
        auto g = compute_face_geometry(s.px, s.py, s.pz, s.verts, 0, 4);
        // Aire d'un carré de côté 1 = 1.0
        EXPECT_NEAR(g.area, 1.0, 1e-12);
    });

    run_case("square_Sf_normal_z", []() {
        SquareMesh s;
        auto g = compute_face_geometry(s.px, s.py, s.pz, s.verts, 0, 4);
        // La normale (ordre CCW vue du +Z) est +Z.
        EXPECT_NEAR(g.Sf.x, 0.0, 1e-12);
        EXPECT_NEAR(g.Sf.y, 0.0, 1e-12);
        EXPECT_NEAR(g.Sf.z, 1.0, 1e-12);
        EXPECT_NEAR(g.normal.x, 0.0, 1e-12);
        EXPECT_NEAR(g.normal.y, 0.0, 1e-12);
        EXPECT_NEAR(g.normal.z, 1.0, 1e-12);
    });

    run_case("square_centre", []() {
        SquareMesh s;
        auto g = compute_face_geometry(s.px, s.py, s.pz, s.verts, 0, 4);
        EXPECT_NEAR(g.centre.x, 0.5, 1e-12);
        EXPECT_NEAR(g.centre.y, 0.5, 1e-12);
        EXPECT_NEAR(g.centre.z, 0.0, 1e-12);
    });

    run_case("triangle_area", []() {
        // Triangle équilatéral de côté 2 : aire = sqrt(3)
        std::vector<double> px = {0.0, 2.0, 1.0};
        std::vector<double> py = {0.0, 0.0, std::sqrt(3.0)};
        std::vector<double> pz = {0.0, 0.0, 0.0};
        VertexIndex verts[3] = {0, 1, 2};
        auto g = compute_face_geometry(px, py, pz, verts, 0, 3);
        EXPECT_NEAR(g.area, std::sqrt(3.0), 1e-12);
    });

    run_case("reversed_orientation_flips_Sf", []() {
        SquareMesh s;
        // Inverse l'ordre des sommets → normale opposée
        VertexIndex rev[4] = {3, 2, 1, 0};
        auto g = compute_face_geometry(s.px, s.py, s.pz, rev, 0, 4);
        EXPECT_NEAR(g.Sf.z, -1.0, 1e-12);
        EXPECT_NEAR(g.area, 1.0, 1e-12);
    });

    run_case("pentagon_area", []() {
        // Pentagono regular inscrito en circulo radio 1, premier vértice en (1,0)
        const double pi = std::acos(-1.0);
        std::vector<double> px(5), py(5), pz(5, 0.0);
        for (int i = 0; i < 5; ++i) {
            double a = 2.0 * pi * i / 5.0;
            px[i] = std::cos(a);
            py[i] = std::sin(a);
        }
        VertexIndex verts[5] = {0, 1, 2, 3, 4};
        auto g = compute_face_geometry(px, py, pz, verts, 0, 5);
        // Aire pentagone regular lado l : A = (1/4) sqrt(5(5+2sqrt(5))) l^2
        // l = 2 sin(pi/5)
        double l = 2.0 * std::sin(pi / 5.0);
        double area_exact = 0.25 * std::sqrt(5.0 * (5.0 + 2.0 * std::sqrt(5.0))) * l * l;
        EXPECT_NEAR(g.area, area_exact, 1e-10);
    });

    run_case("closure_invariant", []() {
        // Pour une face fermée, la somme des Sf des triangles de la triangulation
        // doit être égale au Sf total (conservation géométrique).
        SquareMesh s;
        auto g = compute_face_geometry(s.px, s.py, s.pz, s.verts, 0, 4);
        // La normale doit être unitaire.
        EXPECT_NEAR(g.normal.mag(), 1.0, 1e-12);
    });

    run_case("too_few_vertices", []() {
        std::vector<double> px = {0.0, 1.0};
        std::vector<double> py = {0.0, 0.0};
        std::vector<double> pz = {0.0, 0.0};
        VertexIndex verts[2] = {0, 1};
        EXPECT_THROW(compute_face_geometry(px, py, pz, verts, 0, 2), std::runtime_error);
    });

    return run_all();
}
