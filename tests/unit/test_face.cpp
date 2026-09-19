// M0.1-T02 — Tests for FaceConnectivity (CSR)

#include "cfdx/core/mesh/face.h"
#include "test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("empty", []() {
        FaceConnectivity fc;
        EXPECT_TRUE(fc.n_faces() == 0);
        EXPECT_TRUE(fc.n_vertices() == 0);
        EXPECT_TRUE(fc.is_consistent());
    });

    run_case("push_single_face", []() {
        FaceConnectivity fc;
        fc.push_face({0, 1, 2, 3});
        EXPECT_TRUE(fc.n_faces() == 1);
        EXPECT_TRUE(fc.n_vertices() == 4);
        EXPECT_TRUE(fc.face_size(0) == 4);
        EXPECT_TRUE(fc.face_offset(0) == 0);
        EXPECT_TRUE(fc.is_consistent());
    });

    run_case("push_multiple_faces", []() {
        FaceConnectivity fc;
        fc.push_face({0, 1, 2, 3});
        fc.push_face({4, 5, 6});
        fc.push_face({7, 8, 9, 10, 11});

        EXPECT_TRUE(fc.n_faces() == 3);
        EXPECT_TRUE(fc.n_vertices() == 12);
        EXPECT_TRUE(fc.face_size(0) == 4);
        EXPECT_TRUE(fc.face_size(1) == 3);
        EXPECT_TRUE(fc.face_size(2) == 5);
        EXPECT_TRUE(fc.face_offset(0) == 0);
        EXPECT_TRUE(fc.face_offset(1) == 4);
        EXPECT_TRUE(fc.face_offset(2) == 7);
        EXPECT_TRUE(fc.is_consistent());

        const auto& v = fc.vertices();
        EXPECT_TRUE(v[0] == 0 && v[4] == 4 && v[7] == 7);
    });

    run_case("indices_valid", []() {
        FaceConnectivity fc;
        fc.push_face({0, 1, 2, 3});
        fc.push_face({4, 5, 6});
        EXPECT_TRUE(fc.indices_valid(10));
        EXPECT_FALSE(fc.indices_valid(5));
    });

    run_case("face_index_out_of_range", []() {
        FaceConnectivity fc;
        fc.push_face({0, 1, 2});
        EXPECT_THROW(fc.face_size(5), std::out_of_range);
    });

    run_case("reserve", []() {
        FaceConnectivity fc;
        fc.reserve(100, 400);
        fc.push_face({0, 1, 2, 3});
        EXPECT_TRUE(fc.n_faces() == 1);
        EXPECT_TRUE(fc.is_consistent());
    });

    run_case("bulk_data_access", []() {
        FaceConnectivity fc;
        fc.push_face({0, 1, 2, 3});
        fc.push_face({4, 5, 6});

        const auto* vd = fc.vertices_data();
        const auto* od = fc.offsets_data();
        EXPECT_TRUE(vd[0] == 0 && vd[4] == 4);
        EXPECT_TRUE(od[0] == 0 && od[1] == 4 && od[2] == 7);
    });

    run_case("clear", []() {
        FaceConnectivity fc;
        fc.push_face({0, 1, 2});
        fc.clear();
        EXPECT_TRUE(fc.n_faces() == 0);
        EXPECT_TRUE(fc.n_vertices() == 0);
        EXPECT_TRUE(fc.is_consistent());
    });

    return run_all();
}