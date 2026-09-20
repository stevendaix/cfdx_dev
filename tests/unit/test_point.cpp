// M0.1-T01 — Tests for PointCloud
// Spécification CFDX v0.7 §10, §18, §19

#include "cfdx/core/mesh/point.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("default_empty", []() {
        PointCloud pc;
        EXPECT_TRUE(pc.size() == 0);
        EXPECT_TRUE(pc.empty());
    });

    run_case("explicit_size", []() {
        PointCloud pc(5);
        EXPECT_TRUE(pc.size() == 5);
        EXPECT_FALSE(pc.empty());
        for (std::size_t i = 0; i < 5; ++i) {
            EXPECT_TRUE(pc.x(i) == 0.0);
            EXPECT_TRUE(pc.y(i) == 0.0);
            EXPECT_TRUE(pc.z(i) == 0.0);
        }
    });

    run_case("set_get", []() {
        PointCloud pc(3);
        pc.set(0, 1.0, 2.0, 3.0);
        pc.set(1, 4.0, 5.0, 6.0);
        pc.set(2, 7.0, 8.0, 9.0);
        EXPECT_TRUE(pc.x(0) == 1.0);
        EXPECT_TRUE(pc.y(0) == 2.0);
        EXPECT_TRUE(pc.z(0) == 3.0);
        EXPECT_TRUE(pc.x(2) == 7.0);
        EXPECT_TRUE(pc.z(2) == 9.0);
    });

    run_case("index_out_of_range", []() {
        PointCloud pc(1);
        EXPECT_THROW(pc.x(5), std::out_of_range);
    });

    run_case("resize", []() {
        PointCloud pc;
        pc.resize(4);
        EXPECT_TRUE(pc.size() == 4);
        pc.resize(10);
        EXPECT_TRUE(pc.size() == 10);
    });

    run_case("bulk_access", []() {
        PointCloud pc(3);
        pc.set(0, 1.0, 2.0, 3.0);
        pc.set(1, 4.0, 5.0, 6.0);
        pc.set(2, 7.0, 8.0, 9.0);
        const double* xd = pc.x_data();
        const double* yd = pc.y_data();
        const double* zd = pc.z_data();
        EXPECT_TRUE(xd[0] == 1.0 && xd[1] == 4.0 && xd[2] == 7.0);
        EXPECT_TRUE(yd[0] == 2.0 && yd[1] == 5.0 && yd[2] == 8.0);
        EXPECT_TRUE(zd[0] == 3.0 && zd[1] == 6.0 && zd[2] == 9.0);
    });

    run_case("valid_finite", []() {
        PointCloud pc(3);
        pc.set(0, 1.0, 2.0, 3.0);
        pc.set(1, 4.0, 5.0, 6.0);
        pc.set(2, 7.0, 8.0, 9.0);
        EXPECT_TRUE(pc.is_valid());
    });

    run_case("invalid_nan", []() {
        PointCloud pc(2);
        pc.set(0, 1.0, 2.0, 3.0);
        pc.set(1, std::nan(""), 5.0, 6.0);
        EXPECT_FALSE(pc.is_valid());
    });

    run_case("invalid_inf", []() {
        PointCloud pc(2);
        pc.set(0, 1.0, 2.0, 3.0);
        pc.set(1, std::numeric_limits<double>::infinity(), 5.0, 6.0);
        EXPECT_FALSE(pc.is_valid());
    });

    return run_all();
}