// M0.1-T01 — Tests for PointCloud
// Spécification CFDX v0.7 §10, §18, §19

#include "cfdx/core/mesh/point.h"
#include "test_harness.h"
#include <cmath>
#include <limits>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    // --- Test 1: default empty ---
    run_case("default_empty", []() {
        PointCloud pc;
        bool ok = (pc.size() == 0) && pc.empty();
        EXPECT(ok);
    });

    // --- Test 2: explicit size ---
    run_case("explicit_size", []() {
        PointCloud pc(5);
        bool ok = (pc.size() == 5) && !pc.empty();
        for (std::size_t i = 0; ok && i < 5; ++i)
            ok = (pc.x(i) == 0.0 && pc.y(i) == 0.0 && pc.z(i) == 0.0);
        EXPECT(ok);
    });

    // --- Test 3: set/get ---
    run_case("set_get", []() {
        PointCloud pc(3);
        pc.set(0, 1.0, 2.0, 3.0);
        pc.set(1, 4.0, 5.0, 6.0);
        pc.set(2, 7.0, 8.0, 9.0);
        bool ok = (pc.x(0) == 1.0 && pc.y(0) == 2.0 && pc.z(0) == 3.0);
        ok = ok && (pc.x(2) == 7.0 && pc.z(2) == 9.0);
        EXPECT(ok);
    });

    // --- Test 4: index out of range ---
    run_case("index_out_of_range", []() {
        PointCloud pc(1);
        bool threw = false;
        try { pc.x(5); }
        catch (const std::out_of_range&) { threw = true; }
        EXPECT(threw);
    });

    // --- Test 5: resize ---
    run_case("resize", []() {
        PointCloud pc;
        pc.resize(4);
        bool ok = (pc.size() == 4);
        pc.resize(10);
        ok = ok && (pc.size() == 10);
        EXPECT(ok);
    });

    // --- Test 6: bulk access ---
    run_case("bulk_access", []() {
        PointCloud pc(3);
        pc.set(0, 1.0, 2.0, 3.0);
        pc.set(1, 4.0, 5.0, 6.0);
        pc.set(2, 7.0, 8.0, 9.0);
        const double* xd = pc.x_data();
        const double* yd = pc.y_data();
        const double* zd = pc.z_data();
        bool ok = (xd[0] == 1.0 && xd[1] == 4.0 && xd[2] == 7.0);
        ok = ok && (yd[0] == 2.0 && yd[1] == 5.0 && yd[2] == 8.0);
        ok = ok && (zd[0] == 3.0 && zd[1] == 6.0 && zd[2] == 9.0);
        EXPECT(ok);
    });

    // --- Test 7: valid finite ---
    run_case("valid_finite", []() {
        PointCloud pc(3);
        pc.set(0, 1.0, 2.0, 3.0);
        pc.set(1, 4.0, 5.0, 6.0);
        pc.set(2, 7.0, 8.0, 9.0);
        EXPECT(pc.is_valid());
    });

    // --- Test 8: invalid NaN ---
    run_case("invalid_nan", []() {
        PointCloud pc(2);
        pc.set(0, 1.0, 2.0, 3.0);
        pc.set(1, std::nan(""), 5.0, 6.0);
        EXPECT(!pc.is_valid());
    });

    // --- Test 9: invalid Inf ---
    run_case("invalid_inf", []() {
        PointCloud pc(2);
        pc.set(0, 1.0, 2.0, 3.0);
        pc.set(1, std::numeric_limits<double>::infinity(), 5.0, 6.0);
        EXPECT(!pc.is_valid());
    });

    return run_all();
}