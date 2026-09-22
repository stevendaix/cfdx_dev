// M0.8-T02 — Tests for Vector

#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("vector_default", []() {
        Vector v;
        EXPECT_TRUE(v.size() == 0);
        EXPECT_TRUE(v.empty());
    });

    run_case("vector_explicit", []() {
        Vector v(3, 0.0);
        EXPECT_TRUE(v.size() == 3);
        for (std::size_t i = 0; i < 3; ++i) {
            EXPECT_TRUE(v(i) == 0.0);
        }
    });

    run_case("vector_set_get", []() {
        Vector v(3);
        v(0) = 1.0;
        v(1) = 2.0;
        v(2) = 3.0;
        EXPECT_TRUE(v(0) == 1.0);
        EXPECT_TRUE(v(1) == 2.0);
        EXPECT_TRUE(v(2) == 3.0);
    });

    run_case("vector_index_out_of_range", []() {
        Vector v(1);
        EXPECT_THROW(v(5), std::out_of_range);
    });

    run_case("vector_fill", []() {
        Vector v(3);
        v.fill(42.0);
        for (std::size_t i = 0; i < 3; ++i) {
            EXPECT_TRUE(v(i) == 42.0);
        }
    });

    run_case("vector_norm1", []() {
        Vector v(3);
        v(0) = 1.0;
        v(1) = -2.0;
        v(2) = 3.0;
        EXPECT_TRUE(v.norm1() == 6.0);
    });

    run_case("vector_norm2", []() {
        Vector v(3);
        v(0) = 3.0;
        v(1) = 4.0;
        v(2) = 0.0;
        EXPECT_TRUE(v.norm2() == 5.0);
    });

    run_case("vector_norm_inf", []() {
        Vector v(3);
        v(0) = 1.0;
        v(1) = -5.0;
        v(2) = 3.0;
        EXPECT_TRUE(v.norm_inf() == 5.0);
    });

    run_case("vector_dot", []() {
        Vector a(3);
        a(0) = 1.0; a(1) = 2.0; a(2) = 3.0;
        Vector b(3);
        b(0) = 4.0; b(1) = 5.0; b(2) = 6.0;
        EXPECT_TRUE(a.dot(b) == 32.0);
    });

    run_case("vector_dot_dimension_mismatch", []() {
        Vector a(3);
        Vector b(2);
        EXPECT_THROW(a.dot(b), std::runtime_error);
    });

    run_case("vector_add", []() {
        Vector a(2);
        a(0) = 1.0; a(1) = 2.0;
        Vector b(2);
        b(0) = 3.0; b(1) = 4.0;
        Vector c = a + b;
        EXPECT_TRUE(c(0) == 4.0);
        EXPECT_TRUE(c(1) == 6.0);
    });

    run_case("vector_subtract", []() {
        Vector a(2);
        a(0) = 5.0; a(1) = 7.0;
        Vector b(2);
        b(0) = 2.0; b(1) = 3.0;
        Vector c = a - b;
        EXPECT_TRUE(c(0) == 3.0);
        EXPECT_TRUE(c(1) == 4.0);
    });

    run_case("vector_scale", []() {
        Vector a(2);
        a(0) = 1.0; a(1) = 2.0;
        Vector b = a * 3.0;
        EXPECT_TRUE(b(0) == 3.0);
        EXPECT_TRUE(b(1) == 6.0);
    });

    run_case("vector_is_valid", []() {
        Vector v(2);
        v(0) = 1.0; v(1) = 2.0;
        EXPECT_TRUE(v.is_valid());
    });

    run_case("vector_is_invalid_nan", []() {
        Vector v(2);
        v(0) = 1.0; v(1) = std::nan("");
        EXPECT_FALSE(v.is_valid());
    });

    run_case("vector_resize", []() {
        Vector v;
        v.resize(5);
        EXPECT_TRUE(v.size() == 5);
    });

    run_case("vector_clear", []() {
        Vector v(3);
        v.clear();
        EXPECT_TRUE(v.size() == 0);
        EXPECT_TRUE(v.empty());
    });

    run_case("vector_dimension_mismatch_addition_and_subtraction", []() {
        Vector a(3);
        Vector b(2);
        EXPECT_THROW(a += b, std::runtime_error);
        EXPECT_THROW(a -= b, std::runtime_error);
        EXPECT_THROW(a + b, std::runtime_error);
        EXPECT_THROW(a - b, std::runtime_error);
    });

    run_case("vector_norms_are_consistent", []() {
        Vector v(3);
        v(0) = -3.0;
        v(1) = 4.0;
        v(2) = -12.0;
        EXPECT_NEAR(v.norm1(), 19.0, 1e-14);
        EXPECT_NEAR(v.norm2(), 13.0, 1e-14);
        EXPECT_NEAR(v.norm_inf(), 12.0, 1e-14);
        EXPECT_NEAR(v.dot(v), 169.0, 1e-14);
    });

    return run_all();
}