// M0.4-T01 — Tests for Field<T, Location>

#include "cfdx/core/field/field.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("default_empty", []() {
        ScalarCellField f;
        EXPECT_TRUE(f.size() == 0);
        EXPECT_TRUE(f.empty());
        EXPECT_TRUE(f.dimension() == 1);
    });

    run_case("explicit_size_scalar", []() {
        ScalarCellField f(5, "p", "Pa", 1);
        EXPECT_TRUE(f.size() == 5);
        EXPECT_TRUE(f.dimension() == 1);
        for (std::size_t i = 0; i < 5; ++i) {
            EXPECT_TRUE(f(i) == 0.0);
        }
    });

    run_case("set_get_scalar", []() {
        ScalarCellField f(3, "p", "Pa", 1);
        f(0) = 1.0;
        f(1) = 2.0;
        f(2) = 3.0;
        EXPECT_TRUE(f(0) == 1.0);
        EXPECT_TRUE(f(1) == 2.0);
        EXPECT_TRUE(f(2) == 3.0);
    });

    run_case("fill", []() {
        ScalarCellField f(4, "p", "Pa", 1);
        f.fill(42.0);
        for (std::size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(f(i) == 42.0);
        }
    });

    run_case("index_out_of_range", []() {
        ScalarCellField f(1);
        EXPECT_THROW(f(5), std::out_of_range);
    });

    run_case("metadata", []() {
        ScalarCellField f(3, "p", "Pa", 1);
        EXPECT_TRUE(f.name() == "p");
        EXPECT_TRUE(f.loc() == Location::CELL);
        EXPECT_TRUE(f.metadata().unit == "Pa");
    });

    run_case("valid_finite", []() {
        ScalarCellField f(3, "p", "Pa", 1);
        f(0) = 1.0; f(1) = 2.0; f(2) = 3.0;
        EXPECT_TRUE(f.is_valid());
    });

    run_case("invalid_nan", []() {
        ScalarCellField f(2, "p", "Pa", 1);
        f(0) = 1.0;
        f(1) = std::nan("");
        EXPECT_FALSE(f.is_valid());
    });

    run_case("resize", []() {
        ScalarCellField f;
        f.resize(4);
        EXPECT_TRUE(f.size() == 4);
    });

    run_case("clear", []() {
        ScalarCellField f(3);
        f.clear();
        EXPECT_TRUE(f.size() == 0);
        EXPECT_TRUE(f.empty());
    });

    run_case("bulk_data_access", []() {
        ScalarCellField f(3, "p", "Pa", 1);
        f(0) = 1.0; f(1) = 2.0; f(2) = 3.0;
        const double* d = f.component_data(0);
        EXPECT_TRUE(d[0] == 1.0 && d[1] == 2.0 && d[2] == 3.0);
    });

    run_case("location_from_string", []() {
        EXPECT_TRUE(location_from_string("cell") == Location::CELL);
        EXPECT_TRUE(location_from_string("face") == Location::FACE);
        EXPECT_TRUE(location_from_string("point") == Location::POINT);
        EXPECT_TRUE(location_from_string("boundary") == Location::BOUNDARY);
        EXPECT_TRUE(location_from_string("unknown") == Location::UNKNOWN);
    });

    return run_all();
}