// M0.5-T01 — Tests for BoundaryField / PatchField

#include "cfdx/core/boundary/boundary_field.h"
#include "cfdx/core/mesh/boundary.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("patchfield_default", []() {
        PatchField pf;
        EXPECT_TRUE(pf.size() == 0);
        EXPECT_TRUE(pf.empty());
    });

    run_case("patchfield_explicit", []() {
        PatchField pf("inlet", 5, "fixedValue");
        EXPECT_TRUE(pf.size() == 5);
        EXPECT_TRUE(pf.patch_name() == "inlet");
        EXPECT_TRUE(pf.type() == "fixedValue");
        for (std::size_t i = 0; i < 5; ++i) {
            EXPECT_TRUE(pf(i) == 0.0);
        }
    });

    run_case("patchfield_set_get", []() {
        PatchField pf("wall", 3, "zero");
        pf(0) = 1.0;
        pf(1) = 2.0;
        pf(2) = 3.0;
        EXPECT_TRUE(pf(0) == 1.0);
        EXPECT_TRUE(pf(2) == 3.0);
    });

    run_case("patchfield_robin_evaluation", []() {
        PatchField pf("wall", 1, "robin");
        pf.set_robin_coefficients(2.0, 1.0, 10.0);
        EXPECT_NEAR(pf.evaluate(0, 4.0, 0.5), 4.0, 1e-12);
        EXPECT_THROW(pf.evaluate(0, 4.0, 0.0), std::invalid_argument);
    });

    run_case("patchfield_fixed_and_zero_evaluation", []() {
        PatchField fixed("inlet", 1, "fixedValue");
        fixed(0) = 7.0;
        EXPECT_NEAR(fixed.evaluate(0, 2.0, 0.5), 7.0, 1e-12);
        PatchField zero("outlet", 1, "zeroGradient");
        EXPECT_NEAR(zero.evaluate(0, 2.0, 0.5), 2.0, 1e-12);
    });

    run_case("patchfield_fill", []() {
        PatchField pf("inlet", 4, "fixedValue");
        pf.fill(42.0);
        for (std::size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(pf(i) == 42.0);
        }
    });

    run_case("patchfield_index_out_of_range", []() {
        PatchField pf("inlet", 1);
        EXPECT_THROW(pf(5), std::out_of_range);
    });

    run_case("boundaryfield_default_empty", []() {
        BoundaryField bf;
        EXPECT_TRUE(bf.n_patches() == 0);
    });

    run_case("boundaryfield_add_patch", []() {
        BoundaryField bf;
        std::size_t idx = bf.add_patch("inlet", 5, "fixedValue");
        EXPECT_TRUE(idx == 0);
        EXPECT_TRUE(bf.n_patches() == 1);
        EXPECT_TRUE(bf.has_patch("inlet"));
        EXPECT_TRUE(bf["inlet"].size() == 5);
    });

    run_case("boundaryfield_duplicate_patch", []() {
        BoundaryField bf;
        bf.add_patch("inlet", 5, "fixedValue");
        EXPECT_THROW(bf.add_patch("inlet", 3, "zero"), std::runtime_error);
    });

    run_case("boundaryfield_find", []() {
        BoundaryField bf;
        bf.add_patch("inlet", 5, "fixedValue");
        bf.add_patch("outlet", 3, "inlet");
        EXPECT_TRUE(bf.find("inlet") == 0);
        EXPECT_TRUE(bf.find("outlet") == 1);
        EXPECT_TRUE(bf.find("wall") == 2);
    });

    run_case("boundaryfield_init_from_patches", []() {
        BoundaryPatches bp;
        bp.add_patch("inlet", PatchType::INLET);
        bp.add_patch("outlet", PatchType::OUTLET);
        bp.add_patch("wall", PatchType::WALL);

        BoundaryField bf;
        bf.init_from_patches(bp, "zero");
        EXPECT_TRUE(bf.n_patches() == 3);
        EXPECT_TRUE(bf.has_patch("inlet"));
        EXPECT_TRUE(bf.has_patch("outlet"));
        EXPECT_TRUE(bf.has_patch("wall"));
        EXPECT_TRUE(bf["inlet"].type() == "zero");
    });

    run_case("boundaryfield_consistent", []() {
        BoundaryPatches bp;
        bp.add_patch("inlet", PatchType::INLET);
        bp.add_patch("outlet", PatchType::OUTLET);

        BoundaryField bf;
        bf.init_from_patches(bp, "zero");
        EXPECT_TRUE(bf.is_consistent(bp));
    });

    run_case("boundaryfield_inconsistent_size", []() {
        BoundaryPatches bp;
        bp.add_patch("inlet", PatchType::INLET);
        bp.add_patch("outlet", PatchType::OUTLET);

        BoundaryField bf;
        bf.add_patch("inlet", 5, "zero");  // mauvaise taille
        bf.add_patch("outlet", 3, "zero");
        EXPECT_FALSE(bf.is_consistent(bp));
    });

    run_case("boundaryfield_remove_patch", []() {
        BoundaryField bf;
        bf.add_patch("inlet", 5, "fixedValue");
        bf.add_patch("outlet", 3, "zero");
        bf.remove_patch(0);
        EXPECT_TRUE(bf.n_patches() == 1);
        EXPECT_TRUE(bf.has_patch("outlet"));
    });

    run_case("boundaryfield_clear", []() {
        BoundaryField bf;
        bf.add_patch("inlet", 5, "fixedValue");
        bf.clear();
        EXPECT_TRUE(bf.n_patches() == 0);
    });

    run_case("boundaryfield_patch_data_access", []() {
        BoundaryField bf;
        bf.add_patch("inlet", 3, "fixedValue");
        bf["inlet"](0) = 10.0;
        bf["inlet"](1) = 20.0;
        bf["inlet"](2) = 30.0;
        const double* d = bf["inlet"].data();
        EXPECT_TRUE(d[0] == 10.0 && d[1] == 20.0 && d[2] == 30.0);
    });

    return run_all();
}