// M0.9-T03 — Round-trip validation for field only
#include "cfdx/io/hdf5/hdf5_writer.h"
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include "common/test_harness.h"
#include <cstdio>
#include <string>

using namespace cfdx::core;
using namespace cfdx::testing;
using namespace cfdx::io;

// Construit un champ simple
ScalarCellField make_test_field() {
    return ScalarCellField(3, "p", "Pa", 1);
}

int main() {
    run_case("field_write_read_roundtrip", []() {
        Mesh m;
        m.points().resize(8);
        m.points().set(0, 0.0, 0.0, 0.0);
        m.points().set(1, 1.0, 0.0, 0.0);
        m.points().set(2, 1.0, 1.0, 0.0);
        m.points().set(3, 0.0, 1.0, 0.0);
        m.points().set(4, 0.0, 0.0, 1.0);
        m.points().set(5, 1.0, 0.0, 1.0);
        m.points().set(6, 1.0, 1.0, 1.0);
        m.points().set(7, 0.0, 1.0, 1.0);
        m.faces().push_face({0, 3, 2, 1});
        m.faces().push_face({4, 5, 6, 7});
        m.faces().push_face({0, 1, 5, 4});
        m.faces().push_face({3, 7, 6, 2});
        m.faces().push_face({0, 4, 7, 3});
        m.faces().push_face({1, 2, 6, 5});
        m.ownership().resize(6);
        for (std::size_t i = 0; i < 6; ++i) {
            m.ownership().set_owner(i, 0);
            m.ownership().set_neighbour(i, 0);
        }
        m.cells().push_cell({0, 1, 2, 3, 4, 5});

        BoundaryPatches bp;
        const char* names[6] = {"bottom", "top", "front", "back", "left", "right"};
        const PatchType types[6] = {PatchType::WALL, PatchType::WALL, PatchType::INLET,
                                     PatchType::OUTLET, PatchType::SYMMETRY, PatchType::WALL};
        for (int i = 0; i < 6; ++i) {
            Patch p;
            p.name = names[i];
            p.type = types[i];
            p.face_ids = {static_cast<std::uint32_t>(i)};
            bp.add_patch(p);
        }
        m.set_boundary(bp);

        const std::string filename = "/tmp/cfdx_field_only.h5";
        std::remove(filename.c_str());

        ScalarCellField original = make_test_field();
        original(0) = 100.0;
        original(1) = 200.0;
        original(2) = 300.0;

        EXPECT_TRUE(write_mesh_hdf5(filename, m));
        EXPECT_TRUE(write_field_hdf5(filename, original));

        ScalarCellField loaded;
        EXPECT_TRUE(read_field_hdf5(filename, loaded));

        EXPECT_TRUE(loaded.size() == original.size());
        EXPECT_TRUE(loaded.name() == "p");
        EXPECT_TRUE(loaded.metadata().unit == "Pa");
        EXPECT_TRUE(loaded.metadata().dimension == 1);
        for (std::size_t i = 0; i < original.size(); ++i) {
            EXPECT_TRUE(loaded(i) == original(i));
        }
        std::remove(filename.c_str());
    });

    return run_all();
}
