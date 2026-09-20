// M0.9-T01 — Tests for HDF5 writer

#include "cfdx/io/hdf5/hdf5_writer.h"
#include "cfdx/core/mesh/mesh.h"
#include "test_harness.h"
#include <cstdio>
#include <string>

using namespace cfdx::core;
using namespace cfdx::testing;
using namespace cfdx::io;

// Construit un cube unité [0,1]^3 — 1 cellule, 6 faces, 8 sommets.
Mesh make_unit_cube() {
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
        m.ownership().set_neighbour(i, FaceOwnership::BOUNDARY);
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

    return m;
}

int main() {
    run_case("write_mesh_hdf5_creates_file", []() {
        Mesh m = make_unit_cube();
        const std::string filename = "/tmp/cfdx_test_mesh.h5";
        std::remove(filename.c_str());
        bool ok = write_mesh_hdf5(filename, m);
        EXPECT_TRUE(ok);
        // Le fichier doit exister.
        FILE* f = std::fopen(filename.c_str(), "rb");
        EXPECT_TRUE(f != nullptr);
        if (f) std::fclose(f);
        std::remove(filename.c_str());
    });

    run_case("write_mesh_hdf5_empty_mesh", []() {
        Mesh m;
        const std::string filename = "/tmp/cfdx_test_mesh_empty.h5";
        std::remove(filename.c_str());
        bool ok = write_mesh_hdf5(filename, m);
        EXPECT_TRUE(ok);
        std::remove(filename.c_str());
    });

    run_case("write_field_hdf5_creates_file", []() {
        Mesh m = make_unit_cube();
        const std::string filename = "/tmp/cfdx_test_field.h5";
        std::remove(filename.c_str());
        write_mesh_hdf5(filename, m);

        ScalarCellField f(1, "p", "Pa", 1);
        f(0) = 42.0;
        bool ok = write_field_hdf5(filename, f);
        EXPECT_TRUE(ok);
        std::remove(filename.c_str());
    });

    run_case("write_field_hdf5_metadata", []() {
        Mesh m = make_unit_cube();
        const std::string filename = "/tmp/cfdx_test_field_meta.h5";
        std::remove(filename.c_str());
        write_mesh_hdf5(filename, m);

        ScalarCellField f(3, "T", "K", 1);
        f(0) = 300.0; f(1) = 310.0; f(2) = 320.0;
        bool ok = write_field_hdf5(filename, f);
        EXPECT_TRUE(ok);
        std::remove(filename.c_str());
    });

    return run_all();
}