// M0.9-T03 — Round-trip validation

#include "cfdx/io/hdf5/hdf5_writer.h"
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/core/mesh/mesh.h"
#include "common/test_harness.h"
#include <cstdio>
#include <string>
#include <hdf5.h>

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
    run_case("hdf5_schema_and_hash_attributes", []() {
        Mesh original = make_unit_cube();
        const std::string filename = "/tmp/cfdx_schema.h5";
        std::remove(filename.c_str());
        EXPECT_TRUE(write_mesh_hdf5(filename, original));
        hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        EXPECT_TRUE(file >= 0);
        const char* names[] = {"format_version", "schema_version", "cfdx_version", "topology_hash", "mesh_hash"};
        for (const char* name : names) {
            hid_t attr = H5Aopen(file, name, H5P_DEFAULT);
            EXPECT_TRUE(attr >= 0);
            if (attr >= 0) H5Aclose(attr);
        }
        H5Fclose(file);
        std::remove(filename.c_str());
    });

    run_case("roundtrip_mesh_points", []() {
        Mesh original = make_unit_cube();
        const std::string filename = "/tmp/cfdx_roundtrip.h5";
        std::remove(filename.c_str());
        EXPECT_TRUE(write_mesh_hdf5(filename, original));

        Mesh loaded;
        EXPECT_TRUE(read_mesh_hdf5(filename, loaded));

        EXPECT_TRUE(loaded.n_points() == original.n_points());
        EXPECT_TRUE(loaded.n_faces() == original.n_faces());
        EXPECT_TRUE(loaded.n_cells() == original.n_cells());

        for (std::size_t i = 0; i < original.n_points(); ++i) {
            EXPECT_TRUE(loaded.points().x(i) == original.points().x(i));
            EXPECT_TRUE(loaded.points().y(i) == original.points().y(i));
            EXPECT_TRUE(loaded.points().z(i) == original.points().z(i));
        }
        std::remove(filename.c_str());
    });

    run_case("roundtrip_mesh_topology", []() {
        Mesh original = make_unit_cube();
        const std::string filename = "/tmp/cfdx_roundtrip_topo.h5";
        std::remove(filename.c_str());
        EXPECT_TRUE(write_mesh_hdf5(filename, original));

        Mesh loaded;
        EXPECT_TRUE(read_mesh_hdf5(filename, loaded));

        // Vérification topologique.
        auto result = loaded.topo_validate();
        if (!result.ok) {
            for (const auto& e : result.errors) {
                std::fprintf(stderr, "  error: %s\n", e.c_str());
            }
        }
        EXPECT_TRUE(result.ok);
        std::remove(filename.c_str());
    });

    run_case("roundtrip_field_values", []() {
        Mesh m = make_unit_cube();
        const std::string filename = "/tmp/cfdx_roundtrip_field.h5";
        std::remove(filename.c_str());
        EXPECT_TRUE(write_mesh_hdf5(filename, m));

        ScalarCellField original(3, "p", "Pa", 1);
        original(0) = 100.0;
        original(1) = 200.0;
        original(2) = 300.0;
        EXPECT_TRUE(write_field_hdf5(filename, original));

        ScalarCellField loaded;
        EXPECT_TRUE(read_field_hdf5(filename, loaded));

        EXPECT_TRUE(loaded.size() == original.size());
        EXPECT_TRUE(loaded(0) == original(0));
        EXPECT_TRUE(loaded(1) == original(1));
        EXPECT_TRUE(loaded(2) == original(2));
        EXPECT_TRUE(loaded.name() == "p");
        EXPECT_TRUE(loaded.metadata().unit == "Pa");
        std::remove(filename.c_str());
    });

    run_case("roundtrip_field_metadata", []() {
        Mesh m = make_unit_cube();
        const std::string filename = "/tmp/cfdx_roundtrip_meta.h5";
        std::remove(filename.c_str());
        EXPECT_TRUE(write_mesh_hdf5(filename, m));

        Field<double, Location::CELL> original(2, "U", "m/s", 3);
        original.set(0, 1.0, 2.0, 3.0);
        original.set(1, 4.0, 5.0, 6.0);
        EXPECT_TRUE(write_field_hdf5(filename, original));

        Field<double, Location::CELL> loaded;
        EXPECT_TRUE(read_field_hdf5(filename, loaded));

        EXPECT_TRUE(loaded.size() == 2);
        EXPECT_TRUE(loaded.dimension() == 3);
        EXPECT_TRUE(loaded.name() == "U");
        EXPECT_TRUE(loaded.metadata().unit == "m/s");
        double x, y, z;
        loaded.get(0, x, y, z);
        EXPECT_TRUE(x == 1.0 && y == 2.0 && z == 3.0);
        loaded.get(1, x, y, z);
        EXPECT_TRUE(x == 4.0 && y == 5.0 && z == 6.0);
        std::remove(filename.c_str());
    });

    run_case("roundtrip_empty_mesh", []() {
        Mesh original;
        const std::string filename = "/tmp/cfdx_roundtrip_empty.h5";
        std::remove(filename.c_str());
        EXPECT_TRUE(write_mesh_hdf5(filename, original));

        Mesh loaded;
        EXPECT_TRUE(read_mesh_hdf5(filename, loaded));
        EXPECT_TRUE(loaded.n_points() == 0);
        EXPECT_TRUE(loaded.n_faces() == 0);
        EXPECT_TRUE(loaded.n_cells() == 0);
        std::remove(filename.c_str());
    });

    run_case("read_mesh_hdf5_rejects_missing_core_dataset", []() {
        Mesh original = make_unit_cube();
        const std::string filename = "/tmp/cfdx_roundtrip_corrupt.h5";
        std::remove(filename.c_str());
        EXPECT_TRUE(write_mesh_hdf5(filename, original));

        hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        EXPECT_TRUE(file >= 0);
        EXPECT_TRUE(H5Ldelete(file, "owner", H5P_DEFAULT) >= 0);
        H5Fclose(file);

        Mesh loaded;
        EXPECT_FALSE(read_mesh_hdf5(filename, loaded));
        std::remove(filename.c_str());
    });

    run_case("read_mesh_hdf5_missing_file", []() {
        Mesh m;
        EXPECT_FALSE(read_mesh_hdf5("/tmp/cfdx_nonexistent_file.h5", m));
    });

    run_case("read_field_hdf5_missing_group", []() {
        Mesh m = make_unit_cube();
        const std::string filename = "/tmp/cfdx_roundtrip_badgroup.h5";
        std::remove(filename.c_str());
        write_mesh_hdf5(filename, m);
        ScalarCellField f;
        EXPECT_FALSE(read_field_hdf5(filename, f));
        std::remove(filename.c_str());
    });

    return run_all();
}