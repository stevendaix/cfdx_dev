// M0.10-T04: Case-level HDF5 I/O tests
//
// Tests read_case_cfdx_h5() and write_case_cfdx_h5() for round-trip fidelity,
// Python-writer compatibility, source metadata, CaseSetup JSON, GapAnalysis
// JSON, field read/write, and error handling.

#include "cfdx/io/hdf5/case_hdf5_io.h"
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/io/hdf5/hdf5_writer.h"
#include "cfdx/io/hdf5/mini_json.h"
#include "cfdx/io/cfdx_io/io_interface.h"
#include "cfdx/io/cfdx_io/case_schema.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include "common/test_harness.h"
#include <cstdio>
#include <string>
#include <hdf5.h>

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::testing;
using namespace cfdx::io::mini_json;

// --- Helper: unit cube mesh ---
static Mesh make_unit_cube() {
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

static CaseSetup make_default_setup() {
    CaseSetup setup;
    setup.source.solver = "TestSolver";
    setup.source.version = "1.0";
    setup.source.case_path = "/tmp/test_case";
    setup.source.case_name = "test_case";
    setup.source.format = "test";
    setup.mesh_info.n_vertices = 8;
    setup.mesh_info.n_faces = 6;
    setup.mesh_info.n_cells = 1;
    setup.mesh_info.dimension = 3;
    setup.mesh_info.mesh_type = "unstructured";
    setup.mesh_info.cell_types = {"hex"};
    setup.physics_model = "incompressible_laminar";
    setup.turbulence_model = "laminar";
    setup.energy_model = "isothermal";
    setup.numerics.momentum_scheme = "second_order_upwind";
    setup.numerics.max_iterations = 500;

    MaterialSpec mat;
    mat.name = "air";
    mat.density = 1.225;
    setup.materials.push_back(mat);

    BoundarySpec bc;
    bc.patch_name = "inlet";
    bc.type = BCType::INLET;
    bc.value_type = BCValueType::FIXED;
    bc.velocity_magnitude = 10.0;
    setup.boundary_conditions.push_back(bc);

    return setup;
}

static SourceInfo make_test_source() {
    SourceInfo src;
    src.solver = "TestSolver";
    src.version = "1.0";
    src.case_path = "/tmp/test_case";
    src.case_name = "test_case";
    src.format = "test";
    return src;
}

static GapAnalysis make_test_gap() {
    GapAnalysis gap;
    gap.supported("mesh", "points", "Parsed 8 vertices");
    gap.approximated("boundary_condition", "wall", "Wall mapped with approximation");
    gap.unsupported_nonblocking("features", "turbomachinery",
                                "Not supported", "Use base configuration");
    return gap;
}

// ===========================================================================
// Tests
// ===========================================================================

int main() {

    run_case("case_hdf5_roundtrip_mesh_and_metadata", []() {
        const std::string filename = "/tmp/cfdx_case_roundtrip.h5";
        std::remove(filename.c_str());

        Mesh mesh = make_unit_cube();
        SourceInfo source = make_test_source();
        CaseSetup setup = make_default_setup();
        GapAnalysis gap = make_test_gap();

        EXPECT_TRUE(write_case_cfdx_h5(filename, mesh, source, setup, gap));

        Mesh loaded_mesh;
        SourceInfo loaded_source;
        CaseSetup loaded_setup;
        GapAnalysis loaded_gap;

        EXPECT_TRUE(read_case_cfdx_h5(filename, loaded_mesh, loaded_source, loaded_setup, loaded_gap));

        // Mesh topology
        EXPECT_TRUE(loaded_mesh.n_points() == mesh.n_points());
        EXPECT_TRUE(loaded_mesh.n_faces() == mesh.n_faces());
        EXPECT_TRUE(loaded_mesh.n_cells() == mesh.n_cells());

        for (std::size_t i = 0; i < mesh.n_points(); ++i) {
            EXPECT_TRUE(loaded_mesh.points().x(i) == mesh.points().x(i));
            EXPECT_TRUE(loaded_mesh.points().y(i) == mesh.points().y(i));
            EXPECT_TRUE(loaded_mesh.points().z(i) == mesh.points().z(i));
        }

        // Source info
        EXPECT_TRUE(loaded_source.solver == source.solver);
        EXPECT_TRUE(loaded_source.version == source.version);
        EXPECT_TRUE(loaded_source.case_path == source.case_path);

        // CaseSetup
        EXPECT_TRUE(loaded_setup.physics_model == setup.physics_model);
        EXPECT_TRUE(loaded_setup.turbulence_model == setup.turbulence_model);
        EXPECT_TRUE(loaded_setup.numerics.momentum_scheme == setup.numerics.momentum_scheme);
        EXPECT_TRUE(loaded_setup.materials.size() == setup.materials.size());
        EXPECT_TRUE(loaded_setup.materials[0].name == "air");
        EXPECT_TRUE(loaded_setup.boundary_conditions.size() == setup.boundary_conditions.size());
        EXPECT_TRUE(loaded_setup.boundary_conditions[0].type == BCType::INLET);

        // GapAnalysis
        EXPECT_TRUE(loaded_gap.findings().size() == gap.findings().size());
        EXPECT_TRUE(loaded_gap.n_supported() == gap.n_supported());
        EXPECT_TRUE(loaded_gap.n_approximated() == gap.n_approximated());

        std::remove(filename.c_str());
    });

    run_case("case_hdf5_readable_by_mesh_hdf5", []() {
        const std::string filename = "/tmp/cfdx_case_mesh_compat.h5";
        std::remove(filename.c_str());

        Mesh mesh = make_unit_cube();
        SourceInfo source = make_test_source();
        CaseSetup setup = make_default_setup();
        GapAnalysis gap = make_test_gap();

        EXPECT_TRUE(write_case_cfdx_h5(filename, mesh, source, setup, gap));

        // The C++ case writer includes integrity hashes, so read_mesh_hdf5 should work
        Mesh loaded;
        EXPECT_TRUE(read_mesh_hdf5(filename, loaded));
        EXPECT_TRUE(loaded.n_points() == mesh.n_points());
        EXPECT_TRUE(loaded.n_faces() == mesh.n_faces());
        EXPECT_TRUE(loaded.n_cells() == mesh.n_cells());

        std::remove(filename.c_str());
    });

    run_case("case_hdf5_missing_file", []() {
        Mesh mesh;
        SourceInfo source;
        CaseSetup setup;
        GapAnalysis gap;
        EXPECT_FALSE(read_case_cfdx_h5("/tmp/cfdx_does_not_exist_12345.h5", mesh, source, setup, gap));
    });

    run_case("case_hdf5_empty_mesh_roundtrip", []() {
        const std::string filename = "/tmp/cfdx_case_empty.h5";
        std::remove(filename.c_str());

        Mesh mesh;
        SourceInfo source = make_test_source();
        CaseSetup setup = make_default_setup();
        GapAnalysis gap;

        EXPECT_TRUE(write_case_cfdx_h5(filename, mesh, source, setup, gap));

        Mesh loaded_mesh;
        SourceInfo loaded_source;
        CaseSetup loaded_setup;
        GapAnalysis loaded_gap;

        EXPECT_TRUE(read_case_cfdx_h5(filename, loaded_mesh, loaded_source, loaded_setup, loaded_gap));
        EXPECT_TRUE(loaded_mesh.n_points() == 0);
        EXPECT_TRUE(loaded_mesh.n_faces() == 0);
        EXPECT_TRUE(loaded_mesh.n_cells() == 0);

        std::remove(filename.c_str());
    });

    run_case("case_hdf5_boundary_patches_roundtrip", []() {
        const std::string filename = "/tmp/cfdx_case_patches.h5";
        std::remove(filename.c_str());

        Mesh mesh = make_unit_cube();
        SourceInfo source = make_test_source();
        CaseSetup setup = make_default_setup();
        GapAnalysis gap;

        EXPECT_TRUE(write_case_cfdx_h5(filename, mesh, source, setup, gap));

        Mesh loaded;
        SourceInfo loaded_source;
        CaseSetup loaded_setup;
        GapAnalysis loaded_gap;

        EXPECT_TRUE(read_case_cfdx_h5(filename, loaded, loaded_source, loaded_setup, loaded_gap));

        EXPECT_TRUE(loaded.boundary().n_patches() == mesh.boundary().n_patches());
        for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i) {
            EXPECT_TRUE(loaded.boundary().patch(i).name == mesh.boundary().patch(i).name);
            EXPECT_TRUE(loaded.boundary().patch(i).type == mesh.boundary().patch(i).type);
            EXPECT_TRUE(loaded.boundary().patch(i).face_ids.size() == mesh.boundary().patch(i).face_ids.size());
        }

        std::remove(filename.c_str());
    });

    run_case("case_hdf5_python_compatible_no_hashes", []() {
        const std::string filename = "/tmp/cfdx_case_python_compat.h5";
        std::remove(filename.c_str());

        // Write a case file using the C++ writer
        Mesh mesh = make_unit_cube();
        SourceInfo source = make_test_source();
        CaseSetup setup = make_default_setup();
        GapAnalysis gap = make_test_gap();

        EXPECT_TRUE(write_case_cfdx_h5(filename, mesh, source, setup, gap));

        // Delete the integrity hash attributes (simulating Python-writer output)
        hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        EXPECT_TRUE(file >= 0);
        H5Adelete(file, "topology_hash");
        H5Adelete(file, "geometry_hash");
        H5Adelete(file, "mesh_hash");
        H5Fclose(file);

        // Should still read successfully (hashes absent → skip validation)
        Mesh loaded;
        SourceInfo loaded_source;
        CaseSetup loaded_setup;
        GapAnalysis loaded_gap;

        EXPECT_TRUE(read_case_cfdx_h5(filename, loaded, loaded_source, loaded_setup, loaded_gap));
        EXPECT_TRUE(loaded.n_points() == mesh.n_points());
        EXPECT_TRUE(loaded.n_faces() == mesh.n_faces());
        EXPECT_TRUE(loaded.n_cells() == mesh.n_cells());

        std::remove(filename.c_str());
    });

    run_case("case_hdf5_gap_report_findings_roundtrip", []() {
        const std::string filename = "/tmp/cfdx_case_gap.h5";
        std::remove(filename.c_str());

        Mesh mesh = make_unit_cube();
        SourceInfo source = make_test_source();
        CaseSetup setup = make_default_setup();
        GapAnalysis gap = make_test_gap();

        EXPECT_TRUE(write_case_cfdx_h5(filename, mesh, source, setup, gap));

        Mesh loaded_mesh;
        SourceInfo loaded_source;
        CaseSetup loaded_setup;
        GapAnalysis loaded_gap;

        EXPECT_TRUE(read_case_cfdx_h5(filename, loaded_mesh, loaded_source, loaded_setup, loaded_gap));

        // Check individual findings
        const auto& findings = loaded_gap.findings();
        EXPECT_TRUE(findings.size() == 3);
        if (findings.size() >= 3) {
            EXPECT_TRUE(findings[0].severity == Severity::SUPPORTED);
            EXPECT_TRUE(findings[0].category == "mesh");
            EXPECT_TRUE(findings[1].severity == Severity::APPROXIMATED);
            EXPECT_TRUE(findings[1].category == "boundary_condition");
            EXPECT_TRUE(findings[2].severity == Severity::UNSUPPORTED_NONBLOCK);
            EXPECT_TRUE(findings[2].category == "features");
        }

        std::remove(filename.c_str());
    });

    run_case("case_hdf5_case_setup_full_roundtrip", []() {
        const std::string filename = "/tmp/cfdx_case_setup.h5";
        std::remove(filename.c_str());

        Mesh mesh = make_unit_cube();
        SourceInfo source = make_test_source();
        CaseSetup setup = make_default_setup();
        GapAnalysis gap;

        EXPECT_TRUE(write_case_cfdx_h5(filename, mesh, source, setup, gap));

        CaseSetup loaded_setup;
        Mesh loaded_mesh;
        SourceInfo loaded_source;
        GapAnalysis loaded_gap;

        EXPECT_TRUE(read_case_cfdx_h5(filename, loaded_mesh, loaded_source, loaded_setup, loaded_gap));

        // Verify physics
        EXPECT_TRUE(loaded_setup.physics_model == "incompressible_laminar");
        EXPECT_TRUE(loaded_setup.turbulence_model == "laminar");
        EXPECT_TRUE(loaded_setup.energy_model == "isothermal");

        // Verify numerics
        EXPECT_TRUE(loaded_setup.numerics.momentum_scheme == "second_order_upwind");
        EXPECT_TRUE(loaded_setup.numerics.max_iterations == 500);

        // Verify mesh metadata
        EXPECT_TRUE(loaded_setup.mesh_info.n_vertices == 8);
        EXPECT_TRUE(loaded_setup.mesh_info.n_cells == 1);
        EXPECT_TRUE(loaded_setup.mesh_info.dimension == 3);
        EXPECT_TRUE(loaded_setup.mesh_info.cell_types.size() == 1);
        EXPECT_TRUE(loaded_setup.mesh_info.cell_types[0] == "hex");

        // Verify boundary conditions
        EXPECT_TRUE(loaded_setup.boundary_conditions.size() == 1);
        EXPECT_TRUE(loaded_setup.boundary_conditions[0].patch_name == "inlet");
        EXPECT_TRUE(loaded_setup.boundary_conditions[0].velocity_magnitude == 10.0);

        // Verify materials
        EXPECT_TRUE(loaded_setup.materials.size() == 1);
        EXPECT_TRUE(loaded_setup.materials[0].name == "air");
        EXPECT_TRUE(loaded_setup.materials[0].density == 1.225);

        std::remove(filename.c_str());
    });

    run_case("case_hdf5_read_scalar_fields", []() {
        const std::string filename = "/tmp/cfdx_case_scalar_fields.h5";
        std::remove(filename.c_str());

        Mesh mesh = make_unit_cube();
        SourceInfo source = make_test_source();
        CaseSetup setup = make_default_setup();
        GapAnalysis gap;

        EXPECT_TRUE(write_case_cfdx_h5(filename, mesh, source, setup, gap));

        // Write scalar fields directly
        std::vector<ScalarCellField> scalars;
        ScalarCellField pressure(1, "pressure", "Pa", 1);
        pressure(0) = 101325.0;
        scalars.push_back(pressure);

        std::vector<Vec3CellField> empty_vecs;
        EXPECT_TRUE(write_fields_hdf5(filename, scalars, empty_vecs));

        // Read back
        std::vector<std::pair<std::string, ScalarCellField>> loaded;
        EXPECT_TRUE(read_scalar_fields_hdf5(filename, loaded));
        EXPECT_TRUE(loaded.size() == 1);
        if (loaded.size() == 1) {
            EXPECT_TRUE(loaded[0].first == "pressure");
            EXPECT_TRUE(loaded[0].second.size() == 1);
            EXPECT_TRUE(loaded[0].second(0) == 101325.0);
        }

        std::remove(filename.c_str());
    });

    run_case("case_hdf5_read_vector_fields", []() {
        const std::string filename = "/tmp/cfdx_case_vec_fields.h5";
        std::remove(filename.c_str());

        Mesh mesh = make_unit_cube();
        SourceInfo source = make_test_source();
        CaseSetup setup = make_default_setup();
        GapAnalysis gap;

        EXPECT_TRUE(write_case_cfdx_h5(filename, mesh, source, setup, gap));

        // Write vector fields directly
        std::vector<ScalarCellField> empty_scalars;
        std::vector<Vec3CellField> vectors;
        Vec3CellField velocity(1, "velocity", "m/s", 3);
        velocity.set(0, 1.0, 2.0, 3.0);
        vectors.push_back(velocity);

        EXPECT_TRUE(write_fields_hdf5(filename, empty_scalars, vectors));

        // Read back
        std::vector<std::pair<std::string, Vec3CellField>> loaded;
        EXPECT_TRUE(read_vector_fields_hdf5(filename, loaded));
        EXPECT_TRUE(loaded.size() == 1);
        if (loaded.size() == 1) {
            EXPECT_TRUE(loaded[0].first == "velocity");
            EXPECT_TRUE(loaded[0].second.size() == 1);
            EXPECT_TRUE(loaded[0].second(0, 0) == 1.0);
            EXPECT_TRUE(loaded[0].second(0, 1) == 2.0);
            EXPECT_TRUE(loaded[0].second(0, 2) == 3.0);
        }

        std::remove(filename.c_str());
    });

    run_case("case_hdf5_rejects_corrupt_hashes", []() {
        const std::string filename = "/tmp/cfdx_case_bad_hash.h5";
        std::remove(filename.c_str());

        Mesh mesh = make_unit_cube();
        SourceInfo source = make_test_source();
        CaseSetup setup = make_default_setup();
        GapAnalysis gap;

        EXPECT_TRUE(write_case_cfdx_h5(filename, mesh, source, setup, gap));

        // Corrupt the topology hash
        hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        EXPECT_TRUE(file >= 0);
        hid_t attr = H5Aopen(file, "topology_hash", H5P_DEFAULT);
        EXPECT_TRUE(attr >= 0);
        if (attr >= 0) {
            hid_t type = H5Aget_type(attr);
            H5Awrite(attr, type, "0000000000000000");
            H5Tclose(type);
            H5Aclose(attr);
        }
        H5Fclose(file);

        Mesh loaded;
        SourceInfo loaded_source;
        CaseSetup loaded_setup;
        GapAnalysis loaded_gap;
        EXPECT_FALSE(read_case_cfdx_h5(filename, loaded, loaded_source, loaded_setup, loaded_gap));

        std::remove(filename.c_str());
    });

    run_case("case_hdf5_mesh_topology_attribute", []() {
        const std::string filename = "/tmp/cfdx_case_topology_attr.h5";
        std::remove(filename.c_str());

        Mesh mesh = make_unit_cube();
        SourceInfo source = make_test_source();
        CaseSetup setup = make_default_setup();
        GapAnalysis gap;

        EXPECT_TRUE(write_case_cfdx_h5(filename, mesh, source, setup, gap));

        hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        EXPECT_TRUE(file >= 0);

        std::string topo_val;
        hid_t attr = H5Aopen(file, "mesh_topology", H5P_DEFAULT);
        EXPECT_TRUE(attr >= 0);
        if (attr >= 0) {
            hid_t atype = H5Aget_type(attr);
            hsize_t sz = H5Tget_size(atype);
            if (sz == 0) sz = 1;
            std::vector<char> buf(sz + 1, '\0');
            H5Aread(attr, atype, buf.data());
            topo_val = buf.data();
            H5Tclose(atype);
            H5Aclose(attr);
        }
        EXPECT_TRUE(topo_val == "cfdx-csr-v1");

        H5Fclose(file);
        std::remove(filename.c_str());
    });

    run_case("mini_json_parse_simple_object", []() {
        std::string json = R"({"key":"value","num":42,"flag":true})";
        value v = value::parse_safe(json);
        EXPECT_TRUE(v.is_object());
        const value* key_val = mini_json::find(v, "key");
        EXPECT_TRUE(key_val != nullptr && key_val->is_string());
        EXPECT_TRUE(key_val->string == "value");
        const value* num_val = mini_json::find(v, "num");
        EXPECT_TRUE(num_val != nullptr && num_val->is_number());
        EXPECT_TRUE(num_val->number == 42.0);
        const value* flag_val = mini_json::find(v, "flag");
        EXPECT_TRUE(flag_val != nullptr && flag_val->is_boolean());
        EXPECT_TRUE(flag_val->boolean == true);
    });

    run_case("mini_json_parse_array", []() {
        std::string json = R"({"items":[1,2,3]})";
        value v = value::parse_safe(json);
        EXPECT_TRUE(v.is_object());
        const value* items = mini_json::find(v, "items");
        EXPECT_TRUE(items != nullptr);
        EXPECT_TRUE(items->is_array());
        EXPECT_TRUE(items->array.size() == 3);
        EXPECT_TRUE(items->array[0].number == 1.0);
        EXPECT_TRUE(items->array[1].number == 2.0);
        EXPECT_TRUE(items->array[2].number == 3.0);
    });

    run_case("mini_json_parse_nested_object", []() {
        std::string json = R"({"nested":{"deep":"found","count":5}})";
        value v = value::parse_safe(json);
        const value* nested = mini_json::find(v, "nested");
        EXPECT_TRUE(nested != nullptr);
        EXPECT_TRUE(nested->is_object());
        EXPECT_TRUE(mini_json::find(*nested, "deep")->string == "found");
        EXPECT_TRUE(mini_json::find(*nested, "count")->number == 5.0);
    });

    run_case("mini_json_escape_and_serialize", []() {
        value v(value_type::object);
        v.object["message"] = std::string("hello world");
        v.object["escaped"] = std::string("with \"quotes\" and \\backslash");
        v.object["num"] = mini_json::value(3.14);

        std::string serialized = v.serialize();
        value parsed = value::parse_safe(serialized);
        EXPECT_TRUE(parsed.is_object());
        EXPECT_TRUE(parsed.object.at("message").string == "hello world");
        EXPECT_TRUE(parsed.object.at("escaped").string == "with \"quotes\" and \\backslash");
        EXPECT_TRUE(parsed.object.at("num").number == 3.14);
    });

    run_case("mini_json_invalid_returns_null", []() {
        value v = value::parse_safe("not valid json");
        EXPECT_TRUE(v.is_null());
    });

    return run_all();
}
