#include "cfdx/io/vtu/vtu_writer.h"
#include "common/test_harness.h"
#include <filesystem>
#include <fstream>
#include <string>

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::testing;

static Mesh make_unit_cube()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i = 0; i < 8; ++i) {
        m.points().set(i, p[i][0], p[i][1], p[i][2]);
    }
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for (std::size_t f = 0; f < 6; ++f) {
        m.ownership().set_owner(f, 0);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,1,2,3,4,5};
    m.boundary().add_patch(wall);
    return m;
}

int main()
{
    run_case("vtu_polyhedron_preserves_cell_topology_and_fields", [] {
        const Mesh m = make_unit_cube();
        ScalarCellField field(1, "cell_value", "", 1);
        field(0) = 42.5;

        const auto path = std::filesystem::temp_directory_path() / "cfdx_vtu_polyhedron_test.vtu";
        VtuWriter writer;
        EXPECT_TRUE(writer.write(path.string(), m, {{"cell_value", field}}));

        std::ifstream in(path);
        const std::string xml((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        std::filesystem::remove(path);

        EXPECT_TRUE(xml.find("NumberOfCells=\"1\"") != std::string::npos);
        EXPECT_TRUE(xml.find("types") != std::string::npos);
        EXPECT_TRUE(xml.find("4.250000000000e+01") != std::string::npos);

        // VTK_POLYHEDRON is type 42 and encodes the six cube faces as
        // [number_of_faces, number_of_points, point_ids...].
        EXPECT_TRUE(xml.find("42") != std::string::npos);
        EXPECT_TRUE(xml.find("6 4 0 3 2 1") != std::string::npos);

        // Dataset-level provenance must remain outside Piece.
        EXPECT_TRUE(writer.write(path.string(), m, {}, {}, {}, 2.5, 17, true));
        std::ifstream with_metadata(path);
        const std::string metadata_xml(
            (std::istreambuf_iterator<char>(with_metadata)),
            std::istreambuf_iterator<char>());
        std::filesystem::remove(path);

        const auto field_data = metadata_xml.find("<FieldData>");
        const auto piece = metadata_xml.find("<Piece ");
        EXPECT_TRUE(field_data != std::string::npos);
        EXPECT_TRUE(piece != std::string::npos);
        EXPECT_TRUE(field_data < piece);
        EXPECT_TRUE(metadata_xml.find("2.500000000000e+00") != std::string::npos);
        EXPECT_TRUE(metadata_xml.find("Name=\"iteration\" NumberOfTuples=\"1\" format=\"ascii\">\n    17\n") != std::string::npos);
    });

    return run_all();
}
