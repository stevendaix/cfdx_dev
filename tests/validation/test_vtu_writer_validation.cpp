#include "cfdx/io/vtu/vtu_writer.h"
#include "common/test_harness.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

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
        EXPECT_TRUE(xml.find("4.250000000000e+01") != std::string::npos);
        EXPECT_TRUE(xml.find("Name=\"types\"") != std::string::npos);
        const auto types = xml.find("Name=\"types\"");
        EXPECT_TRUE(types != std::string::npos);
        EXPECT_TRUE(xml.find("42", types) != std::string::npos);
        const auto connectivity = xml.find("Name=\"connectivity\"");
        EXPECT_TRUE(connectivity != std::string::npos);
        const auto connectivity_end = xml.find("</DataArray>", connectivity);
        EXPECT_TRUE(connectivity_end != std::string::npos);
        const auto connectivity_xml = xml.substr(connectivity, connectivity_end - connectivity);
        for (int point = 0; point < 8; ++point) {
            EXPECT_TRUE(connectivity_xml.find(std::to_string(point)) != std::string::npos);
        }
        EXPECT_TRUE(xml.find("6 4 0 3 2 1 4 4 5 6 7") != std::string::npos);
        EXPECT_TRUE(xml.find("Name=\"faceoffsets\"") != std::string::npos);

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

    // `fields_face` used to be accepted and discarded, so no face quantity
    // could ever reach the file. It is documented as being interpolated to cell
    // centres, and that is now what happens.
    run_case("vtu_face_field_is_interpolated_to_cell_data", [] {
        const Mesh m = make_unit_cube();
        ScalarFaceField face_field(6, "face_value", "", 1);
        // All six faces of the cube have area 1, so the area-weighted mean of a
        // constant field is that constant.
        for (std::size_t f = 0; f < 6; ++f) face_field(f) = 3.25;

        const auto path = std::filesystem::temp_directory_path() / "cfdx_vtu_face_field_test.vtu";
        VtuWriter writer;
        EXPECT_TRUE(writer.write(path.string(), m, {}, {{"face_value", face_field}}));

        std::ifstream in(path);
        const std::string xml((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        std::filesystem::remove(path);

        EXPECT_TRUE(xml.find("Name=\"face_value\"") != std::string::npos);
        EXPECT_TRUE(xml.find("3.250000000000e+00") != std::string::npos);
        // It must land in CellData, which is the section a reader consumes.
        const auto cell_data = xml.find("<CellData>");
        const auto name_at = xml.find("Name=\"face_value\"");
        EXPECT_TRUE(cell_data != std::string::npos);
        EXPECT_TRUE(name_at != std::string::npos);
        EXPECT_TRUE(cell_data < name_at);
    });

    run_case("vtu_face_field_is_area_weighted_across_the_cell_faces", [] {
        // A 2x1x1 box, so the six faces have two distinct areas: faces 0..3 are
        // 2x1 (area 2) and faces 4,5 are 1x1 (area 1). A uniform average would
        // give 4/6, an area-weighted one gives 8/10, so the two are separable.
        Mesh m = make_unit_cube();
        m.points().set(1, 2.0, 0.0, 0.0);
        m.points().set(2, 2.0, 1.0, 0.0);
        m.points().set(5, 2.0, 0.0, 1.0);
        m.points().set(6, 2.0, 1.0, 1.0);

        ScalarFaceField face_field(6, "w", "", 1);
        for (std::size_t f = 0; f < 6; ++f) face_field(f) = (f < 4) ? 1.0 : 0.0;

        const auto path = std::filesystem::temp_directory_path() / "cfdx_vtu_face_weight_test.vtu";
        VtuWriter writer;
        EXPECT_TRUE(writer.write(path.string(), m, {}, {{"w", face_field}}));

        std::ifstream in(path);
        const std::string xml((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        std::filesystem::remove(path);

        const auto name_at = xml.find("Name=\"w\"");
        EXPECT_TRUE(name_at != std::string::npos);
        const auto data_at = xml.find(">", name_at);
        EXPECT_TRUE(data_at != std::string::npos);
        const auto close = xml.find("</DataArray>", data_at);
        EXPECT_TRUE(close != std::string::npos);
        std::istringstream values(xml.substr(data_at + 1, close - data_at - 1));
        std::vector<double> emitted;
        double parsed = 0.0;
        while (values >> parsed) emitted.push_back(parsed);

        EXPECT_TRUE(emitted.size() == 1);
        EXPECT_NEAR(emitted[0], 0.8, 1e-12);
        // Guards the specific failure mode: a uniform mean would pass a
        // constant-field test and silently misplace every varying field.
        EXPECT_TRUE(std::abs(emitted[0] - 4.0 / 6.0) > 1e-6);
    });

    run_case("vtu_rejects_a_face_field_of_the_wrong_length", [] {
        const Mesh m = make_unit_cube();
        ScalarFaceField face_field(5, "bad", "", 1);

        const auto path = std::filesystem::temp_directory_path() / "cfdx_vtu_face_bad.vtu";
        VtuWriter writer;
        // Silently dropping a mis-sized field, as the cell path does, would
        // report a successful export that lost its data.
        EXPECT_THROW(writer.write(path.string(), m, {}, {{"bad", face_field}}), std::runtime_error);
        std::filesystem::remove(path);
    });

    run_case("vtu_rejects_a_non_finite_face_field", [] {
        const Mesh m = make_unit_cube();
        ScalarFaceField face_field(6, "nan_field", "", 1);
        face_field(2) = std::nan("");

        const auto path = std::filesystem::temp_directory_path() / "cfdx_vtu_face_nan.vtu";
        VtuWriter writer;
        EXPECT_THROW(writer.write(path.string(), m, {}, {{"nan_field", face_field}}), std::runtime_error);
        std::filesystem::remove(path);
    });

    run_case("vtu_cell_and_face_fields_coexist", [] {
        const Mesh m = make_unit_cube();
        ScalarCellField cell_field(1, "cell_value", "", 1);
        cell_field(0) = 1.0;
        ScalarFaceField face_field(6, "face_value", "", 1);
        for (std::size_t f=0; f<6; ++f) face_field(f) = 2.0;

        const auto path = std::filesystem::temp_directory_path() / "cfdx_vtu_both.vtu";
        VtuWriter writer;
        EXPECT_TRUE(writer.write(
            path.string(), m, {{"cell_value", cell_field}}, {{"face_value", face_field}}));

        std::ifstream in(path);
        const std::string xml((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        std::filesystem::remove(path);

        EXPECT_TRUE(xml.find("Name=\"cell_value\"") != std::string::npos);
        EXPECT_TRUE(xml.find("Name=\"face_value\"") != std::string::npos);
        // Both must appear inside the same CellData section.
        const auto cell_data = xml.find("<CellData>");
        const auto cell_at = xml.find("Name=\"cell_value\"");
        const auto face_at = xml.find("Name=\"face_value\"");
        const auto close_at = xml.find("</CellData>");
        EXPECT_TRUE(cell_data < cell_at);
        EXPECT_TRUE(cell_at < face_at);
        EXPECT_TRUE(face_at < close_at);
    });

    return run_all();
}
