// Mesh-identity contract for the native DAT restart path (#60, #426).
//
// A native checkpoint used to be validated against the target mesh by cell
// count alone. Two different meshes with the same cell count were therefore
// interchangeable: a checkpoint written on one mesh loaded its field values
// into another without complaint. These cases pin the identity contract that
// replaces that check.
#include "cfdx/io/restart/dat_restart.h"
#include "common/test_harness.h"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::testing;

namespace {

Mesh make_cube(double scale, double dx = 0.0)
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i=0; i<8; ++i)
        m.points().set(i, p[i][0]*scale+dx, p[i][1]*scale, p[i][2]*scale);
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for (std::size_t f=0; f<6; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,1,2,3,4,5};
    m.boundary().add_patch(wall);
    return m;
}

std::filesystem::path temp_dat(const char* stem)
{
    return std::filesystem::temp_directory_path() /
        (std::string("cfdx_dat_identity_") + stem + ".dat");
}

Field<double,Location::CELL> seeded_u(double value)
{
    Field<double,Location::CELL> u(1,"U","m/s",3);
    u.set(0, value, 0.0, 0.0);
    return u;
}

Field<double,Location::CELL> seeded_p(double value)
{
    Field<double,Location::CELL> p(1,"p","Pa",1);
    p(0) = value;
    return p;
}

// A checkpoint in the pre-identity layout, written by hand because the current
// writers no longer produce one.
void write_legacy_v1_dat(const std::filesystem::path& path)
{
    std::ofstream out(path);
    out.precision(17);
    out << "CFDX-DAT 1\n";
    out << "cells 1\n";
    out << "iteration 17\n";
    out << "time 2.5\n";
    out << "field U 3\n";
    out << "0.25 -0.15 0.05\n";
    out << "field p 1\n";
    out << "37.5\n";
}

void run_cases()
{
    // The defect this file exists for: two distinct meshes, one cell each.
    run_case("dat_restart_rejects_checkpoint_from_a_different_mesh", [] {
        const Mesh source = make_cube(1.0);
        const Mesh other = make_cube(2.0);
        EXPECT_TRUE(source.n_cells() == other.n_cells());

        const auto dat = temp_dat("other_mesh");
        write_dat_restart_fields(
            dat.string(), source, seeded_u(0.25), seeded_p(37.5), {}, 17, 2.5);

        Field<double,Location::CELL> u(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        EXPECT_THROW_WITH(
            read_dat_restart_fields(dat.string(), other, u, p, {}),
            std::runtime_error, "belongs to a different mesh");

        std::filesystem::remove(dat);
    });

    run_case("dat_restart_rejects_checkpoint_written_before_the_identity_contract", [] {
        const Mesh m = make_cube(1.0);
        const auto dat = temp_dat("legacy_v1");
        write_legacy_v1_dat(dat);

        Field<double,Location::CELL> u(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        EXPECT_THROW_WITH(
            read_dat_restart_fields(dat.string(), m, u, p, {}),
            std::runtime_error, "predates the mesh identity contract");

        std::filesystem::remove(dat);
    });

    run_case("dat_restart_round_trips_the_same_mesh_identity", [] {
        const Mesh m = make_cube(1.0);
        const auto dat = temp_dat("same_mesh");
        write_dat_restart_fields(
            dat.string(), m, seeded_u(0.25), seeded_p(37.5), {}, 17, 2.5);

        Field<double,Location::CELL> u(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        const auto state = read_dat_restart_fields(dat.string(), m, u, p, {});

        EXPECT_TRUE(state.iteration == std::size_t{17});
        EXPECT_NEAR(state.time, 2.5, 1e-14);
        EXPECT_NEAR(u(0,0), 0.25, 1e-14);
        EXPECT_NEAR(p(0), 37.5, 1e-14);

        std::filesystem::remove(dat);
    });

    run_case("dat_restart_identity_distinguishes_translated_geometry", [] {
        const Mesh a = make_cube(1.0);
        const Mesh b = make_cube(1.0, 5.0);

        const auto id_a = dat_restart_identity(a);
        const auto id_b = dat_restart_identity(b);
        EXPECT_TRUE(id_a.n_cells == id_b.n_cells);
        EXPECT_TRUE(id_a.n_points == id_b.n_points);
        EXPECT_TRUE(id_a.geometry_checksum != id_b.geometry_checksum);

        const auto dat = temp_dat("translated");
        write_dat_restart_fields(dat.string(), a, seeded_u(1.0), seeded_p(2.0), {}, 1, 0.0);

        Field<double,Location::CELL> u(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        EXPECT_THROW_WITH(
            read_dat_restart_fields(dat.string(), b, u, p, {}),
            std::runtime_error, "belongs to a different mesh");

        std::filesystem::remove(dat);
    });

    run_case("dat_restart_identity_is_deterministic", [] {
        const Mesh a = make_cube(1.25);
        const Mesh b = make_cube(1.25);
        EXPECT_TRUE(dat_restart_identity(a).geometry_checksum ==
                    dat_restart_identity(b).geometry_checksum);
    });

    run_case("dat_restart_rejects_a_tampered_identity", [] {
        const Mesh m = make_cube(1.0);
        const auto dat = temp_dat("tampered");
        write_dat_restart_fields(dat.string(), m, seeded_u(1.0), seeded_p(2.0), {}, 1, 0.0);

        const auto id = dat_restart_identity(m);
        std::stringstream buffer;
        {
            std::ifstream in(dat.string());
            buffer << in.rdbuf();
        }
        std::string text = buffer.str();
        const std::string recorded = std::to_string(id.geometry_checksum);
        const auto pos = text.find(recorded);
        EXPECT_TRUE(pos != std::string::npos);
        // Flip one digit so the recorded identity no longer describes the mesh.
        std::string forged = recorded;
        forged[forged.size() - 1] = (forged.back() == '0') ? '1' : '0';
        text.replace(pos, recorded.size(), forged);
        {
            std::ofstream out(dat.string());
            out << text;
        }

        Field<double,Location::CELL> u(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa",1);
        EXPECT_THROW_WITH(
            read_dat_restart_fields(dat.string(), m, u, p, {}),
            std::runtime_error, "belongs to a different mesh");

        std::filesystem::remove(dat);
    });

    run_case("dat_restart_identity_rejects_non_finite_geometry", [] {
        Mesh broken = make_cube(1.0);
        broken.points().set(0, std::nan(""), 0.0, 0.0);
        EXPECT_THROW(dat_restart_identity(broken), std::invalid_argument);
    });
}

} // namespace

int main()
{
    run_cases();
    return run_all();
}