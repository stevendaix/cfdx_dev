#include "cfdx/physics/thermohydraulic_observables.h"
#include "common/test_harness.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {
Mesh cube() {
    Mesh mesh;
    mesh.points().resize(8);
    const double xyz[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i=0;i<8;++i)
        mesh.points().set(i, xyz[i][0], xyz[i][1], xyz[i][2]);
    mesh.faces().push_face({0,3,2,1});
    mesh.faces().push_face({4,5,6,7});
    mesh.faces().push_face({0,1,5,4});
    mesh.faces().push_face({3,7,6,2});
    mesh.faces().push_face({0,4,7,3});
    mesh.faces().push_face({1,2,6,5});
    mesh.ownership().resize(6);
    for (std::size_t f=0;f<6;++f) {
        mesh.ownership().set_owner(f,0);
        mesh.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0,1,2,3,4,5});
    mesh.boundary().add_patch("wall", PatchType::WALL);
    for (std::size_t f=0;f<6;++f)
        mesh.boundary().patch(0).face_ids.push_back(static_cast<FaceIndex>(f));
    return mesh;
}
std::string file_text(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("missing test CSV: " + path.string());
    std::ostringstream out;
    out << input.rdbuf();
    return out.str();
}
}

int main() {
    run_case("zero_net_flux_and_one_sided_wall_shear", [] {
        const auto mesh = cube();
        Field<double,Location::CELL> u(1,"U","m/s",3);
        u.set(0,2,0,0);
        Field<double,Location::CELL> p(1,"p","Pa");
        p(0)=10;
        Field<double,Location::FACE> mdot(6,"mdot","kg/s");
        mdot.fill(0);
        mdot(4)=-2;
        mdot(5)=2;
        VelocityBoundaryConditions bcs;
        bcs["wall"]={VelocityBoundaryCondition::Type::FIXED_VALUE, Vec3{0,0,0}};
        const auto dir=unique_temp_path("_th_observables");
        const auto result=write_thermohydraulic_observables(
            dir, mesh, u, p, mdot, bcs, 1, 1);
        EXPECT_NEAR(result.boundary_net_mass_flux_kg_s,0.0,1e-14);
        EXPECT_NEAR(result.global_mass_imbalance_relative,0.0,1e-14);
        EXPECT_NEAR(result.maximum_cell_mass_imbalance_relative,0.0,1e-14);
        EXPECT_TRUE(result.no_slip_wall_faces==6);
        const auto cells=file_text(dir/"cell_fields.csv");
        const auto faces=file_text(dir/"face_fluxes.csv");
        EXPECT_TRUE(cells.find("0,0.5,0.5,0.5,1,2,0,0,10,10,,,")!=std::string::npos);
        // Face 2 is y=0 with wall-normal distance 0.5 m.
        EXPECT_TRUE(faces.find("2,0,-1,\"wall\"")!=std::string::npos);
        EXPECT_TRUE(faces.find(",0.5,4,0,0,first_order_one_sided_diagnostic")!=std::string::npos);
        EXPECT_TRUE(file_text(dir/"mass_conservation.json").find(
            "\"physical_qualification\": \"NOT_QUALIFIED\"") !=std::string::npos);
        std::filesystem::remove_all(dir);
    });
    run_case("sst_modified_pressure_is_reported_as_static", [] {
        const auto mesh=cube();
        Field<double,Location::CELL> u(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa");
        Field<double,Location::FACE> mdot(6,"mdot","kg/s");
        p(0)=12;
        mdot.fill(0);
        SSTProductionState sst(1);
        sst.k(0)=3;
        sst.omega(0)=1;
        sst.nu_t(0)=0;
        sst.wall_distance(0)=0.5;
        sst.strain_rate(0)=0;
        const auto dir=unique_temp_path("_th_sst");
        write_thermohydraulic_observables(dir,mesh,u,p,mdot,{},2,0.001,&sst);
        // 12 - 2/3*2*3 = 8 Pa; failing to transform corrupts Cp.
        const auto fields=file_text(dir/"cell_fields.csv");
        EXPECT_TRUE(fields.find(",12,8,3,1,0,")!=std::string::npos);
        std::filesystem::remove_all(dir);
    });
    run_case("nonfinite_temperature_rejected", [] {
        const auto mesh=cube();
        Field<double,Location::CELL> u(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa");
        Field<double,Location::FACE> mdot(6,"mdot","kg/s");
        Field<double,Location::CELL> temperature(1,"T","K");
        temperature(0)=std::numeric_limits<double>::quiet_NaN();
        const auto dir=unique_temp_path("_th_nonfinite");
        EXPECT_THROW(write_thermohydraulic_observables(
            dir,mesh,u,p,mdot,{},1,0.001,nullptr,&temperature),std::invalid_argument);
        EXPECT_FALSE(std::filesystem::exists(dir/"cell_fields.csv"));
        EXPECT_FALSE(std::filesystem::exists(dir/"cell_fields.csv.tmp"));
        std::filesystem::remove_all(dir);
    });
    run_case("mass_balance_flags_leak", [] {
        const auto mesh=cube();
        Field<double,Location::CELL> u(1,"U","m/s",3);
        Field<double,Location::CELL> p(1,"p","Pa");
        Field<double,Location::FACE> mdot(6,"mdot","kg/s");
        mdot.fill(0);
        mdot(4)=-2;
        mdot(5)=2.2;
        const auto dir=unique_temp_path("_th_leak");
        const auto result=write_thermohydraulic_observables(
            dir,mesh,u,p,mdot,{},1,0.001);
        EXPECT_TRUE(result.global_mass_imbalance_relative > 0.09);
        EXPECT_TRUE(result.maximum_cell_mass_imbalance_relative > 0.04);
        std::filesystem::remove_all(dir);
    });
    return run_all();
}
