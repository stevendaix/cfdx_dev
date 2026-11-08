#include "naca_paired_restart.h"
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace cfdx::core;
using namespace cfdx::validation;
namespace fs=std::filesystem;
static void check(bool value) {
    if (!value) throw std::runtime_error("paired NACA restart test failed");
}
static bool rejects(const std::function<void()>& f) {
    try { f(); } catch(const std::exception&) { return true; } return false;
}
int main() {
    const auto folder=fs::temp_directory_path()/"cfdx_paired_checkpoint_contract";
    fs::remove_all(folder);
    fs::create_directories(folder);
    const auto manifest=(folder/"last.manifest").string();
    try {
        Mesh mesh;
        mesh.points().resize(8);
        const double pts[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                                {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
        for (std::size_t i=0;i<8;++i) mesh.points().set(i,pts[i][0],pts[i][1],pts[i][2]);
        mesh.faces().push_face({0,3,2,1}); mesh.faces().push_face({4,5,6,7});
        mesh.faces().push_face({0,1,5,4}); mesh.faces().push_face({3,7,6,2});
        mesh.faces().push_face({0,4,7,3}); mesh.faces().push_face({1,2,6,5});
        mesh.ownership().resize(6);
        for (std::size_t f=0;f<6;++f) {
            mesh.ownership().set_owner(f,0);
            mesh.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
        }
        mesh.cells().push_cell({0,1,2,3,4,5});
        const auto geom=cfdx::physics::build_fv_geometry(mesh);
        const std::vector<NacaStateRow> state={
            {geom.cell_centres[0],1.0,2.0,3.0}};
        Field<double,Location::FACE> phi(mesh.n_faces(),"mass_flux","kg/s",1);
        for(std::size_t f=0;f<phi.size();++f) phi(f)=.0025*double(f+1);
        naca_write_paired_checkpoint_atomic(manifest,10,state,mesh,geom,phi,0.42,0.99);
        const auto first=naca_read_paired_checkpoint(manifest,mesh,geom);
        check(first.accepted_iteration==10 && first.cells.size()==1 && first.mass_flux.size()==6);
        for(std::size_t f=0;f<phi.size();++f) check(first.mass_flux[f]==phi(f));
        check(rejects([&]{
            auto altered=geom; altered.face_area_vectors[0].x+=1;
            (void)naca_read_paired_checkpoint(manifest,mesh,altered);
        }));
        // Publish a second accepted iterate, then corrupt its immutable phi:
        // loading MUST fail, never fall back to U-derived interpolation.
        naca_write_paired_checkpoint_atomic(manifest,20,state,mesh,geom,phi,0.42,0.99);
        std::ifstream in(manifest);
        std::string line; std::getline(in,line); std::getline(in,line);
        auto fields=naca_csv_fields(line);
        check(fields.size()==9 && fields[0]=="20");
        const auto flux_path=folder/fields[6];
        { std::ofstream out(flux_path,std::ios::app); out << "\ntampered"; }
        check(rejects([&]{(void)naca_read_paired_checkpoint(manifest,mesh,geom);}));
        // A new complete manifest atomically restores loadability.
        naca_write_paired_checkpoint_atomic(manifest,30,state,mesh,geom,phi,0.42,0.99);
        check(naca_read_paired_checkpoint(manifest,mesh,geom).accepted_iteration==30);
        // Attack: a manifest pointing outside the directory is rejected.
        { std::ofstream out(manifest); out << "CFDX_NACA_PAIRED_V2\n30,1,6,0.42,0.99,../escape.csv,foo.csv,000,000\n"; }
        check(rejects([&]{(void)naca_read_paired_checkpoint(manifest,mesh,geom);}));
    } catch(...) {fs::remove_all(folder); throw;}
    fs::remove_all(folder);
    std::cout << "NACA_PAIRED_RESTART: PASS\n";
}
