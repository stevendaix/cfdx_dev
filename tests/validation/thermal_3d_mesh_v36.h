#pragma once
// Common conforming 3D hex mesh for the V36 energy, CHT and radiation tests.
// Plane interface vertices are shared geometrically across region builders.
#include "cfdx/physics/energy_solver.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
namespace cfdx::testing::thermal_v36 {
using namespace cfdx::core;
constexpr double mesh_pi=3.14159265358979323846;
inline void check_mesh(bool ok,const std::string& s) { if(!ok) throw std::runtime_error(s); }
inline Mesh make_mesh(std::size_t nx,std::size_t ny,std::size_t nz,double warp,
                      double z0=0.0,double z1=0.01) {
    if(!(z1>z0) || !std::isfinite(z0) || !std::isfinite(z1))
        throw std::invalid_argument("invalid thermal box z range");
    check_mesh(nx>=3&&ny>=3&&nz>=3&&warp>=0.0&&warp<=0.04,"invalid 3D MMS geometry");
    Mesh mesh;
    const auto pid=[&](std::size_t i,std::size_t j,std::size_t k){
        return (k*(ny+1)+j)*(nx+1)+i;
    };
    mesh.points().resize((nx+1)*(ny+1)*(nz+1));
    for(std::size_t k=0;k<=nz;++k) for(std::size_t j=0;j<=ny;++j)
        for(std::size_t i=0;i<=nx;++i) {
            const double x=double(i)/double(nx),y=double(j)/double(ny),z=double(k)/double(nz);
            const double wx=warp*std::sin(mesh_pi*x)*std::sin(mesh_pi*y)*std::sin(mesh_pi*z);
            mesh.points().set(pid(i,j,k),x+wx,y,z0+(z1-z0)*z);
        }
    std::map<std::array<std::size_t,4>,std::size_t> lookup;
    std::vector<std::size_t> owners;
    std::vector<std::int64_t> neighbours;
    std::vector<int> kinds;
    const auto add=[&](std::array<std::size_t,4> v,std::size_t c,int kind) {
        auto key=v;
        std::sort(key.begin(),key.end());
        const auto it=lookup.find(key);
        if(it!=lookup.end()) {
            const std::size_t face=it->second;
            check_mesh(neighbours[face]==FaceOwnership::BOUNDARY && owners[face]!=c,
                  "3D nonmanifold interface");
            neighbours[face]=static_cast<std::int64_t>(c);
            kinds[face]=6;
            return face;
        }
        const auto face=mesh.faces().n_faces();
        mesh.faces().push_face(std::vector<FaceIndex>(v.begin(),v.end()));
        lookup.emplace(key,face);
        owners.push_back(c);
        neighbours.push_back(FaceOwnership::BOUNDARY);
        kinds.push_back(kind);
        return face;
    };
    for(std::size_t k=0;k<nz;++k) for(std::size_t j=0;j<ny;++j)
        for(std::size_t i=0;i<nx;++i) {
            const auto c=(k*ny+j)*nx+i;
            const auto a=pid(i,j,k),b=pid(i+1,j,k),d=pid(i,j+1,k),e=pid(i+1,j+1,k);
            const auto A=pid(i,j,k+1),B=pid(i+1,j,k+1),D=pid(i,j+1,k+1),E=pid(i+1,j+1,k+1);
            const auto bottom=add({a,d,e,b},c,4);
            const auto top=add({A,B,E,D},c,5);
            const auto west=add({d,a,A,D},c,0);
            const auto east=add({b,e,E,B},c,1);
            const auto south=add({a,b,B,A},c,2);
            const auto north=add({e,d,D,E},c,3);
            mesh.cells().push_cell({bottom,top,west,east,south,north});
        }
    mesh.ownership().resize(mesh.n_faces());
    for(std::size_t f=0;f<mesh.n_faces();++f) {
        mesh.ownership().set_owner(f,owners[f]);
        mesh.ownership().set_neighbour(f,neighbours[f]);
    }
    for(int kind=0;kind<6;++kind) {
        const char* names[]={"west","east","south","north","bottom","top"};
        Patch patch;patch.name=names[kind];patch.type=PatchType::WALL;
        for(std::size_t f=0;f<mesh.n_faces();++f)
            if(kinds[f]==kind) patch.face_ids.push_back(f);
        check_mesh(!patch.face_ids.empty(),"3D missing patch");
        mesh.boundary().add_patch(patch);
    }
    check_mesh(mesh.n_cells()==nx*ny*nz,"3D cell count mismatch");
    return mesh;
}
}
