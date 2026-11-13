#pragma once
// Real BFS152 mesh topology extracted from original V54 HDF5 proof. Test fields
// are synthetic, not a converged CFD solution.
#include "cfdx/physics/convective_face_audit_v56.h"
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
static void check_fixture_v58(bool c,const char* why){if(!c)throw std::runtime_error(why);}
static Mesh load_bfs152_fixture_v58(const char* file) {
    std::ifstream in(file);
    check_fixture_v58(bool(in),"V56 BFS fixture missing");
    std::string title;std::getline(in,title);
    check_fixture_v58(title.find("real BFS152 mesh topology")!=std::string::npos,"wrong fixture provenance");
    std::size_t np=0,nf=0,nc=0,npatches=0;
    in>>np>>nf>>nc>>npatches;
    check_fixture_v58(np==374 && nf==634 && nc==152 && npatches==8,"wrong BFS152 fixture sizes");
    Mesh m;m.points().resize(np);
    for(std::size_t i=0;i<np;++i){double x,y,z;in>>x>>y>>z;m.points().set(i,x,y,z);}
    m.ownership().resize(nf);
    for(std::size_t f=0;f<nf;++f){
        std::size_t count;in>>count;std::vector<FaceIndex> vertex(count);
        for(auto& v:vertex)in>>v;
        m.faces().push_face(vertex);
        std::size_t owner;std::int64_t neigh;in>>owner>>neigh;
        m.ownership().set_owner(f,owner);
        m.ownership().set_neighbour(f,neigh);
    }
    for(std::size_t c=0;c<nc;++c){
        std::size_t count;in>>count;std::vector<FaceIndex> ids(count);
        for(auto& f:ids)in>>f;
        m.cells().push_cell(ids);
    }
    for(std::size_t p=0;p<npatches;++p){
        std::string name;bool empty;std::size_t count;in>>name>>empty>>count;
        Patch patch;patch.name=name;patch.type=empty?PatchType::EMPTY:PatchType::WALL;
        for(std::size_t j=0;j<count;++j){FaceIndex f;in>>f;patch.face_ids.push_back(f);}
        m.boundary().add_patch(patch);
    }
    check_fixture_v58(bool(in),"truncated BFS fixture");
    return m;
}
