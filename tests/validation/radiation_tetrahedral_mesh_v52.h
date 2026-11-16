#pragma once
// V52 research fixture: consistently oriented Kuhn tetrahedra.
// Independent of the V48/V49 SMM assembly and its gradient reconstruction.
#include "cfdx/core/mesh/mesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

namespace cfdx::testing::radiation_v52 {
inline cfdx::core::Mesh tetrahedral_mesh(std::size_t n) {
    using namespace cfdx::core;
    if(n<3 || n>32)throw std::invalid_argument("V52 tetra: unsupported grid");
    Mesh m;
    const std::size_t nv=n+1;
    const auto vid=[nv](std::size_t i,std::size_t j,std::size_t k){
        return i+nv*(j+nv*k);
    };
    m.points().resize(nv*nv*nv);
    for(std::size_t k=0;k<=n;++k)
      for(std::size_t j=0;j<=n;++j)
        for(std::size_t i=0;i<=n;++i){
            const double x=double(i)/double(n),y=double(j)/double(n);
            const double z=double(k)/double(n);
            m.points().set(vid(i,j,k),x+0.17*y,y+0.09*z,0.7*z);
        }
    const std::array<std::array<unsigned,3>,8> corners{{
        {{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},
        {{0,0,1}},{{1,0,1}},{{1,1,1}},{{0,1,1}}}};
    const std::array<std::array<unsigned,4>,6> tets{{
        {{0,1,2,6}},{{0,1,5,6}},{{0,3,2,6}},
        {{0,3,7,6}},{{0,4,5,6}},{{0,4,7,6}}}};
    std::map<std::array<std::size_t,3>,std::size_t> faces;
    std::vector<std::size_t> owners;
    std::vector<std::int64_t> neighbours;
    const auto point=[&](std::size_t id){
        return Vec3{m.points().x(id),m.points().y(id),m.points().z(id)};
    };
    for(std::size_t k=0;k<n;++k)
      for(std::size_t j=0;j<n;++j)
        for(std::size_t i=0;i<n;++i){
            std::array<std::size_t,8> p{};
            for(unsigned q=0;q<8;++q)
                p[q]=vid(i+corners[q][0],j+corners[q][1],k+corners[q][2]);
            for(const auto& t:tets){
                const std::array<std::size_t,4> v{{p[t[0]],p[t[1]],p[t[2]],p[t[3]]}};
                const auto centre=(point(v[0])+point(v[1])+point(v[2])+point(v[3]))*0.25;
                const std::array<std::array<std::size_t,3>,4> sides{{
                    {{v[0],v[1],v[2]}},{{v[0],v[1],v[3]}},
                    {{v[0],v[2],v[3]}},{{v[1],v[2],v[3]}}}};
                const auto cid=m.n_cells();
                std::vector<std::size_t> local;
                local.reserve(4);
                for(auto s:sides){
                    auto key=s;std::sort(key.begin(),key.end());
                    auto it=faces.find(key);
                    if(it==faces.end()){
                        const auto A=point(s[0]),B=point(s[1]),C=point(s[2]);
                        if((B-A).cross(C-A).dot((A+B+C)*(1./3.)-centre)<0.)
                            std::swap(s[1],s[2]);
                        const auto fid=m.n_faces();
                        m.faces().push_face({s[0],s[1],s[2]});
                        faces.emplace(key,fid);
                        owners.push_back(cid);
                        neighbours.push_back(FaceOwnership::BOUNDARY);
                        local.push_back(fid);
                    }else{
                        const auto fid=it->second;
                        if(neighbours[fid]!=FaceOwnership::BOUNDARY || owners[fid]==cid)
                            throw std::runtime_error("V52 tetra: nonmanifold face");
                        neighbours[fid]=static_cast<std::int64_t>(cid);
                        local.push_back(fid);
                    }
                }
                m.cells().push_cell(local);
            }
        }
    m.ownership().resize(m.n_faces());
    Patch patch;patch.name="radiating_wall";patch.type=PatchType::WALL;
    for(std::size_t f=0;f<m.n_faces();++f){
        m.ownership().set_owner(f,owners[f]);
        m.ownership().set_neighbour(f,neighbours[f]);
        if(neighbours[f]==FaceOwnership::BOUNDARY)patch.face_ids.push_back(f);
    }
    if(patch.face_ids.empty())throw std::runtime_error("V52 tetra: no boundary");
    m.boundary().add_patch(patch);
    if(m.n_cells()!=6*n*n*n)throw std::runtime_error("V52 tetra: cell count");
    return m;
}
} // namespace cfdx::testing::radiation_v52
