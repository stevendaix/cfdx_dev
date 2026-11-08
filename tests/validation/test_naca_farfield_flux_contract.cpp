#include "naca_farfield_flux_contract.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <limits>
using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::validation;
static void require(bool v) {if(!v) throw std::runtime_error("NACA farfield flux contract FAIL");}
static bool rejects(const std::function<void()>& f) {try {f();}catch(const std::exception&){return true;}return false;}
int main() {
    Mesh m;
    m.points().resize(8);
    const double x[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                          {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for(std::size_t i=0;i<8;++i)m.points().set(i,x[i][0],x[i][1],x[i][2]);
    m.faces().push_face({0,3,2,1}); m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4}); m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3}); m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0;f<6;++f){m.ownership().set_owner(f,0);m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);}
    m.cells().push_cell({0,1,2,3,4,5});
    Patch far; far.name="farfield"; far.type=PatchType::INLET; far.face_ids={0,1};
    m.boundary().add_patch(far);
    Patch wall; wall.name="wall"; wall.type=PatchType::WALL; wall.face_ids={2,3,4,5};
    m.boundary().add_patch(wall);
    const auto g=build_fv_geometry(m);
    Field<double,Location::FACE> phi(m.n_faces(),"mass_flux","kg/s",1);
    for(std::size_t f=0;f<6;++f)phi(f)=0.0;
    const Vec3 ubc{0,0,1};
    for(const auto f:far.face_ids)phi(f)=g.face_area_vectors[f].z;
    auto audit=audit_naca_fixed_velocity_farfield(m,g,phi,"farfield",ubc,1.0);
    require(audit.inspected_faces==2 && audit.mismatch_linf==0.0);
    require(naca_kinematic_face_flux_matches(audit));
    phi(0)+=0.05; phi(1)-=0.05;
    audit=audit_naca_fixed_velocity_farfield(m,g,phi,"farfield",ubc,1.0);
    require(std::abs(audit.mismatch_signed)<1e-14);
    require(std::abs(audit.mismatch_l1-0.1)<1e-14);
    require(!naca_kinematic_face_flux_matches(audit));
    phi(0)=std::numeric_limits<double>::quiet_NaN();
    require(rejects([&]{(void)audit_naca_fixed_velocity_farfield(m,g,phi,"farfield",ubc,1.0);}));
    require(rejects([&]{(void)audit_naca_fixed_velocity_farfield(m,g,phi,"missing",ubc,1.0);}));
    require(rejects([&]{(void)audit_naca_fixed_velocity_farfield(m,g,phi,"farfield",ubc,-1.0);}));
    std::cout << "NACA_FARFIELD_FLUX_CONTRACT: PASS\n";
}
