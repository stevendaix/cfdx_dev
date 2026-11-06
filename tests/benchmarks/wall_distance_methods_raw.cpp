// Independent analytical reference on cube/rotated cube/shifted plane.
// Raw algorithms ONLY; never invokes compute_wall_distance_accuracy_checked.
// Compile: c++ -std=c++17 -O2 -Isrc tests/benchmarks/wall_distance_methods_raw.cpp -o wall_distance_methods_raw
#include "cfdx/physics/wall_distance.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
using namespace cfdx::physics;
constexpr double PI=3.14159265358979323846;
static WallSurface cube(){WallSurface s;s.points={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};const std::array<std::array<std::size_t,4>,6> f={{{{0,3,2,1}},{{4,5,6,7}},{{0,1,5,4}},{{1,2,6,5}},{{2,3,7,6}},{{3,0,4,7}}}};for(const auto&x:f){s.triangles.push_back({{x[0],x[1],x[2]}});s.triangles.push_back({{x[0],x[2],x[3]}});}return s;}
static WallDistanceVec3 rotate(WallDistanceVec3 p,double a){double x=p.x-.5,y=p.y-.5;return {.5+std::cos(a)*x-std::sin(a)*y,.5+std::sin(a)*x+std::cos(a)*y,p.z};}
static double analytic_cube(WallDistanceVec3 p){double dx=std::max({0.,-p.x,p.x-1.});double dy=std::max({0.,-p.y,p.y-1.});double dz=std::max({0.,-p.z,p.z-1.});return std::sqrt(dx*dx+dy*dy+dz*dz);}
int main(){
 std::cout<<"geometry,n,method,l2_percent,linf_absolute,invalid,converged,pde_residual,wall_seeds,outer_fixed,geometric_exact\n";
 for(auto geo:{std::string("cube"),std::string("rotated_cube"),std::string("shifted_plane")}) for(int n:{9,17,25}){
 double h=2./(n-1);auto s=cube();double ang=PI/6;
 if(geo=="rotated_cube")for(auto &p:s.points)p=rotate(p,ang);
 if(geo=="shifted_plane") {s.points={{0.13,-9,-9},{0.13,9,-9},{0.13,9,9},{0.13,-9,9}};s.triangles={{{0,1,2}},{{0,2,3}}};}
 auto trans=[&](WallDistanceVec3 p){return geo=="rotated_cube" ? rotate(p,-ang) : p;};
 auto inside=[&](WallDistanceVec3 p){if(geo=="shifted_plane")return p.x<.13; p=trans(p);return p.x>0&&p.x<1&&p.y>0&&p.y<1&&p.z>0&&p.z<1;};
 auto exact=[&](WallDistanceVec3 p){if(geo=="shifted_plane")return std::abs(p.x-.13);return analytic_cube(trans(p));};
 auto g=make_wall_distance_grid(n,n,n,{-.5,-.5,-.5},{h,h,h},inside);
 for(auto method:{WallDistanceMethod::EXACT_GEOMETRIC,WallDistanceMethod::SEARCH_BASED,WallDistanceMethod::MESH_WAVE,WallDistanceMethod::DIRECTIONAL_MESH_WAVE,WallDistanceMethod::POISSON,WallDistanceMethod::EIKONAL,WallDistanceMethod::EIKONAL_SECOND_ORDER,WallDistanceMethod::HAMILTON_JACOBI,WallDistanceMethod::ADVECTION_DIFFUSION,WallDistanceMethod::HYBRID_POISSON_EIKONAL}){
 auto r=compute_wall_distance(method,s,g,1000);
 double sum=0,tot=0,linf=0;size_t invalid=0;
 for(size_t i=0;i<g.points.size();++i)if(!g.solid[i]){auto ex=exact(g.points[i]);if(!std::isfinite(r.distance[i])){++invalid;continue;}double e=r.distance[i]-ex;sum+=e*e;tot+=ex*ex;linf=std::max(linf,std::abs(e));}
 std::cout<<geo<<","<<n<<","<<wall_distance_method_name(method)<<","<<std::setprecision(12)<<100*std::sqrt(sum/tot)<<","<<linf<<","<<invalid<<","<<r.converged<<","<<r.residual_inf<<","<<r.exact_wall_seed_samples<<","<<r.exact_outer_boundary_samples<<","<<r.geometrically_exact_on_input_mesh<<"\n";
 }
 }
}
