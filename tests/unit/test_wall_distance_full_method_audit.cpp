#include "cfdx/physics/wall_distance_precision.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cfdx::physics;

static void check(bool ok,const char* what) {
    if(!ok) throw std::runtime_error(what);
}
static double seg2(WallDistanceVec3 p,WallDistanceVec3 a,WallDistanceVec3 b) {
    auto e=b-a;
    const auto r=wd_norm2(e);
    const auto t=r>0?std::clamp(wd_dot(p-a,e)/r,0.0,1.0):0.0;
    return wd_norm2(p-(a+e*t));
}
static WallSurface cube(double radians=0) {
    WallSurface s;
    s.points={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
              {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    const std::array<std::array<std::size_t,4>,6> faces={{
        {{0,3,2,1}},{{4,5,6,7}},{{0,1,5,4}},
        {{1,2,6,5}},{{2,3,7,6}},{{3,0,4,7}}}};
    for(const auto& f:faces) {
        s.triangles.push_back({{f[0],f[1],f[2]}});
        s.triangles.push_back({{f[0],f[2],f[3]}});
    }
    for(auto& p:s.points) {
        const double x=p.x-.5,y=p.y-.5;
        p.x=.5+std::cos(radians)*x-std::sin(radians)*y;
        p.y=.5+std::sin(radians)*x+std::cos(radians)*y;
    }
    return s;
}
static WallDistanceVec3 unrotate(WallDistanceVec3 p,double radians) {
    const double x=p.x-.5,y=p.y-.5;
    return {.5+std::cos(radians)*x+std::sin(radians)*y,
            .5-std::sin(radians)*x+std::cos(radians)*y,p.z};
}
static double cube_analytic(WallDistanceVec3 p) {
    const double dx=std::max({0.0,-p.x,p.x-1.0});
    const double dy=std::max({0.0,-p.y,p.y-1.0});
    const double dz=std::max({0.0,-p.z,p.z-1.0});
    return std::sqrt(dx*dx+dy*dy+dz*dz);
}
static double l2_error(const WallDistanceGrid& g,const WallDistanceResult& r,
                       const std::function<double(WallDistanceVec3)>& oracle) {
    check(r.distance.size()==g.points.size(),"result must have one value per grid point");
    double err2=0.0,ref2=0.0;
    for(std::size_t i=0;i<g.points.size();++i) if(!g.solid[i]) {
        check(std::isfinite(r.distance[i]),"invalid fluid distance");
        const double exact=oracle(g.points[i]);
        const double e=r.distance[i]-exact;
        err2+=e*e; ref2+=exact*exact;
    }
    return std::sqrt(err2/ref2);
}
int main() {
    // Real pre-patch failure: collinear triangle reported a distance about
    // 0.0158 (squared units) above the closest segment in this configuration.
    const WallDistanceVec3 a{-4.5146360758192197,-0.74922346056342448,-3.4261314561376923};
    const WallDistanceVec3 b{-1.3041288684536463,4.0231776762652878,-0.69343362405746456};
    const WallDistanceVec3 c{4.195146042365784,12.197819295482756,3.9874020371528283};
    const WallDistanceVec3 p{-3.2828215196543651,2.0481703384841161,0.36056870760956006};
    const double expected=std::min({seg2(p,a,b),seg2(p,a,c),seg2(p,b,c)});
    check(std::abs(point_triangle_distance2(p,a,b,c)-expected)<1e-12,
          "degenerate-triangle geometry must reduce to the closest segment");
    WallSurface degenerate;degenerate.points={a,b,c};degenerate.triangles={{{0,1,2}}};
    const WallDistanceBvh degenerate_bvh(degenerate);
    check(std::abs(std::pow(degenerate_bvh.nearest_distance(p),2)-expected)<1e-12,
          "BVH must use the same degenerate geometry as the independent edge oracle");

    // Pre-patch failure: h_min was 0.005 despite h_x=0.2; no seed was found.
    WallSurface plane;
    plane.points={{0.05,-2,-2},{0.05,2,-2},{0.05,2,2},{0.05,-2,2}};
    plane.triangles={{{0,1,2}},{{0,2,3}}};
    auto stretched=make_wall_distance_grid(7,3,3,{0.1,-0.005,-0.005},
        {0.2,0.005,0.005},[](auto){return false;});
    const auto anisotropic=compute_wall_distance(WallDistanceMethod::EIKONAL,plane,stretched,100);
    check(anisotropic.converged &&
          std::abs(anisotropic.distance[stretched.index(0,1,1)]-0.05)<1e-13 &&
          std::abs(anisotropic.distance[stretched.index(4,1,1)]-0.85)<1e-12,
          "Eikonal failed to seed or propagate on a stretched grid");
    const auto anisotropic_o2=compute_wall_distance(WallDistanceMethod::EIKONAL_SECOND_ORDER,plane,stretched,100);
    check(anisotropic_o2.converged &&
          std::abs(anisotropic_o2.distance[stretched.index(4,1,1)]-0.85)<1e-10,
          "second-order Eikonal is not dispatched on stretched grids");

    // Check invalid geometry before allocation and index overflow.
    bool invalid=false;
    try { (void)make_wall_distance_grid(5,5,5,{0,0,0},{1,0,1},[](auto){return false;}); }
    catch(const std::invalid_argument&) { invalid=true; }
    check(invalid,"zero spacing accepted");
    invalid=false;
    try { (void)make_wall_distance_grid(5,5,5,{0,0,0},
            {std::numeric_limits<double>::quiet_NaN(),1,1},[](auto){return false;}); }
    catch(const std::invalid_argument&) { invalid=true; }
    check(invalid,"NaN spacing accepted");
    invalid=false;
    try { (void)make_wall_distance_grid(std::numeric_limits<std::size_t>::max(),5,5,
            {0,0,0},{1,1,1},[](auto){return false;}); }
    catch(const std::invalid_argument&) { invalid=true; }
    check(invalid,"grid index overflow accepted");

    // A rotated cube and its independent inverse-rotation analytical oracle
    // expose graph anisotropy; tests do not reference any BVH distance oracle.
    constexpr double angle=3.14159265358979323846/6;
    double old_error=1.0;
    for(const std::size_t n:std::array<std::size_t,3>{9,17,25}) {
        const double h=2.0/static_cast<double>(n-1);
        const auto surface=cube(angle);
        const auto g=make_wall_distance_grid(n,n,n,{-0.5,-0.5,-0.5},
            {h,h,h},[&](WallDistanceVec3 q){
                q=unrotate(q,angle);
                return q.x>0 && q.x<1 && q.y>0 && q.y<1 && q.z>0 && q.z<1;
            });
        auto oracle=[&](WallDistanceVec3 q){return cube_analytic(unrotate(q,angle));};
        const auto mesh=compute_wall_distance(WallDistanceMethod::MESH_WAVE,surface,g,100);
        const auto directional=compute_wall_distance(WallDistanceMethod::DIRECTIONAL_MESH_WAVE,surface,g,100);
        const auto ex=compute_wall_distance(WallDistanceMethod::EXACT_GEOMETRIC,surface,g,100);
        const auto bvh=compute_wall_distance(WallDistanceMethod::SEARCH_BASED,surface,g,100);
        const auto o2=compute_wall_distance(WallDistanceMethod::EIKONAL_SECOND_ORDER,surface,g,100);
        const double e2=l2_error(g,o2,oracle);
        check(o2.converged && o2.residual_inf<1e-6 && e2<0.025,
              "second-order Eikonal must converge on each rotated-cube refinement");
        const double em=l2_error(g,mesh,oracle);
        const double ed=l2_error(g,directional,oracle);
        check(directional.stopping_reason=="graph_complete_nonstandard_directional_metric",
              "nonstandard directional metric must be explicitly disclosed");
        check(em<0.025 && em<=old_error+1e-12,
              "wall-point wave has lost mesh-refinement accuracy on rotated cube");
        check(ed<0.06 && ed>em,
              "directional wave must retain an independently identifiable metric");
        check(l2_error(g,ex,oracle)<1e-12 && l2_error(g,bvh,oracle)<1e-12,
              "direct and BVH geometry must agree with analytic rotated cube");
        check(ex.geometrically_exact_on_input_mesh && bvh.geometrically_exact_on_input_mesh &&
              !mesh.geometrically_exact_on_input_mesh && !directional.geometrically_exact_on_input_mesh,
              "geometric provenance flags must distinguish approximate methods");
        old_error=em;
    }

    // On the cube, Poisson is an auxiliary POTENTIAL model. Even a converged
    // true b-A*phi residual cannot certify the reconstructed distance.
    const auto surface=cube();
    const auto g=make_wall_distance_grid(17,17,17,{-0.5,-0.5,-0.5},
        {.125,.125,.125},[](auto q){return q.x>0&&q.x<1&&q.y>0&&q.y<1&&q.z>0&&q.z<1;});
    const auto poisson=compute_wall_distance(WallDistanceMethod::POISSON,surface,g,1000);
    check(poisson.converged && poisson.residual_inf<1e-8 &&
          poisson.poisson_distance_l2_error>0.5 &&
          !poisson.geometrically_exact_on_input_mesh,
          "Poisson PDE convergence was confused with geometric distance accuracy");
    const auto second=compute_wall_distance(WallDistanceMethod::EIKONAL_SECOND_ORDER,surface,g,100);
    const auto first=compute_wall_distance(WallDistanceMethod::EIKONAL,surface,g,100);
    const double second_error=l2_error(g,second,[](auto q){return cube_analytic(q);});
    const double first_error=l2_error(g,first,[](auto q){return cube_analytic(q);});
    check(second.converged && second.residual_inf<1e-6 && second_error<.03 &&
          second_error<first_error*.5,"second-order Eikonal dispatch/accuracy regression");
    check(!second.geometrically_exact_on_input_mesh,
          "PDE Eikonal must not be advertised as geometrically exact");

    // A shifted plane has d=x-a. Diffusive-distance solvers previously forced
    // a homogeneous-Neumann outflow derivative, corrupting this exact affine
    // solution. Now only the artificial OUTER boundary is fixed by BVH.
    WallSurface shifted;
    shifted.points={{.13,-9,-9},{.13,9,-9},{.13,9,9},{.13,-9,9}};
    shifted.triangles={{{0,1,2}},{{0,2,3}}};
    const auto sg=make_wall_distance_grid(17,17,17,{-0.5,-0.5,-0.5},
        {.125,.125,.125},[](auto q){return q.x<.13;});
    const WallDistanceBvh sbvh(shifted);
    const auto planar_oracle=[](WallDistanceVec3 q){return std::abs(q.x-.13);};
    for(auto method:{WallDistanceMethod::HAMILTON_JACOBI,
                     WallDistanceMethod::ADVECTION_DIFFUSION,
                     WallDistanceMethod::HYBRID_POISSON_EIKONAL}) {
        const auto result=compute_wall_distance(method,shifted,sg,1000);
        const auto error=l2_error(sg,result,planar_oracle);
        check(result.converged && error<1e-9,
              "transport method failed the exact planar solution after boundary fix");
    }
    const auto hybrid=compute_wall_distance(WallDistanceMethod::HYBRID_POISSON_EIKONAL,shifted,sg,1000);
    std::vector<unsigned char> empty_fixed(sg.points.size(),0);
    check(std::isinf(wall_distance_pde_residual_inf(
        WallDistanceMethod::HYBRID_POISSON_EIKONAL,hybrid.distance,sg,&empty_fixed)),
        "generic residual must reject hybrid without its Poisson velocity");

    // Previous V2 default silently hid bad distance algorithms behind BVH.
    const auto policy_result=compute_wall_distance_accuracy_checked(
        WallDistanceMethod::POISSON,surface,g,1000);
    check(!policy_result.accuracy_qualified && policy_result.corrected_samples==0 &&
          policy_result.violations>0 && !policy_result.accepted.converged,
          "default accuracy check must NOT replace Poisson with BVH");
    std::cout<<"PASS: geometry degeneracy, anisotropic seeds, exact/BVH, graph metrics, "
                "Poisson truthfulness, Eikonal-2 dispatch, H-J/AD/hybrid BC and fail-closed gate\n";
}
