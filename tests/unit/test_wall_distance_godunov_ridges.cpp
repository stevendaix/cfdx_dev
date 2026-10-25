#include "cfdx/physics/wall_distance.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
int main() {
    WallSurface wall;
    for(double x:{0.0,2.0}) {
        const std::size_t b=wall.points.size();
        wall.points.insert(wall.points.end(),{{x,-4,-4},{x,4,-4},{x,4,4},{x,-4,4}});
        wall.triangles.push_back({{b,b+1,b+2}});
        wall.triangles.push_back({{b,b+2,b+3}});
    }
    const auto grid=make_wall_distance_grid(9,3,3,{0,-1,-1},{.25,1,1},
        [](const WallDistanceVec3& p){return p.x<=0.0 || p.x>=2.0;});
    const WallDistanceBvh bvh(wall);
    const WallDistanceCutFaceOffsetCache cache(bvh,grid);
    const auto expected=exact_reference(wall,grid);
    const auto result=compute_wall_distance(WallDistanceMethod::EIKONAL,wall,grid,50);
    const std::size_t ridge=grid.index(4,1,1);
    const auto require=[](bool ok,const char* message){if(!ok)throw std::runtime_error(message);};
    require(std::abs(result.distance[ridge]-1.0)<1e-13,"channel ridge distance is not exact");
    for(const std::size_t i:{2u,4u,6u}) {
        const std::size_t id=grid.index(i,1,1);
        require(std::abs(godunov_gradient_at(expected,grid,id)-1.0)<1e-13,
                "Godunov gradient double-counted opposite directions");
        require(std::abs(godunov_gradient_at(bvh,expected,grid,id)-1.0)<1e-13,
                "BVH Godunov gradient double-counted opposite directions");
        require(std::abs(godunov_gradient_at(bvh,expected,grid,id,&cache)-1.0)<1e-13,
                "cached Godunov gradient is inconsistent");
    }
    const auto metrics=compare_wall_distance(grid,expected,result.distance,.5);
    require(metrics.l2_relative<1e-13,"channel distance differs from independent two-wall oracle");
    require(result.converged && result.residual_inf<1e-13,
            "exact channel distance is incorrectly reported as unconverged");
    std::cout << "Godunov ridge and monotone branches: PASS, L2=" << metrics.l2_relative
              << " residual=" << result.residual_inf << '\n';
}
