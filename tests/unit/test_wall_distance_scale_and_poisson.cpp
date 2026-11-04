// Regression for two numerical defects isolated by independent geometry and
// floating-point counterexamples. No BVH fallback substitutes a PDE result.
#include "cfdx/physics/wall_distance.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace cfdx::physics;

static void demand(bool ok,const char* message) {
    if (!ok) throw std::runtime_error(message);
}

static_assert(!std::is_constructible<WallDistanceBvh,WallSurface&&>::value,
              "BVH cannot retain a temporary WallSurface");
static_assert(!std::is_constructible<WallDistanceBvh,const WallSurface&&>::value,
              "BVH cannot retain a const temporary WallSurface");

int main() {
    // A ray crosses the INTERIOR of a triangle whose size changes over ten
    // decades. The expected answer is analytic; nearest_distance() is only
    // checked to ensure the surface itself is geometrically well formed.
    for (double scale : {1.0,1e-4,1e-6,1e-8,1e-10}) {
        WallSurface wall;
        wall.points={{0,0,0},{scale,0,0},{0,scale,0}};
        wall.triangles.push_back({{{0,1,2}}});
        const WallDistanceBvh bvh(wall);
        const WallDistanceVec3 origin{scale/4,scale/4,-scale};
        const double ray=bvh.ray_distance(origin,{0,0,1},2*scale);
        demand(std::isfinite(ray) && std::abs(ray/scale-1.0)<1e-10,
               "BVH ray misses a small nondegenerate triangle");
        demand(std::abs(bvh.nearest_distance(origin)/scale-1.0)<1e-10,
               "nearest-point reference is wrong");
    }

    // The positive solution d=sqrt(g*g+2*phi)-g loses precision for small phi.
    demand(std::abs(poisson_distance_from_potential(1e-20,1.0)/1e-20-1.0)<1e-12,
           "Poisson reconstructed distance cancelled to zero");
    demand(std::abs(poisson_distance_from_potential(0.5,0.0)-1.0)<1e-14,
           "Poisson reconstruction at a medial-axis extremum is wrong");
    demand(poisson_distance_from_potential(0.0,1.0)==0.0,
           "zero wall potential did not produce zero distance");

    // Exercise the public reconstruction, not only its scalar helper.
    WallSurface wall;
    wall.points={{0,-10,-10},{0,10,-10},{0,0,10}};
    wall.triangles.push_back({{{0,1,2}}});
    const WallDistanceBvh bvh(wall);
    auto grid=make_wall_distance_grid(3,3,3,{-1,-1,-1},{1,1,1},
                                     [](const WallDistanceVec3&){return false;});
    std::vector<double> phi(grid.points.size());
    for (std::size_t i=0;i<phi.size();++i) phi[i]=grid.points[i].x;
    const auto center=grid.index(1,1,1);
    phi[center]=1e-20;
    const double d=poisson_reconstructed_distance(bvh,grid,phi,center);
    demand(std::abs(d/1e-20-1.0)<1e-12,
           "full Poisson gradient/reconstruction did not preserve a tiny distance");

    std::cout << "PASS: ray scale invariance, Poisson cancellation, medial-axis reconstruction\n";
}
