#include "cfdx/physics/wall_face_reconstruction_quadratic.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using cfdx::core::Vec3;
using cfdx::physics::WallQuadraticSample;
using cfdx::physics::reconstruct_wall_face_quadratic;

static double pressure(Vec3 x) {
    return 2.0 + 0.7*x.x - 0.4*x.y + 0.3*x.z +
        0.8*x.x*x.x - 0.2*x.x*x.y + 0.1*x.x*x.z +
        0.5*x.y*x.y + 0.15*x.y*x.z - 0.35*x.z*x.z;
}
static Vec3 velocity(Vec3 x) {
    return {x.x*x.x + 2*x.x*x.y + 0.3*x.z*x.z,
            0.5*x.y*x.y - x.x*x.z + 0.2*x.y*x.z,
            0.4*x.z*x.z + x.x*x.y - 0.1*x.x*x.z};
}
static std::array<Vec3,3> exact_gradient(Vec3 x) {
    return {{{2*x.x+2*x.y, 2*x.x, 0.6*x.z},
             {-x.z, x.y+0.2*x.z, -x.x+0.2*x.y},
             {x.y-0.1*x.z, x.x, 0.8*x.z-0.1*x.x}}};
}
static Vec3 cubic_velocity(Vec3 x) {
    return {x.z + 0.8*x.x*x.z*x.z,
            2.0*x.z + 0.4*x.y*x.z*x.z,
            -0.5*x.z + 0.5*x.z*x.z*x.z};
}
int main() {
    const Vec3 xf{0.17,-0.11,0.23}, sf{0.8,-0.3,0.5};
    std::vector<WallQuadraticSample> samples;
    // Non-coplanar, deliberately irregular 3-D stencil.
    const Vec3 offsets[] = {
        {-1.,0.,0.},{1.,0.,0.},{0.,-1.,0.},{0.,1.,0.},
        {0.,0.,-1.},{0.,0.,1.},{1.,1.,0.},{-1.,0.,1.},
        {0.,-1.,1.},{1.,0.,-1.},{-1.,-1.,0.},{1.,-1.,1.}
    };
    for (const auto& d : offsets) {
        const Vec3 x = xf + d*0.37;
        samples.push_back({x,pressure(x),velocity(x)-velocity(xf)});
    }
    // Enforce an actual stationary no-slip wall value at the face.
    samples.push_back({xf,pressure(xf),Vec3{0.0,0.0,0.0}});
    const auto result = reconstruct_wall_face_quadratic(samples,xf,sf,0.7);
    const auto exact = exact_gradient(xf);
    for (std::size_t i=0;i<3;++i)
        if ((result.velocity_gradient[i]-exact[i]).mag() > 2e-11)
            throw std::runtime_error("quadratic wall gradient is not exact");
    if (std::abs(result.pressure-pressure(xf)) > 2e-11)
        throw std::runtime_error("quadratic wall pressure is not exact");
    const auto n = sf.normalized();
    const auto& g = exact;
    const Vec3 tau_n{
        0.7*(2*g[0].x*n.x+(g[0].y+g[1].x)*n.y+(g[0].z+g[2].x)*n.z),
        0.7*((g[1].x+g[0].y)*n.x+2*g[1].y*n.y+(g[1].z+g[2].y)*n.z),
        0.7*((g[2].x+g[0].z)*n.x+(g[2].y+g[1].z)*n.y+2*g[2].z*n.z)};
    const Vec3 expected_force = n*pressure(xf)-tau_n;
    if ((result.force_per_area-expected_force).mag() > 3e-11)
        throw std::runtime_error("quadratic wall traction is not exact");
    bool rejected=false;
    try {
        std::vector<WallQuadraticSample> planar;
        for (int i=0;i<13;++i) {
            const Vec3 x{double(i%4),double(i/4),0.0};
            planar.push_back({x,pressure(x),velocity(x)});
        }
        (void)reconstruct_wall_face_quadratic(planar,{0.,0.,0.},{1.,0.,0.},1.);
    } catch (const std::invalid_argument&) { rejected=true; }
    if (!rejected) throw std::runtime_error("rank-deficient stencil was not rejected");
    // Cubic manufactured fields expose the asymptotic error of a quadratic
    // fit. The one-sided cloud is representative of a wall stencil; the
    // exact wall value is zero and is included explicitly at the face.
    std::array<double,3> errors{};
    const double h_values[] = {0.4, 0.2, 0.1};
    const std::array<Vec3,3> exact_cubic_gradient{{
        Vec3{0.0,0.0,1.0}, Vec3{0.0,0.0,2.0}, Vec3{0.0,0.0,-0.5}}};
    for (std::size_t level=0; level<3; ++level) {
        const double h = h_values[level];
        std::vector<WallQuadraticSample> cloud;
        for (int i=-2; i<=2; ++i)
            for (int j=-2; j<=2; ++j)
                for (int k=0; k<5; ++k) {
                    const Vec3 x{h*static_cast<double>(i),
                                 h*static_cast<double>(j),
                                 h*(static_cast<double>(k)+0.5)};
                    cloud.push_back({x, pressure(x), cubic_velocity(x)});
                }
        cloud.push_back({Vec3{0.0,0.0,0.0}, pressure(Vec3{0.0,0.0,0.0}),
                         Vec3{0.0,0.0,0.0}});
        const auto fit = reconstruct_wall_face_quadratic(
            cloud, Vec3{0.0,0.0,0.0}, Vec3{0.0,0.0,1.0}, 0.7);
        for (std::size_t component=0; component<3; ++component)
            errors[level] += (fit.velocity_gradient[component] -
                              exact_cubic_gradient[component]).mag2();
        errors[level] = std::sqrt(errors[level]);
    }
    const double order_01 = std::log(errors[0]/errors[1]) / std::log(2.0);
    const double order_12 = std::log(errors[1]/errors[2]) / std::log(2.0);
    if (!(std::isfinite(order_01) && std::isfinite(order_12) &&
          order_01 > 1.7 && order_12 > 1.7))
        throw std::runtime_error("quadratic wall reconstruction did not show second-order convergence");
    std::cout << "WALL_FACE_QUADRATIC_3D_MMS PASS polynomial_reproduction=2"
              << " no_slip=PASS observed_orders=" << order_01 << "," << order_12 << '\n';
}
