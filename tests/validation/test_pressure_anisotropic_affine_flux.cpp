// Phase C independent manufactured-solution certificate for the production
// anisotropic pressure-face flux decomposition. This is NOT a full coupled
// Jacobian qualification on skewed meshes.
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("affine_pressure_response_rotated_stretched_faces", [] {
        struct Case { Vec3 face; Vec3 center_distance; Vec3 response; Vec3 gradient; };
        const std::array<Case, 5> cases{{
            {{1.0,0.0,0.0}, {0.01,0.0,0.0}, {1e-6,1.0,1.0}, {2.0,-3.0,0.5}},
            {{1.0,0.5,0.0}, {0.01,0.0008,0.0}, {1e-6,1e-2,1.0}, {2.0,-3.0,0.5}},
            {{0.6,1.0,0.1}, {0.002,0.02,0.001}, {1e-4,1000.0,2.0}, {-6.0,2.0,1.0}},
            {{1.0,0.2,-0.15}, {1e-3,2e-4,-1e-4}, {1e-6,1e2,1e4}, {0.5,-4.0,3.0}},
            {{0.1,0.0,1.0}, {1e-4,5e-5,1e-3}, {1e2,1e-6,1.0}, {8.0,1.0,-2.0}}
        }};
        for (const auto& c : cases) {
            const double jump = c.gradient.dot(c.center_distance);
            const double coefficient = pressure_face_orthogonal_coefficient(
                c.face,c.center_distance,c.response);
            const double explicit_part = pressure_face_nonorthogonal_flux(
                c.face,c.center_distance,c.response,c.gradient);
            const double recovered = coefficient*jump + explicit_part;
            const double exact =
                c.face.x*c.response.x*c.gradient.x +
                c.face.y*c.response.y*c.gradient.y +
                c.face.z*c.response.z*c.gradient.z;
            EXPECT_NEAR(recovered,exact,1e-10*std::max(1.0,std::abs(exact)));
            EXPECT_TRUE(std::isfinite(coefficient));
            EXPECT_TRUE(coefficient>0.0);
        }
    });
    run_case("owner_neighbour_face_signs_conserve_internal_exchange", [] {
        const Vec3 area{0.9,0.6,0.2}, d{0.002,0.001,0.0002};
        const Vec3 response{1e-4,100.0,3.0}, g{1.0,-2.0,0.5};
        const double orth=pressure_face_orthogonal_coefficient(area,d,response);
        const double deferred=pressure_face_nonorthogonal_flux(area,d,response,g);
        const double owner=orth*g.dot(d)+deferred;
        const double neighbour=-owner;
        EXPECT_NEAR(owner+neighbour,0.0,1e-14);
        EXPECT_NEAR(owner,response.x*g.x*area.x+response.y*g.y*area.y+
                     response.z*g.z*area.z,1e-12);
    });
    run_case("zero_distance_and_inconsistent_normal_fail_closed", [] {
        EXPECT_THROW(pressure_face_orthogonal_coefficient(
            {1,0,0},{0,0,0},{1,1,1}),std::runtime_error);
        EXPECT_THROW(pressure_face_orthogonal_coefficient(
            {-1,0,0},{1,0,0},{1,1,1}),std::runtime_error);
        EXPECT_THROW(pressure_face_nonorthogonal_flux(
            {1,0,0},{0,0,0},{1,1,1},{1,0,0}),std::runtime_error);
    });
    return run_all();
}