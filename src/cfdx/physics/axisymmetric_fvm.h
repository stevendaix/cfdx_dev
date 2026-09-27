#pragma once

// CFDX native axisymmetric finite-volume kernel for no-swirl incompressible flow.
//
// Geometry is represented in the meridional (x,r) plane. The revolution factor
// 2*pi*r is part of the operator, not an a-posteriori scaling of a planar solve.
// This module is deliberately independent of the existing 3-D Mesh class so
// that the axisymmetric metric can be verified before wiring it into the
// general solver.

#include <algorithm>
#include "cfdx/physics/forces.h"
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics::axisymmetric {

inline constexpr double pi = 3.141592653589793238462643383279502884;

struct Point { double x=0.0, r=0.0; };
struct Face {
    Point centre;
    double ds=0.0;
    bool axis=false;
};
struct Cell {
    Point centre;
    double area=0.0;       // meridional polygon area
    double volume=0.0;     // physical 3-D control-volume measure
};
struct StructuredMesh {
    std::size_t nx=0, nr=0;
    double xmin=0.0, xmax=0.0, rmax=0.0;
    std::vector<double> x;
    std::vector<double> r;
    std::vector<Cell> cells;
};

inline double revolution_volume(double meridional_area, double mean_radius)
{
    if (!(meridional_area > 0.0) || mean_radius < 0.0 ||
        !std::isfinite(meridional_area) || !std::isfinite(mean_radius))
        throw std::invalid_argument("axisymmetric: invalid volume metric");
    return 2.0*pi*mean_radius*meridional_area;
}

inline double revolution_face_measure(double ds, double radius)
{
    if (!(ds > 0.0) || radius < 0.0 ||
        !std::isfinite(ds) || !std::isfinite(radius))
        throw std::invalid_argument("axisymmetric: invalid face metric");
    return 2.0*pi*radius*ds;
}

inline double cylindrical_divergence(double du_dx, double u_r, double r)
{
    if (!(r > 0.0) || !std::isfinite(r))
        throw std::invalid_argument("axisymmetric: divergence evaluated at r<=0");
    return du_dx + u_r/r;
}

inline double radial_laplacian(double d2ur_dr2, double dur_dr, double ur, double r)
{
    if (!(r > 0.0) || !std::isfinite(r))
        throw std::invalid_argument("axisymmetric: radial Laplacian evaluated at r<=0");
    return d2ur_dr2 + dur_dr/r - ur/(r*r);
}

struct MomentumSource {
    double axial=0.0;
    double radial=0.0;
};

inline MomentumSource cylindrical_momentum_viscous(
    double nu, double d2ux_dx2, double d2ux_dr2, double dux_dr,
    double d2ur_dx2, double d2ur_dr2, double dur_dr, double ur, double r)
{
    if (!(nu >= 0.0) || !std::isfinite(nu))
        throw std::invalid_argument("axisymmetric: invalid viscosity");
    if (!(r > 0.0) || !std::isfinite(r))
        throw std::invalid_argument("axisymmetric: viscous term evaluated at r<=0");

    // No swirl, cylindrical coordinates:
    // ∇²u_x = u_x,xx + u_x,rr + (1/r)u_x,r
    // (∇²u)_r = u_r,xx + u_r,rr + (1/r)u_r,r - u_r/r²
    return {
        nu * (d2ux_dx2 + d2ux_dr2 + dux_dr/r),
        nu * radial_laplacian(d2ur_dr2, dur_dr, ur, r) + nu*d2ur_dx2
    };
}

inline void check_axis_regularity(double ur, double radial_flux, double r)
{
    if (!std::isfinite(ur) || !std::isfinite(radial_flux) || !std::isfinite(r))
        throw std::invalid_argument("axisymmetric: non-finite axis state");
    if (r < 0.0)
        throw std::invalid_argument("axisymmetric: negative radius");
    if (r == 0.0 && (std::abs(ur) > 1e-12 || std::abs(radial_flux) > 1e-12))
        throw std::runtime_error("axisymmetric: axis regularity violated");
}

struct BoundaryCondition {
    enum class Type { AXIS, WALL, INLET, OUTLET, FARFIELD };
    Type type;
    double ux=0.0;
    double ur=0.0;
    double pressure=0.0;
};

struct PressureFace {
    std::size_t owner=0;
    std::size_t neighbour=0;
    double area=0.0;
    double distance=0.0;
    double radius=0.0;
};

inline double pressure_correction_coefficient(
    const PressureFace& f, double rho, double dAU)
{
    if (!(rho > 0.0) || !(dAU > 0.0) || !(f.area > 0.0) || !(f.distance > 0.0) ||
        !(f.radius >= 0.0))
        throw std::invalid_argument("axisymmetric: invalid pressure coefficient");
    return rho * dAU * f.area / f.distance;
}

// Force extraction is provided by cfdx/physics/forces.h. This kernel only
// verifies the axisymmetric metric and momentum operators; it intentionally
// does not maintain a second force-integration implementation.

struct VMFL036Definition {
    double diameter=1.0;
    double rho=1.0;
    double velocity=1.0;
    double viscosity=0.01;
    double outer_radius=50.0;
    double re() const { return rho*velocity*diameter/viscosity; }
    double reference_area() const { return pi*diameter*diameter/4.0; }
    double reference_cd() const { return 1.0895; }
};

inline VMFL036Definition make_vmfl036_re100()
{
    return {};
}

inline VMFL036Definition make_vmfl036_fluent_exact()
{
    VMFL036Definition c;
    c.viscosity=0.02;
    return c;
}

inline void validate_definition(const VMFL036Definition& c)
{
    if (!(c.diameter>0.0) || !(c.rho>0.0) || !(c.velocity>0.0) ||
        !(c.viscosity>0.0) || !(c.outer_radius>=50.0*c.diameter))
        throw std::invalid_argument("axisymmetric VMFL036: invalid definition");
}

} // namespace cfdx::physics::axisymmetric
