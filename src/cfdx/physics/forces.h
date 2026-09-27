#pragma once

// Surface force and moment integration utilities.
//
// The evaluator is deliberately solver-agnostic: it consumes surface pressure,
// viscous traction/stress and geometry, and returns dimensional forces, moments,
// and non-dimensional coefficients. It is intended for walls, bodies, porous
// interfaces, axisymmetric surfaces and future moving/deforming meshes.
//
// Sign convention:
//   n points from the body into the fluid.
//   sigma = -p I + tau.
//   fluid force on the body = sigma . n.
// This convention makes a positive streamwise drag correspond to the force
// exerted by the fluid on a stationary body in the +drag direction.

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics::forces {

struct Vec3 {
    double x{0.0}, y{0.0}, z{0.0};

    Vec3& operator+=(const Vec3& o) noexcept {
        x += o.x; y += o.y; z += o.z; return *this;
    }
};

inline Vec3 operator+(Vec3 a, const Vec3& b) noexcept { return a += b; }
inline Vec3 operator-(Vec3 a, const Vec3& b) noexcept {
    return {a.x-b.x, a.y-b.y, a.z-b.z};
}
inline Vec3 operator*(double a, const Vec3& v) noexcept {
    return {a*v.x, a*v.y, a*v.z};
}
inline Vec3 operator*(const Vec3& v, double a) noexcept { return a*v; }

inline double dot(const Vec3& a, const Vec3& b) noexcept {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
inline Vec3 cross(const Vec3& a, const Vec3& b) noexcept {
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}
inline double norm(const Vec3& v) noexcept { return std::sqrt(dot(v,v)); }

inline Vec3 normalized(Vec3 v) {
    const double n=norm(v);
    if (!(n>0.0) || !std::isfinite(n))
        throw std::invalid_argument("forces: zero or non-finite direction");
    return (1.0/n)*v;
}

struct SurfaceSample {
    Vec3 position{};
    Vec3 normal{};          // unit normal: body -> fluid
    double area{0.0};       // physical surface measure
    double pressure{0.0};   // dimensional pressure [Pa]
    Vec3 viscous_traction{}; // tau . n [Pa]
};

struct Reference {
    double rho_inf{1.0};
    double u_inf{1.0};
    double area{1.0};
    double length{1.0};
    double pressure_ref{0.0};
    Vec3 drag_direction{1.0,0.0,0.0};
    Vec3 lift_direction{0.0,1.0,0.0};
    Vec3 side_direction{0.0,0.0,1.0};
    Vec3 moment_center{};
};

struct Result {
    Vec3 pressure_force{};
    Vec3 viscous_force{};
    Vec3 total_force{};
    Vec3 pressure_moment{};
    Vec3 viscous_moment{};
    Vec3 total_moment{};
    double drag_pressure{0.0};
    double drag_viscous{0.0};
    double drag{0.0};
    double lift_pressure{0.0};
    double lift_viscous{0.0};
    double lift{0.0};
    double side_pressure{0.0};
    double side_viscous{0.0};
    double side{0.0};
    double cd_pressure{0.0};
    double cd_viscous{0.0};
    double cd{0.0};
    double cl_pressure{0.0};
    double cl_viscous{0.0};
    double cl{0.0};
    double cs_pressure{0.0};
    double cs_viscous{0.0};
    double cs{0.0};
    double cmx{0.0};
    double cmy{0.0};
    double cmz{0.0};
};

inline void validate_reference(const Reference& r) {
    if (!std::isfinite(r.rho_inf) || !std::isfinite(r.u_inf) || !std::isfinite(r.area) ||
        !std::isfinite(r.length) ||
        !(r.rho_inf>0.0) || !(r.u_inf>0.0) || !(r.area>0.0) || !(r.length>0.0))
        throw std::invalid_argument("forces: invalid reference density, velocity, area or length");
    if (!std::isfinite(r.pressure_ref))
        throw std::invalid_argument("forces: non-finite reference pressure");
    (void)normalized(r.drag_direction);
    (void)normalized(r.lift_direction);
    (void)normalized(r.side_direction);
}

inline Result integrate(const std::vector<SurfaceSample>& samples,
                        const Reference& ref)
{
    validate_reference(ref);
    const Vec3 eD=normalized(ref.drag_direction);
    const Vec3 eL=normalized(ref.lift_direction);
    const Vec3 eS=normalized(ref.side_direction);
    const double q=0.5*ref.rho_inf*ref.u_inf*ref.u_inf;
    if (!(q>0.0) || !std::isfinite(q))
        throw std::invalid_argument("forces: invalid dynamic pressure");

    Result out{};
    for(std::size_t i=0;i<samples.size();++i){
        const auto& s=samples[i];
        const double nn=norm(s.normal);
        if (!(s.area>=0.0) || !std::isfinite(s.area) || !(std::abs(nn-1.0)<1e-10) ||
            !std::isfinite(s.pressure) || !std::isfinite(nn) ||
            !std::isfinite(s.position.x) || !std::isfinite(s.position.y) ||
            !std::isfinite(s.position.z) || !std::isfinite(s.viscous_traction.x) ||
            !std::isfinite(s.viscous_traction.y) || !std::isfinite(s.viscous_traction.z))
            throw std::invalid_argument("forces: invalid surface sample at index "+std::to_string(i));

        const Vec3 fp=(-(s.pressure-ref.pressure_ref))*s.normal*s.area;
        const Vec3 fv=s.viscous_traction*s.area;
        const Vec3 r=s.position-ref.moment_center;
        out.pressure_force += fp;
        out.viscous_force += fv;
        out.pressure_moment += cross(r,fp);
        out.viscous_moment += cross(r,fv);
    }

    out.total_force=out.pressure_force+out.viscous_force;
    out.total_moment=out.pressure_moment+out.viscous_moment;
    out.drag_pressure=dot(out.pressure_force,eD);
    out.drag_viscous=dot(out.viscous_force,eD);
    out.drag=dot(out.total_force,eD);
    out.lift_pressure=dot(out.pressure_force,eL);
    out.lift_viscous=dot(out.viscous_force,eL);
    out.lift=dot(out.total_force,eL);
    out.side_pressure=dot(out.pressure_force,eS);
    out.side_viscous=dot(out.viscous_force,eS);
    out.side=dot(out.total_force,eS);
    const double qA=q*ref.area;
    out.cd_pressure=out.drag_pressure/qA;
    out.cd_viscous=out.drag_viscous/qA;
    out.cd=out.drag/qA;
    out.cl_pressure=out.lift_pressure/qA;
    out.cl_viscous=out.lift_viscous/qA;
    out.cl=out.lift/qA;
    out.cs_pressure=out.side_pressure/qA;
    out.cs_viscous=out.side_viscous/qA;
    out.cs=out.side/qA;
    const double qAL=q*ref.area*ref.length;
    out.cmx=dot(out.total_moment,{1.0,0.0,0.0})/qAL;
    out.cmy=dot(out.total_moment,{0.0,1.0,0.0})/qAL;
    out.cmz=dot(out.total_moment,{0.0,0.0,1.0})/qAL;
    return out;
}

// Native axisymmetric surface sample. The meridional normal points from the
// body into the fluid. Revolution measure is 2*pi*r*ds.
struct AxisymmetricSample {
    double x{0.0};
    double r{0.0};
    double nx{1.0};
    double nr{0.0};
    double ds{0.0};
    double pressure{0.0};
    double tau_xx{0.0};
    double tau_xr{0.0};
    double tau_rr{0.0};
};

struct AxisymmetricResult {
    double pressure_force{0.0};
    double viscous_force{0.0};
    double total_force{0.0};
    double cd_pressure{0.0};
    double cd_viscous{0.0};
    double cd_total{0.0};
};

inline AxisymmetricResult integrate_axisymmetric(
    const std::vector<AxisymmetricSample>& samples,
    double rho_inf, double u_inf, double reference_area, double pressure_ref = 0.0)
{
    if (!std::isfinite(rho_inf) || !std::isfinite(u_inf) || !std::isfinite(reference_area) ||
        !(rho_inf>0.0) || !(u_inf>0.0) || !(reference_area>0.0))
        throw std::invalid_argument("forces: invalid axisymmetric reference");
    const double q=0.5*rho_inf*u_inf*u_inf;
    if (!std::isfinite(q) || !(q>0.0))
        throw std::invalid_argument("forces: invalid axisymmetric dynamic pressure");
    if(!std::isfinite(pressure_ref))
        throw std::invalid_argument("forces: invalid axisymmetric reference pressure");
    AxisymmetricResult out{};
    for(std::size_t i=0;i<samples.size();++i){
        const auto& s=samples[i];
        const double nn=std::hypot(s.nx,s.nr);
        if (!(s.r>=0.0) || !(s.ds>=0.0) || !(std::abs(nn-1.0)<1e-10) ||
            !std::isfinite(s.x) || !std::isfinite(s.r) || !std::isfinite(s.ds) ||
            !std::isfinite(s.pressure) || !std::isfinite(s.tau_xx) ||
            !std::isfinite(s.tau_xr) || !std::isfinite(s.tau_rr))
            throw std::invalid_argument("forces: invalid axisymmetric surface sample at index "+std::to_string(i));
        const double dS=2.0*3.141592653589793238462643383279502884*s.r*s.ds;
        if(!std::isfinite(dS) || dS<0.0)
            throw std::invalid_argument("forces: invalid axisymmetric surface measure at index "+std::to_string(i));
        // sigma . n, axial component. sigma = -p I + tau.
        const double fp=-(s.pressure-pressure_ref)*s.nx*dS;
        const double fv=(s.tau_xx*s.nx+s.tau_xr*s.nr)*dS;
        out.pressure_force += fp;
        out.viscous_force += fv;
    }
    out.total_force=out.pressure_force+out.viscous_force;
    const double qA=q*reference_area;
    out.cd_pressure=out.pressure_force/qA;
    out.cd_viscous=out.viscous_force/qA;
    out.cd_total=out.total_force/qA;
    return out;
}

} // namespace cfdx::physics::forces
