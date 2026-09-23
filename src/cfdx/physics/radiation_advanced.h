#pragma once

#include "cfdx/physics/radiation.h"
#include "cfdx/physics/radiation_solver.h"
#include "cfdx/physics/finite_volume_transport.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

// -----------------------------------------------------------------------------
// Radiation material properties
// -----------------------------------------------------------------------------

struct RadiationOpticalProperties {
    double absorption = 0.0;   // 1/m
    double scattering = 0.0;   // 1/m
    double refractive_index = 1.0;
    double emissivity = 1.0;
    double anisotropy = 0.0;   // Henyey-Greenstein g, [-1,1]

    void validate() const {
        if (!std::isfinite(absorption) || !std::isfinite(scattering) ||
            !std::isfinite(refractive_index) || !std::isfinite(emissivity) ||
            !std::isfinite(anisotropy) ||
            absorption < 0.0 || scattering < 0.0 ||
            refractive_index <= 0.0 || emissivity < 0.0 || emissivity > 1.0 ||
            anisotropy < -1.0 || anisotropy > 1.0) {
            throw std::invalid_argument("invalid radiation optical properties");
        }
    }

    double extinction() const { return absorption + scattering; }

    double optical_thickness(double length) const {
        validate();
        if (!std::isfinite(length) || length < 0.0)
            throw std::invalid_argument("radiation path length must be finite and non-negative");
        return extinction() * length;
    }
};

struct RadiationOpticalPropertyField {
    std::vector<RadiationOpticalProperties> cells;

    void validate() const {
        for (const auto& p : cells) p.validate();
    }

    RadiationOpticalProperties& operator[](std::size_t i) { return cells.at(i); }
    const RadiationOpticalProperties& operator[](std::size_t i) const { return cells.at(i); }
};

// -----------------------------------------------------------------------------
// Spectral / non-gray infrastructure
// -----------------------------------------------------------------------------

struct RadiationBand {
    double wavelength_min = 0.0; // m
    double wavelength_max = 0.0; // m
    double weight = 1.0;
    RadiationOpticalProperties properties;

    void validate() const {
        if (!std::isfinite(wavelength_min) || !std::isfinite(wavelength_max) ||
            !std::isfinite(weight) || wavelength_min < 0.0 ||
            wavelength_max <= wavelength_min || weight < 0.0) {
            throw std::invalid_argument("invalid radiation spectral band");
        }
        properties.validate();
    }
};

inline void validate_radiation_bands(const std::vector<RadiationBand>& bands)
{
    if (bands.empty())
        throw std::invalid_argument("at least one radiation spectral band is required");
    double weight_sum = 0.0;
    double previous = 0.0;
    for (const auto& b : bands) {
        b.validate();
        if (b.wavelength_min < previous)
            throw std::invalid_argument("radiation spectral bands must be ordered");
        previous = b.wavelength_max;
        weight_sum += b.weight;
    }
    if (!(weight_sum > 0.0) || !std::isfinite(weight_sum))
        throw std::invalid_argument("radiation spectral band weights must have positive sum");
}

inline double weighted_band_absorption(
    const std::vector<RadiationBand>& bands,
    double temperature)
{
    validate_radiation_bands(bands);
    if (!std::isfinite(temperature) || temperature <= 0.0)
        throw std::invalid_argument("invalid temperature for spectral radiation");
    double numerator = 0.0, denominator = 0.0;
    for (const auto& b : bands) {
        numerator += b.weight * b.properties.absorption;
        denominator += b.weight;
    }
    return numerator / denominator;
}

// -----------------------------------------------------------------------------
// Diffuse-gray DOM wall operator
// -----------------------------------------------------------------------------

inline double diffuse_gray_wall_intensity(
    double emissivity,
    double temperature,
    double incoming_irradiation)
{
    return gray_diffuse_wall_intensity(emissivity, temperature,
                                       incoming_irradiation);
}

inline double hemispherical_irradiation(
    const std::vector<DiscreteDirection>& directions,
    const std::vector<double>& intensities,
    const std::array<double,3>& normal)
{
    if (directions.size() != intensities.size())
        throw std::invalid_argument("DOM direction/intensity size mismatch");
    const double nn = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1] +
                                normal[2]*normal[2]);
    if (!std::isfinite(nn) || nn <= 0.0)
        throw std::invalid_argument("wall normal must be finite and non-zero");
    double G = 0.0;
    for (std::size_t m = 0; m < directions.size(); ++m) {
        const auto& d = directions[m];
        const double mu = (d.dx*normal[0] + d.dy*normal[1] + d.dz*normal[2]) / nn;
        if (mu < 0.0) {
            if (!std::isfinite(intensities[m]) || intensities[m] < 0.0)
                throw std::invalid_argument("incoming DOM intensity must be finite and non-negative");
            G += directions[m].weight * (-mu) * intensities[m];
        }
    }
    return G;
}

inline void apply_diffuse_gray_wall(
    const std::vector<DiscreteDirection>& directions,
    const std::vector<double>& incoming_intensities,
    const std::array<double,3>& outward_normal,
    double emissivity,
    double temperature,
    std::vector<double>& boundary_intensities)
{
    validate_discrete_directions(directions);
    if (incoming_intensities.size() != directions.size() ||
        boundary_intensities.size() != directions.size())
        throw std::invalid_argument("DOM wall direction/intensity size mismatch");
    const double G = hemispherical_irradiation(
        directions, incoming_intensities, outward_normal);
    const double Iw = gray_diffuse_wall_intensity(emissivity, temperature, G);
    const double nn = std::sqrt(outward_normal[0]*outward_normal[0] +
                                outward_normal[1]*outward_normal[1] +
                                outward_normal[2]*outward_normal[2]);
    for (std::size_t m = 0; m < directions.size(); ++m) {
        const double mu = (directions[m].dx*outward_normal[0] +
                           directions[m].dy*outward_normal[1] +
                           directions[m].dz*outward_normal[2]) / nn;
        if (mu > 0.0)
            boundary_intensities[m] = Iw;
    }
}

// -----------------------------------------------------------------------------
// View-factor matrix construction for a supplied geometric kernel.
//
// The geometry backend deliberately operates on a generic differential-patch
// representation so it can be fed by any mesh/face accelerator. The kernel
// uses the centroid differential-area approximation and is conservative only
// after closure. It is therefore a deterministic geometry estimator, not a
// replacement for a visibility/ray-tracing backend.
// -----------------------------------------------------------------------------

struct ViewFactorPatch {
    std::array<double,3> center{0.0,0.0,0.0};
    std::array<double,3> normal{0.0,0.0,1.0};
    double area = 0.0;
};

inline double patch_pair_view_factor(
    const ViewFactorPatch& source,
    const ViewFactorPatch& target)
{
    const double rx = target.center[0]-source.center[0];
    const double ry = target.center[1]-source.center[1];
    const double rz = target.center[2]-source.center[2];
    const double r2 = rx*rx+ry*ry+rz*rz;
    if (!std::isfinite(r2) || r2 <= 0.0)
        return 0.0;
    const double r = std::sqrt(r2);
    const double ns = std::sqrt(source.normal[0]*source.normal[0] +
                                 source.normal[1]*source.normal[1] +
                                 source.normal[2]*source.normal[2]);
    const double nt = std::sqrt(target.normal[0]*target.normal[0] +
                                 target.normal[1]*target.normal[1] +
                                 target.normal[2]*target.normal[2]);
    if (!(source.area > 0.0) || !(target.area > 0.0) ||
        !(ns > 0.0) || !(nt > 0.0))
        throw std::invalid_argument("invalid view-factor patch");
    const double cs = (source.normal[0]*rx + source.normal[1]*ry +
                       source.normal[2]*rz)/(ns*r);
    const double ct = -(target.normal[0]*rx + target.normal[1]*ry +
                        target.normal[2]*rz)/(nt*r);
    if (cs <= 0.0 || ct <= 0.0)
        return 0.0;
    return std::max(0.0, std::min(1.0, target.area * cs * ct /
                                        (M_PI * r2)));
}

inline std::vector<double> estimate_view_factor_matrix(
    const std::vector<ViewFactorPatch>& patches)
{
    if (patches.empty())
        throw std::invalid_argument("view-factor patch set is empty");
    const std::size_t n = patches.size();
    std::vector<double> F(n*n, 0.0);
    for (std::size_t i=0; i<n; ++i)
        for (std::size_t j=0; j<n; ++j)
            if (i != j)
                F[i*n+j] = patch_pair_view_factor(patches[i], patches[j]);

    // Enforce enclosure closure while preserving the physically required
    // reciprocity as far as the centroid approximation permits.
    for (std::size_t i=0; i<n; ++i) {
        double sum = 0.0;
        for (std::size_t j=0; j<n; ++j) sum += F[i*n+j];
        if (sum > 1.0) {
            for (std::size_t j=0; j<n; ++j) F[i*n+j] /= sum;
        }
    }
    return F;
}

// -----------------------------------------------------------------------------
// Radiation diagnostics / conservation
// -----------------------------------------------------------------------------

struct RadiationBalance {
    double emitted = 0.0;
    double absorbed = 0.0;
    double net = 0.0;
    double relative_error = 0.0;
};

inline RadiationBalance radiation_balance(double emitted, double absorbed)
{
    if (!std::isfinite(emitted) || !std::isfinite(absorbed) ||
        emitted < 0.0 || absorbed < 0.0)
        throw std::invalid_argument("invalid radiation balance");
    RadiationBalance b;
    b.emitted = emitted;
    b.absorbed = absorbed;
    b.net = emitted - absorbed;
    b.relative_error = std::abs(b.net) /
        std::max(1.0, std::max(emitted, absorbed));
    return b;
}


// -----------------------------------------------------------------------------
// Band-wise gray/fvDOM radiation. Each band is transported independently and
// the Planck-band weights are used to reconstruct the total gray-equivalent
// irradiation and source. This is the extensible non-gray foundation.
// -----------------------------------------------------------------------------

struct SpectralRadiationSolveResult {
    bool converged = false;
    std::size_t iterations = 0;
    std::vector<RadiationSolveResult> band_results;
};

inline SpectralRadiationSolveResult solve_spectral_dom(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& irradiation,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& radiation_source,
    const std::vector<DiscreteDirection>& directions,
    const std::vector<RadiationBand>& bands,
    RadiationTransportControls controls = {},
    const ScalarBoundaryConditions& wall_intensity_bcs = {})
{
    validate_radiation_bands(bands);
    if(temperature.size()!=mesh.n_cells() ||
       irradiation.size()!=mesh.n_cells() ||
       radiation_source.size()!=mesh.n_cells())
        throw std::invalid_argument("spectral radiation field size mismatch");

    double wsum=0.0;
    for(const auto& b:bands) wsum+=b.weight;
    irradiation.fill(0.0);
    radiation_source.fill(0.0);

    SpectralRadiationSolveResult result;
    result.converged=true;
    for(const auto& band:bands) {
        RadiationOpticalPropertyField p;
        p.cells.resize(mesh.n_cells(),band.properties);
        cfdx::core::Field<double,cfdx::core::Location::CELL> G(
            mesh.n_cells(),"G_band","W/m2",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> Q(
            mesh.n_cells(),"Q_band","W/m3",1);
        auto br=solve_participating_radiation_variable_properties(
            mesh,geometry,temperature,G,Q,directions,p,controls,wall_intensity_bcs);
        result.band_results.push_back(br);
        result.converged=result.converged && br.converged;
        result.iterations=std::max(result.iterations,br.iterations);
        const double w=band.weight/wsum;
        for(std::size_t c=0;c<mesh.n_cells();++c) {
            irradiation(c)+=w*G(c);
            radiation_source(c)+=w*Q(c);
        }
    }
    return result;
}


// -----------------------------------------------------------------------------
// Deterministic Monte-Carlo S2S visibility kernel.
//
// This is the geometry-driven counterpart to the centroid estimator above.
// Triangles are used as the minimal mesh-independent representation. Rays are
// cosine-weighted, visibility is tested against blockers, and the result is a
// directly measurable view-factor estimate. Increasing samples is a controlled
// accuracy refinement.
// -----------------------------------------------------------------------------

struct RadiationTriangle {
    std::array<double,3> a{0,0,0};
    std::array<double,3> b{0,0,0};
    std::array<double,3> c{0,0,0};
};

inline std::array<double,3> radiation_cross(
    const std::array<double,3>& a,const std::array<double,3>& b)
{
    return {a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]};
}

inline double radiation_dot(
    const std::array<double,3>& a,const std::array<double,3>& b)
{
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline std::array<double,3> radiation_sub(
    const std::array<double,3>& a,const std::array<double,3>& b)
{
    return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};
}

inline std::array<double,3> radiation_add(
    const std::array<double,3>& a,const std::array<double,3>& b)
{
    return {a[0]+b[0],a[1]+b[1],a[2]+b[2]};
}

inline std::array<double,3> radiation_scale(
    const std::array<double,3>& a,double s)
{
    return {s*a[0],s*a[1],s*a[2]};
}

inline double radiation_norm(const std::array<double,3>& a)
{
    return std::sqrt(radiation_dot(a,a));
}

inline double radiation_triangle_area(const RadiationTriangle& t)
{
    return 0.5*radiation_norm(radiation_cross(
        radiation_sub(t.b,t.a),radiation_sub(t.c,t.a)));
}

inline std::array<double,3> radiation_triangle_normal(const RadiationTriangle& t)
{
    const auto n=radiation_cross(radiation_sub(t.b,t.a),radiation_sub(t.c,t.a));
    const double l=radiation_norm(n);
    if(!(l>0.0)) throw std::invalid_argument("degenerate radiation triangle");
    return radiation_scale(n,1.0/l);
}

inline bool radiation_ray_triangle_hit(
    const std::array<double,3>& origin,
    const std::array<double,3>& direction,
    const RadiationTriangle& tri,
    double& distance)
{
    constexpr double eps=1e-11;
    const auto e1=radiation_sub(tri.b,tri.a);
    const auto e2=radiation_sub(tri.c,tri.a);
    const auto h=radiation_cross(direction,e2);
    const double det=radiation_dot(e1,h);
    if(std::abs(det)<eps) return false;
    const double inv=1.0/det;
    const auto s=radiation_sub(origin,tri.a);
    const double u=inv*radiation_dot(s,h);
    if(u<0.0 || u>1.0) return false;
    const auto q=radiation_cross(s,e1);
    const double v=inv*radiation_dot(direction,q);
    if(v<0.0 || u+v>1.0) return false;
    const double t=inv*radiation_dot(e2,q);
    if(t<=eps) return false;
    distance=t;
    return true;
}

inline double estimate_view_factor_ray_traced(
    const std::vector<RadiationTriangle>& source,
    const std::vector<RadiationTriangle>& target,
    const std::vector<RadiationTriangle>& blockers,
    std::size_t samples = 4096)
{
    if(source.empty() || target.empty() || samples==0)
        throw std::invalid_argument("ray-traced view factor requires non-empty surfaces");
    double source_area=0.0;
    for(const auto& t:source) source_area+=radiation_triangle_area(t);
    if(!(source_area>0.0)) throw std::invalid_argument("source surface has zero area");

    std::vector<double> cumulative;
    cumulative.reserve(source.size());
    double accum=0.0;
    for(const auto& t:source) {
        accum+=radiation_triangle_area(t);
        cumulative.push_back(accum);
    }

    std::size_t visible=0;
    for(std::size_t s=0;s<samples;++s) {
        const double u=(static_cast<double>(s)+0.5)/static_cast<double>(samples);
        const double v=std::fmod((static_cast<double>(s)*0.6180339887498949)+0.5,1.0);
        std::size_t ti=0;
        while(ti+1<cumulative.size() && u*source_area>cumulative[ti]) ++ti;
        const auto& tri=source[ti];
        const double su=std::sqrt(v);
        const double sv=std::fmod((static_cast<double>(s)*0.7548776662466927)+0.5,1.0);
        const auto p=radiation_add(
            radiation_add(radiation_scale(tri.a,1.0-su),
                          radiation_scale(tri.b,su*(1.0-sv))),
            radiation_scale(tri.c,su*sv));

        const auto n=radiation_triangle_normal(tri);
        std::array<double,3> ref{0,0,1};
        if(std::abs(n[2])>0.9) ref={1,0,0};
        auto tangent=radiation_cross(ref,n);
        tangent=radiation_scale(tangent,1.0/radiation_norm(tangent));
        auto bitangent=radiation_cross(n,tangent);

        const double r1=(static_cast<double>(s)+0.5)/static_cast<double>(samples);
        const double r2=std::fmod((static_cast<double>(s)*0.569840296)+0.25,1.0);
        const double phi=2.0*M_PI*r2;
        const double z=std::sqrt(1.0-r1);
        const double rho=std::sqrt(r1);
        auto dir=radiation_add(
            radiation_add(radiation_scale(tangent,rho*std::cos(phi)),
                          radiation_scale(bitangent,rho*std::sin(phi))),
            radiation_scale(n,z));
        const double dn=radiation_norm(dir);
        dir=radiation_scale(dir,1.0/dn);

        double nearest=std::numeric_limits<double>::infinity();
        bool target_hit=false;
        for(const auto& tt:target) {
            double d=0.0;
            if(radiation_ray_triangle_hit(radiation_add(p,radiation_scale(n,1e-9)),
                                           dir,tt,d) && d<nearest) {
                nearest=d; target_hit=true;
            }
        }
        if(!target_hit) continue;

        bool blocked=false;
        for(const auto& bt:blockers) {
            double d=0.0;
            if(radiation_ray_triangle_hit(radiation_add(p,radiation_scale(n,1e-9)),
                                           dir,bt,d) && d<nearest-1e-9) {
                blocked=true; break;
            }
        }
        if(!blocked) ++visible;
    }
    return static_cast<double>(visible)/static_cast<double>(samples);
}

} // namespace cfdx::physics

// -----------------------------------------------------------------------------
// Nonlinear Rosseland energy solve
// -----------------------------------------------------------------------------

struct RosselandSolveControls {
    std::size_t max_iterations = 100;
    double tolerance = 1e-8;
    double relaxation = 0.8;
    double minimum_temperature = 1.0;
    double absorption_floor = 1e-12;
};

struct RosselandSolveResult {
    bool converged = false;
    std::size_t iterations = 0;
    std::vector<double> temperature_residuals;
    std::vector<double> energy_balance_residuals;
};

inline RosselandSolveResult solve_rosseland_energy(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& source,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& absorption,
    const EnergySolverControls& energy_controls,
    const RosselandSolveControls& controls = {},
    const ScalarBoundaryConditions& bcs = {})
{
    if (temperature.size()!=mesh.n_cells() || source.size()!=mesh.n_cells() ||
        absorption.size()!=mesh.n_cells())
        throw std::invalid_argument("Rosseland field size mismatch");
    if (controls.max_iterations==0 || controls.tolerance<=0.0 ||
        controls.relaxation<=0.0 || controls.relaxation>1.0 ||
        controls.minimum_temperature<=0.0 ||
        controls.absorption_floor<=0.0)
        throw std::invalid_argument("invalid Rosseland controls");
    for (std::size_t c=0;c<mesh.n_cells();++c) {
        if (!std::isfinite(temperature(c)) || temperature(c)<controls.minimum_temperature ||
            !std::isfinite(absorption(c)) || absorption(c)<=controls.absorption_floor)
            throw std::invalid_argument("invalid Rosseland cell state");
    }

    RosselandSolveResult result;
    cfdx::core::Field<double,cfdx::core::Location::CELL> old=temperature;
    cfdx::core::Field<double,cfdx::core::Location::CELL> conductivity(
        mesh.n_cells(),"k_rad","W/m/K",1);
    cfdx::core::Field<double,cfdx::core::Location::CELL> su(
        mesh.n_cells(),"rosseland_source","W/m3",1);
    cfdx::core::Field<double,cfdx::core::Location::CELL> sp(
        mesh.n_cells(),"rosseland_sp","W/m3/K",1);
    cfdx::core::Field<double,cfdx::core::Location::FACE> zero_flux(
        mesh.n_faces(),"mass_flux","kg/s",1);
    zero_flux.fill(0.0);

    for (std::size_t iter=1;iter<=controls.max_iterations;++iter) {
        for(std::size_t c=0;c<mesh.n_cells();++c) {
            const double T=std::max(controls.minimum_temperature,temperature(c));
            conductivity(c)=rosseland_conductivity(T,absorption(c));
            su(c)=source(c);
            sp(c)=0.0;
        }

        std::vector<double> transient_diag(mesh.n_cells(),0.0);
        std::vector<double> transient_rhs(mesh.n_cells(),0.0);
        if (energy_controls.dt>0.0) {
            for(std::size_t c=0;c<mesh.n_cells();++c) {
                transient_diag[c]=energy_controls.density*energy_controls.cp*
                    geometry.cell_volumes[c]/energy_controls.dt;
                transient_rhs[c]=transient_diag[c]*old(c);
            }
        }

        auto eq=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,su,sp,bcs,true,nullptr,
            &transient_diag,&transient_rhs,&conductivity);

        cfdx::core::Vector candidate(mesh.n_cells(),0.0);
        for(std::size_t c=0;c<mesh.n_cells();++c) candidate(c)=temperature(c);
        ScalarSolveControls sc;
        sc.max_iterations=2000;
        sc.tolerance=energy_controls.tolerance;
        sc.relaxation=controls.relaxation;
        const auto lr=solve_scalar_equation(eq,candidate,sc);
        if(lr.status!=cfdx::core::SolverStatus::CONVERGED)
            throw std::runtime_error("Rosseland nonlinear iteration linear solve did not converge");

        double max_delta=0.0;
        for(std::size_t c=0;c<mesh.n_cells();++c) {
            const double bounded=std::max(controls.minimum_temperature,candidate(c));
            max_delta=std::max(max_delta,std::abs(bounded-temperature(c)));
            temperature(c)=bounded;
        }
        double scale=1.0;
        for(std::size_t c=0;c<mesh.n_cells();++c)
            scale=std::max(scale,std::abs(temperature(c)));
        const double rel=max_delta/scale;
        result.temperature_residuals.push_back(rel);
        result.energy_balance_residuals.push_back(energy_balance_relative(
            mesh,geometry,mass_flux,temperature,old,source,energy_controls,bcs));
        result.iterations=iter;
        if(rel<=controls.tolerance &&
           result.energy_balance_residuals.back()<=controls.tolerance) {
            result.converged=true;
            break;
        }
    }
    return result;
}

// -----------------------------------------------------------------------------
// Spatially varying gray participating-media transport.
// -----------------------------------------------------------------------------

inline RadiationSolveResult solve_participating_radiation_variable_properties(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& irradiation,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& radiation_source,
    const std::vector<DiscreteDirection>& directions,
    const RadiationOpticalPropertyField& properties,
    RadiationTransportControls controls = {},
    const ScalarBoundaryConditions& wall_intensity_bcs = {})
{
    if(properties.cells.size()!=mesh.n_cells())
        throw std::invalid_argument("radiation property field size mismatch");
    properties.validate();
    validate_discrete_directions(directions);
    if(temperature.size()!=mesh.n_cells() ||
       irradiation.size()!=mesh.n_cells() ||
       radiation_source.size()!=mesh.n_cells())
        throw std::invalid_argument("radiation field size mismatch");

    const std::size_t nc=mesh.n_cells();
    std::vector<cfdx::core::Field<double,cfdx::core::Location::CELL>> I;
    for(std::size_t m=0;m<directions.size();++m)
        I.emplace_back(nc,"I_"+std::to_string(m),"W/m2/sr",1);

    RadiationSolveResult result;
    for(std::size_t iter=1;iter<=controls.max_iterations;++iter) {
        auto old_source=radiation_source;
        std::vector<double> G(nc,0.0);
        for(std::size_t m=0;m<directions.size();++m)
            for(std::size_t c=0;c<nc;++c)
                G[c]+=directions[m].weight*I[m](c);

        double max_rel=0.0;
        for(std::size_t m=0;m<directions.size();++m) {
            cfdx::core::Field<double,cfdx::core::Location::FACE> directional_flux(
                mesh.n_faces(),"sI","m2",1);
            const auto& d=directions[m];
            for(std::size_t f=0;f<mesh.n_faces();++f)
                directional_flux(f)=d.dx*geometry.face_area_vectors[f].x+
                                    d.dy*geometry.face_area_vectors[f].y+
                                    d.dz*geometry.face_area_vectors[f].z;
            cfdx::core::Field<double,cfdx::core::Location::CELL> su(
                nc,"I_source","W/m3/sr",1);
            cfdx::core::Field<double,cfdx::core::Location::CELL> sp(
                nc,"I_sp","1/m",1);
            for(std::size_t c=0;c<nc;++c) {
                const auto& p=properties[c];
                su(c)=p.absorption*blackbody_intensity(temperature(c))+
                      p.scattering*G[c]/(4.0*M_PI);
                sp(c)=-(p.absorption+p.scattering);
            }
            auto eq=assemble_scalar_equation(mesh,geometry,directional_flux,0.0,
                                             su,sp,wall_intensity_bcs,true);
            cfdx::core::Vector candidate(nc,0.0);
            for(std::size_t c=0;c<nc;++c) candidate(c)=I[m](c);
            ScalarSolveControls sc;
            sc.max_iterations=controls.linear_max_iterations;
            sc.tolerance=controls.linear_tolerance;
            sc.relaxation=controls.intensity_relaxation;
            const auto lr=solve_scalar_equation(eq,candidate,sc);
            if(lr.status!=cfdx::core::SolverStatus::CONVERGED)
                throw std::runtime_error("variable-property DOM linear solve did not converge");
            for(std::size_t c=0;c<nc;++c) {
                if(!std::isfinite(candidate(c)) || candidate(c)<-1e-10)
                    throw std::runtime_error("variable-property DOM produced non-physical intensity");
                const double next=std::max(0.0,candidate(c));
                max_rel=std::max(max_rel,std::abs(next-I[m](c))/
                                      std::max(1.0,std::abs(next)));
                I[m](c)=next;
            }
        }

        for(std::size_t c=0;c<nc;++c) {
            double g=0.0;
            for(std::size_t m=0;m<directions.size();++m)
                g+=directions[m].weight*I[m](c);
            irradiation(c)=g;
            radiation_source(c)=properties[c].absorption*
                (4.0*M_PI*blackbody_intensity(temperature(c))-g);
        }
        double src_rel=0.0;
        for(std::size_t c=0;c<nc;++c)
            src_rel=std::max(src_rel,
                std::abs(radiation_source(c)-old_source(c))/
                std::max(1.0,std::abs(radiation_source(c))));
        result.history.push_back({iter,max_rel,src_rel});
        result.iterations=iter;
        if(max_rel<=controls.tolerance && src_rel<=controls.tolerance) {
            result.converged=true;
            break;
        }
    }
    return result;
}

