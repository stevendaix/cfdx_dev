#pragma once

#include "cfdx/physics/radiation.h"
#include "cfdx/physics/radiation_solver.h"
#include "cfdx/physics/radiation_models.h"
#include "cfdx/physics/finite_volume_transport.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
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

inline RadiationSolveResult solve_participating_radiation_variable_properties(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& irradiation,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& radiation_source,
    const std::vector<DiscreteDirection>& directions,
    const RadiationOpticalPropertyField& properties,
    RadiationTransportControls controls,
    const ScalarBoundaryConditions& wall_intensity_bcs);

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


// -----------------------------------------------------------------------------
// CPU BVH acceleration for S2S ray traversal.
// -----------------------------------------------------------------------------

struct RadiationAabb {
    std::array<double,3> lo{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    std::array<double,3> hi{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()};
};

inline RadiationAabb radiation_triangle_aabb(const RadiationTriangle& t)
{
    RadiationAabb b;
    for (int k=0;k<3;++k) {
        if (!std::isfinite(t.a[k]) || !std::isfinite(t.b[k]) ||
            !std::isfinite(t.c[k]))
            throw std::invalid_argument("radiation triangle coordinates must be finite");
        b.lo[k]=std::min({t.a[k],t.b[k],t.c[k]});
        b.hi[k]=std::max({t.a[k],t.b[k],t.c[k]});
    }
    return b;
}

inline RadiationAabb radiation_merge_aabb(const RadiationAabb& a,
                                          const RadiationAabb& b)
{
    RadiationAabb r;
    for (int k=0;k<3;++k) {
        r.lo[k]=std::min(a.lo[k],b.lo[k]);
        r.hi[k]=std::max(a.hi[k],b.hi[k]);
    }
    return r;
}

inline std::array<double,3> radiation_aabb_centroid(const RadiationAabb& b)
{
    return {(b.lo[0]+b.hi[0])*0.5,(b.lo[1]+b.hi[1])*0.5,
            (b.lo[2]+b.hi[2])*0.5};
}

inline bool radiation_ray_aabb_hit(
    const std::array<double,3>& origin,
    const std::array<double,3>& direction,
    const RadiationAabb& box,
    double max_distance)
{
    double tmin=0.0;
    double tmax=max_distance;
    for (int k=0;k<3;++k) {
        if (std::abs(direction[k])<1e-15) {
            if (origin[k]<box.lo[k] || origin[k]>box.hi[k]) return false;
            continue;
        }
        const double inv=1.0/direction[k];
        double t1=(box.lo[k]-origin[k])*inv;
        double t2=(box.hi[k]-origin[k])*inv;
        if (t1>t2) std::swap(t1,t2);
        tmin=std::max(tmin,t1);
        tmax=std::min(tmax,t2);
        if (tmin>tmax) return false;
    }
    return tmax>=std::max(tmin,1e-11);
}

class RadiationBvh {
public:
    struct Node {
        RadiationAabb bounds;
        std::uint32_t begin=0;
        std::uint32_t count=0;
        std::int32_t left=-1;
        std::int32_t right=-1;
        bool leaf() const { return left<0; }
    };

    RadiationBvh() = default;
    explicit RadiationBvh(const std::vector<RadiationTriangle>& triangles) {
        build(triangles);
    }
    RadiationBvh(std::vector<RadiationTriangle>&&) = delete;

    void build(const std::vector<RadiationTriangle>& triangles)
    {
        if (triangles.size() > std::numeric_limits<std::uint32_t>::max())
            throw std::length_error("radiation BVH triangle index overflow");
        for (const auto& triangle:triangles)
            (void)radiation_triangle_aabb(triangle);
        triangles_=&triangles;
        indices_.resize(triangles.size());
        for (std::size_t i=0;i<indices_.size();++i)
            indices_[i]=static_cast<std::uint32_t>(i);
        nodes_.clear();
        if (!indices_.empty()) build_node(0,indices_.size());
    }
    void build(std::vector<RadiationTriangle>&&) = delete;

    bool nearest_hit(const std::array<double,3>& origin,
                     const std::array<double,3>& direction,
                     double& distance) const
    {
        double direction_norm_squared=0.0;
        for (int k=0;k<3;++k) {
            if (!std::isfinite(origin[k]) || !std::isfinite(direction[k]))
                throw std::invalid_argument("radiation BVH ray must be finite");
            direction_norm_squared+=direction[k]*direction[k];
        }
        if (!(direction_norm_squared>0.0) || !std::isfinite(direction_norm_squared))
            throw std::invalid_argument("radiation BVH ray direction must be non-zero");
        distance=std::numeric_limits<double>::infinity();
        if (nodes_.empty()) return false;
        std::vector<std::int32_t> stack{0};
        bool hit=false;
        while (!stack.empty()) {
            const auto ni=stack.back();
            stack.pop_back();
            const auto& node=nodes_[ni];
            if (!radiation_ray_aabb_hit(origin,direction,node.bounds,distance))
                continue;
            if (node.leaf()) {
                for (std::uint32_t k=0;k<node.count;++k) {
                    const auto ti=indices_[node.begin+k];
                    double d=0.0;
                    if (radiation_ray_triangle_hit(origin,direction,
                                                   (*triangles_)[ti],d) &&
                        d<distance) {
                        distance=d;
                        hit=true;
                    }
                }
            } else {
                stack.push_back(node.left);
                stack.push_back(node.right);
            }
        }
        return hit;
    }

private:
    const std::vector<RadiationTriangle>* triangles_=nullptr;
    std::vector<std::uint32_t> indices_;
    std::vector<Node> nodes_;

    std::int32_t build_node(std::size_t begin,std::size_t end)
    {
        const auto id=static_cast<std::int32_t>(nodes_.size());
        nodes_.push_back({});
        RadiationAabb bounds;
        RadiationAabb centroids;
        for (std::size_t k=begin;k<end;++k) {
            const auto b=radiation_triangle_aabb((*triangles_)[indices_[k]]);
            bounds=radiation_merge_aabb(bounds,b);
            const auto c=radiation_aabb_centroid(b);
            centroids=radiation_merge_aabb(centroids,RadiationAabb{c,c});
        }
        nodes_[id].bounds=bounds;
        nodes_[id].begin=static_cast<std::uint32_t>(begin);
        nodes_[id].count=static_cast<std::uint32_t>(end-begin);
        constexpr std::size_t leaf_size=8;
        if (end-begin<=leaf_size) return id;

        const std::array<double,3> extent{
            centroids.hi[0]-centroids.lo[0],
            centroids.hi[1]-centroids.lo[1],
            centroids.hi[2]-centroids.lo[2]};
        int axis=0;
        if (extent[1]>extent[axis]) axis=1;
        if (extent[2]>extent[axis]) axis=2;
        const auto mid=begin+(end-begin)/2;
        std::nth_element(indices_.begin()+begin,indices_.begin()+mid,
                         indices_.begin()+end,[&](std::uint32_t a,
                                                   std::uint32_t b) {
            const auto ca=radiation_aabb_centroid(
                radiation_triangle_aabb((*triangles_)[a]));
            const auto cb=radiation_aabb_centroid(
                radiation_triangle_aabb((*triangles_)[b]));
            return ca[axis]<cb[axis];
        });
        nodes_[id].left=build_node(begin,mid);
        nodes_[id].right=build_node(mid,end);
        nodes_[id].count=0;
        return id;
    }
};

inline double estimate_view_factor_ray_traced(
    const std::vector<RadiationTriangle>& source,
    const std::vector<RadiationTriangle>& target,
    const std::vector<RadiationTriangle>& blockers,
    std::size_t samples = 4096)
{
    if(source.empty() || target.empty() || samples==0)
        throw std::invalid_argument("ray-traced view factor requires non-empty surfaces");

    RadiationBvh target_bvh(target);
    RadiationBvh blocker_bvh(blockers);
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

        const auto ray_origin=radiation_add(p,radiation_scale(n,1e-9));
        double nearest=std::numeric_limits<double>::infinity();
        if(!target_bvh.nearest_hit(ray_origin,dir,nearest))
            continue;
        double blocker_distance=std::numeric_limits<double>::infinity();
        const bool blocked=blocker_bvh.nearest_hit(
            ray_origin,dir,blocker_distance) &&
            blocker_distance < nearest-1e-9;
        if(!blocked) ++visible;
    }
    return static_cast<double>(visible)/static_cast<double>(samples);
}


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
            conductivity(c)=energy_controls.conductivity+
                rosseland_conductivity(T,absorption(c));
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

        std::vector<double> conductivity_values(mesh.n_cells(),0.0);
        for(std::size_t c=0;c<mesh.n_cells();++c)
            conductivity_values[c]=conductivity(c);
        const auto enthalpy_flux=
            enthalpy_face_flux(mass_flux,energy_controls.cp);
        auto eq=assemble_scalar_equation(
            mesh,geometry,enthalpy_flux,0.0,su,sp,bcs,true,nullptr,
            &transient_diag,&transient_rhs,&conductivity_values);

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
            mesh,geometry,mass_flux,temperature,old,source,energy_controls,bcs,
            nullptr,&conductivity_values));
        result.iterations=iter;
        if(rel<=controls.tolerance &&
           result.energy_balance_residuals.back()<=controls.tolerance) {
            result.converged=true;
            break;
        }
    }
    return result;
}

// Rosseland solve using the reduced transport extinction
// kappa_R = kappa_a + sigma_s (1-g). The isotropic-scattering default is g=0.
// This wrapper keeps the existing absorption-only API unambiguous while making
// the transport-opacity convention explicit for scattering media.
inline RosselandSolveResult solve_rosseland_energy_with_scattering(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& source,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& absorption,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& scattering,
    double asymmetry_factor,
    const EnergySolverControls& energy_controls,
    const RosselandSolveControls& controls = {},
    const ScalarBoundaryConditions& bcs = {})
{
    if (scattering.size()!=mesh.n_cells())
        throw std::invalid_argument("Rosseland scattering field size mismatch");
    if (!std::isfinite(asymmetry_factor) ||
        asymmetry_factor < -1.0 || asymmetry_factor > 1.0)
        throw std::invalid_argument("Rosseland asymmetry factor must be in [-1,1]");

    cfdx::core::Field<double,cfdx::core::Location::CELL> transport_opacity(
        mesh.n_cells(),"kappa_R","1/m",1);
    for (std::size_t c=0;c<mesh.n_cells();++c) {
        if (!std::isfinite(absorption(c)) || absorption(c)<=controls.absorption_floor ||
            !std::isfinite(scattering(c)) || scattering(c)<0.0)
            throw std::invalid_argument("invalid Rosseland absorption/scattering state");
        const double reduced_scattering=scattering(c)*(1.0-asymmetry_factor);
        transport_opacity(c)=absorption(c)+reduced_scattering;
        if (!(transport_opacity(c)>controls.absorption_floor) ||
            !std::isfinite(transport_opacity(c)))
            throw std::invalid_argument("invalid Rosseland transport opacity");
    }

    return solve_rosseland_energy(
        mesh,geometry,mass_flux,temperature,source,transport_opacity,
        energy_controls,controls,bcs);
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



} // namespace cfdx::physics
