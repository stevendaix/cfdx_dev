#pragma once
// V33 standalone verification oracle only (NOT integrated FV discretisation).
// A future CPU/GPU production port needs preallocated CSR stencils, fixed-size
// QR workspaces, and tests of full face and near-wall coupled operators.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct SSTV33Vec3 { double x=0, y=0, z=0; };
inline double v33_dot(SSTV33Vec3 a,SSTV33Vec3 b) {
    return a.x*b.x+a.y*b.y+a.z*b.z;
}
inline double v33_norm(SSTV33Vec3 a) {return std::sqrt(v33_dot(a,a));}
inline bool v33_finite(SSTV33Vec3 a) {
    return std::isfinite(a.x)&&std::isfinite(a.y)&&std::isfinite(a.z);
}
struct SSTV33WlsSample {
    SSTV33Vec3 offset;
    double value_delta=0.0;
    double weight=1.0;
};
struct SSTV33WlsResult {
    SSTV33Vec3 gradient;
    double relative_pivot=0.0;
};
// Column-scaled, twice-reorthogonalised modified Gram-Schmidt.
// Avoids squaring the condition number as normal equations do.
// Do not force 3D rank onto an intrinsically 2D/extruded mesh.
inline SSTV33WlsResult sst_v33_scaled_qr_gradient(
    const std::vector<SSTV33WlsSample>& samples,int dimension)
{
    if ((dimension!=2&&dimension!=3) ||
        samples.size()<static_cast<std::size_t>(dimension))
        throw std::invalid_argument("V33 WLS-QR: invalid intrinsic dimension or stencil");
    std::array<double,3> scales{0.,0.,0.};
    const auto dim=static_cast<std::size_t>(dimension);
    const auto m=samples.size();
    for (const auto& s:samples) {
        if (!v33_finite(s.offset) || !std::isfinite(s.value_delta) ||
            !std::isfinite(s.weight) || s.weight<=0)
            throw std::invalid_argument("V33 WLS-QR: nonfinite or invalid sample");
        std::array<double,3> d{s.offset.x,s.offset.y,s.offset.z};
        for (std::size_t j=0;j<dim;++j) scales[j]=std::max(scales[j],std::abs(d[j]));
    }
    for (std::size_t j=0;j<dim;++j)
        if (scales[j]<=0.)
            throw std::domain_error("V33 WLS-QR: missing coordinate rank");
    std::vector<std::array<double,3>> q(m);
    std::vector<double> rhs(m);
    for (std::size_t i=0;i<m;++i) {
        const auto w=std::sqrt(samples[i].weight);
        std::array<double,3> d{samples[i].offset.x,samples[i].offset.y,samples[i].offset.z};
        for(std::size_t j=0;j<dim;++j) q[i][j]=w*d[j]/scales[j];
        rhs[i]=w*samples[i].value_delta;
    }
    double r[3][3]{};
    double min_pivot=std::numeric_limits<double>::infinity(),max_pivot=0.;
    for (std::size_t j=0;j<dim;++j) {
        for (int pass=0;pass<2;++pass) for (std::size_t k=0;k<j;++k) {
            double projection=0.;
            for (std::size_t i=0;i<m;++i) projection+=q[i][j]*q[i][k];
            r[k][j]+=projection;
            for (std::size_t i=0;i<m;++i) q[i][j]-=projection*q[i][k];
        }
        double sq=0.;
        for (std::size_t i=0;i<m;++i) sq+=q[i][j]*q[i][j];
        r[j][j]=std::sqrt(sq);
        if (!(r[j][j]>1e-12)||!std::isfinite(r[j][j]))
            throw std::domain_error("V33 WLS-QR: rank-deficient stencil; extend neighbours");
        min_pivot=std::min(min_pivot,r[j][j]);max_pivot=std::max(max_pivot,r[j][j]);
        for (std::size_t i=0;i<m;++i) q[i][j]/=r[j][j];
    }
    if (min_pivot/max_pivot<=1e-10)
        throw std::domain_error("V33 WLS-QR: near-degenerate stencil");
    std::array<double,3> proj{},z{};
    for (std::size_t j=0;j<dim;++j)
        for (std::size_t i=0;i<m;++i) proj[j]+=q[i][j]*rhs[i];
    for (std::size_t j=dim; j-- > 0;) {
        double b=proj[j];
        for (std::size_t k=j+1;k<dim;++k) b-=r[j][k]*z[k];
        z[j]=b/r[j][j];
    }
    SSTV33WlsResult out;
    out.gradient={z[0]/scales[0],z[1]/scales[1],
                  dimension==3 ? z[2]/scales[2]:0.};
    out.relative_pivot=min_pivot/max_pivot;
    if (!v33_finite(out.gradient))
        throw std::domain_error("V33 WLS-QR: invalid reconstructed gradient");
    return out;
}
struct SSTV33DiffusionFace {
    double owner_flux=0.;
    double neighbour_flux=0.;
    double two_point_conductance=0.;
    double nonorthogonal_correction=0.;
};
// Shares one signed flux between owner and neighbour; only a mathematical
// oracle. Do not drop the skew correction in the production solver.
inline SSTV33DiffusionFace sst_v33_conservative_diffusion(
    SSTV33Vec3 center_delta,SSTV33Vec3 face_area,
    SSTV33Vec3 face_gradient,double owner_value,
    double neighbour_value,double diffusivity,double min_cosine=1e-6)
{
    if (!v33_finite(center_delta)||!v33_finite(face_area)||
        !v33_finite(face_gradient)||!std::isfinite(owner_value)||
        !std::isfinite(neighbour_value)||!std::isfinite(diffusivity)||
        diffusivity<0||!std::isfinite(min_cosine)||min_cosine<=0||min_cosine>1)
        throw std::invalid_argument("V33 diffusion: invalid face data");
    const auto n=v33_norm(center_delta),a=v33_norm(face_area);
    const auto projected=v33_dot(face_area,center_delta);
    if (n<=0||a<=0||projected<=0||projected/(n*a)<min_cosine)
        throw std::domain_error("V33 diffusion: inverted or too non-orthogonal face");
    const auto two_point=diffusivity*a*a/projected;
    const auto projected_delta=(a*a)/projected;
    const auto tangent=SSTV33Vec3{
        face_area.x-projected_delta*center_delta.x,
        face_area.y-projected_delta*center_delta.y,
        face_area.z-projected_delta*center_delta.z};
    SSTV33DiffusionFace out;
    out.two_point_conductance=two_point;
    out.nonorthogonal_correction=-diffusivity*v33_dot(face_gradient,tangent);
    out.owner_flux=two_point*(owner_value-neighbour_value)+out.nonorthogonal_correction;
    out.neighbour_flux=-out.owner_flux;
    if (!std::isfinite(out.owner_flux)||!std::isfinite(out.two_point_conductance))
        throw std::domain_error("V33 diffusion: overflow");
    return out;
}
}
