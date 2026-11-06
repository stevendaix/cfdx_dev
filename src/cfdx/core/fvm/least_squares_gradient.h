#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/numerics/gradient_stencil.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::core {

enum class GradientWeighting {
    UNIFORM,
    INVERSE_DISTANCE,
    INVERSE_DISTANCE_SQUARED
};

// Forward declaration: the legacy LS entry point delegates to the same
// rank-revealing one-sided Jacobi factorization as the robust WLS method.
// requested_dimension=0 is private automatic-rank mode (rank 1/2/3), needed
// for legitimate lower-dimensional stencils embedded in 3-D.
inline Vec3 weighted_least_squares_gradient_robust(
    const Vec3& centre,
    double value,
    const std::vector<Vec3>& neighbour_centres,
    const std::vector<double>& neighbour_values,
    GradientWeighting weighting,
    int requested_dimension,
    StencilQuality* quality);

/**
 * Minimum-norm cell-centred least-squares gradient, including valid rank-1
 * and rank-2 stencils in any 3-D orientation. Previous normal-equation
 * elimination set a missing Cartesian component to zero, so a rotated 2-D
 * stencil (z=x) incorrectly returned (2,2,0) instead of (1,2,1).
 *
 * Use inverse-distance-squared row weights as before, but solve the weighted
 * DESIGN matrix by one-sided Jacobi rather than rank-deficient normal
 * equations. No singular direction is assigned an arbitrary Cartesian axis.
 * Rank and conditioning remain observable through StencilQuality.
 */
inline Vec3 least_squares_gradient(
    const Vec3& centre,
    double value,
    const std::vector<Vec3>& neighbour_centres,
    const std::vector<double>& neighbour_values,
    StencilQuality* quality)
{
    return weighted_least_squares_gradient_robust(
        centre, value, neighbour_centres, neighbour_values,
        GradientWeighting::INVERSE_DISTANCE_SQUARED, 0, quality);
}

inline Vec3 least_squares_gradient(
    const Vec3& centre,
    double value,
    const std::vector<Vec3>& neighbour_centres,
    const std::vector<double>& neighbour_values)
{
    return least_squares_gradient(centre, value, neighbour_centres, neighbour_values, nullptr);
}

// Explicit robust WLS: one-sided Jacobi decomposition of the weighted design
// matrix (not normal equations). Dimension=2 returns the minimum-norm
// tangent-plane gradient on an embedded planar mesh; dimension=3 requires a
// full 3-D stencil. Dimension=0 is reserved for the legacy LS wrapper and
// deliberately infers rank 1/2/3 to preserve embedded-mesh compatibility.
inline Vec3 weighted_least_squares_gradient_robust(
    const Vec3& centre,
    double value,
    const std::vector<Vec3>& neighbour_centres,
    const std::vector<double>& neighbour_values,
    GradientWeighting weighting,
    int requested_dimension,
    StencilQuality* quality = nullptr)
{
    if (requested_dimension != 0 && requested_dimension != 2 &&
        requested_dimension != 3)
        throw std::invalid_argument("robust WLS: dimension must be 0, 2 or 3");
    if (neighbour_centres.size() != neighbour_values.size())
        throw std::invalid_argument("robust WLS: size mismatch");
    StencilQuality q;
    const bool automatic_rank = requested_dimension == 0;
    q.dimension = automatic_rank ? 3 : requested_dimension;
    q.samples = neighbour_centres.size();
    auto fail = [&](const char* message) -> void {
        if (quality) *quality = q;
        throw std::runtime_error(message);
    };
    if (!std::isfinite(value) || !std::isfinite(centre.x) ||
        !std::isfinite(centre.y) || !std::isfinite(centre.z)) {
        q.finite = false;
        fail("weighted_least_squares_gradient_planar: non-finite input");
    }

    std::vector<std::array<double, 3>> design;
    std::vector<double> rhs;
    design.reserve(q.samples);
    rhs.reserve(q.samples);
    for (std::size_t i = 0; i < q.samples; ++i) {
        const Vec3 d = neighbour_centres[i] - centre;
        const double r2 = d.mag2();
        const double dv = neighbour_values[i] - value;
        if (!(r2 > 0.0) || !std::isfinite(r2) || !std::isfinite(dv)) {
            q.finite = false;
            fail("weighted_least_squares_gradient_planar: non-finite or duplicate sample");
        }
        double factor = 1.0;
        if (weighting == GradientWeighting::INVERSE_DISTANCE)
            factor = 1.0 / std::sqrt(std::sqrt(r2));
        else if (weighting == GradientWeighting::INVERSE_DISTANCE_SQUARED)
            factor = 1.0 / std::sqrt(r2);
        if (!std::isfinite(factor)) {
            q.finite = false;
            fail("weighted_least_squares_gradient_planar: invalid weight");
        }
        design.push_back({factor*d.x, factor*d.y, factor*d.z});
        rhs.push_back(factor*dv);
    }

    // Right-singular vector matrix. Orthogonalize all three design columns
    // by one-sided Jacobi sweeps; no squared condition number is formed.
    std::array<std::array<double,3>,3> V{{
        {{1.0,0.0,0.0}}, {{0.0,1.0,0.0}}, {{0.0,0.0,1.0}}}};
    for (int sweep = 0; sweep < 40; ++sweep) {
        bool changed = false;
        for (int p = 0; p < 2; ++p) {
            for (int k = p+1; k < 3; ++k) {
                double alpha = 0.0, beta = 0.0, gamma = 0.0;
                for (const auto& row : design) {
                    alpha += row[p]*row[p];
                    beta += row[k]*row[k];
                    gamma += row[p]*row[k];
                }
                if (!(alpha > 0.0) || !(beta > 0.0)) continue;
                if (std::abs(gamma) <= 8.0 *
                    std::numeric_limits<double>::epsilon() * std::sqrt(alpha*beta))
                    continue;
                const double tau = (beta-alpha)/(2.0*gamma);
                const double t = std::copysign(1.0, tau) /
                    (std::abs(tau) + std::sqrt(1.0+tau*tau));
                const double c = 1.0/std::sqrt(1.0+t*t);
                const double sn = c*t;
                for (auto& row : design) {
                    const double old_p=row[p], old_k=row[k];
                    row[p] = c*old_p-sn*old_k;
                    row[k] = sn*old_p+c*old_k;
                }
                for (auto& row : V) {
                    const double old_p=row[p], old_k=row[k];
                    row[p] = c*old_p-sn*old_k;
                    row[k] = sn*old_p+c*old_k;
                }
                changed = true;
            }
        }
        if (!changed) break;
    }
    std::array<double,3> sigma{};
    for (int j=0;j<3;++j) {
        for (const auto& row:design) sigma[j] += row[j]*row[j];
        sigma[j] = std::sqrt(sigma[j]);
    }
    std::array<int,3> order{{0,1,2}};
    std::sort(order.begin(),order.end(),
              [&](int a,int b){return sigma[a]>sigma[b];});
    q.scale = sigma[order[0]];
    q.max_pivot = q.scale;
    const double cutoff = 128.0*std::numeric_limits<double>::epsilon()*q.scale;
    for (int j : order)
        if (sigma[j] > cutoff) ++q.rank;
    // Rounding of polyhedral centres can produce a tiny spurious third
    // singular value on a geometrically planar 2-D extrusion. Do not reject
    // an otherwise full-rank 2-D stencil for that numerical noise. True 3-D
    // stencils are still rejected explicitly (fail-closed) in 2-D mode.
    if (requested_dimension == 2 && q.rank == 3 &&
        sigma[order[2]] <= 1e-10 * sigma[order[0]])
        q.rank = 2;
    const int solve_rank = automatic_rank ? q.rank : requested_dimension;
    q.full_rank = (q.rank == q.dimension);
    q.min_pivot = q.full_rank && solve_rank > 0
        ? sigma[order[solve_rank-1]] : 0.0;
    q.condition_estimate = q.full_rank && q.min_pivot > 0.0
        ? q.max_pivot/q.min_pivot
        : std::numeric_limits<double>::infinity();
    if (quality) *quality = q;
    if (!automatic_rank && q.rank != requested_dimension)
        fail("robust WLS: stencil rank differs from requested dimension");
    std::array<double,3> g{{0.0,0.0,0.0}};
    for (int oi=0;oi<solve_rank;++oi) {
        const int j=order[oi];
        double proj=0.0;
        for (std::size_t row=0;row<design.size();++row)
            proj += design[row][j]*rhs[row];
        const double coeff = proj/(sigma[j]*sigma[j]);
        for (int k=0;k<3;++k) g[k] += V[k][j]*coeff;
    }
    if (!std::isfinite(g[0]) || !std::isfinite(g[1]) || !std::isfinite(g[2]))
        fail("robust WLS: non-finite gradient");
    return {g[0],g[1],g[2]};
}

inline Vec3 weighted_least_squares_gradient_planar(
    const Vec3& centre, double value,
    const std::vector<Vec3>& neighbour_centres,
    const std::vector<double>& neighbour_values,
    GradientWeighting weighting, StencilQuality* quality = nullptr)
{
    return weighted_least_squares_gradient_robust(
        centre, value, neighbour_centres, neighbour_values, weighting, 2, quality);
}

inline Vec3 weighted_least_squares_gradient(
    const Vec3& centre,
    double value,
    const std::vector<Vec3>& neighbour_centres,
    const std::vector<double>& neighbour_values,
    GradientWeighting weighting,
    StencilQuality* quality = nullptr,
    int requested_dimension = 3)
{
    if (requested_dimension == 2)
        return weighted_least_squares_gradient_planar(
            centre, value, neighbour_centres, neighbour_values, weighting, quality);
    if (requested_dimension != 3)
        throw std::invalid_argument("weighted_least_squares_gradient: dimension must be 2 or 3");
    if (neighbour_centres.size() != neighbour_values.size())
        throw std::invalid_argument("weighted_least_squares_gradient: neighbour dimensions do not match");

    StencilQuality q;
    q.samples = neighbour_centres.size();
    q.dimension = 3;
    double a[3][3] = {};
    double b[3] = {};

    for (std::size_t k = 0; k < neighbour_centres.size(); ++k) {
        const Vec3 d = neighbour_centres[k] - centre;
        const double r2 = d.mag2();
        if (!(r2 > 0.0) || !std::isfinite(r2)) {
            q.finite = false;
            continue;
        }
        const double r = std::sqrt(r2);
        double w = 1.0;
        if (weighting == GradientWeighting::INVERSE_DISTANCE) w = 1.0 / r;
        if (weighting == GradientWeighting::INVERSE_DISTANCE_SQUARED) w = 1.0 / r2;
        const double dv = neighbour_values[k] - value;
        a[0][0]+=w*d.x*d.x; a[0][1]+=w*d.x*d.y; a[0][2]+=w*d.x*d.z;
        a[1][0]+=w*d.y*d.x; a[1][1]+=w*d.y*d.y; a[1][2]+=w*d.y*d.z;
        a[2][0]+=w*d.z*d.x; a[2][1]+=w*d.z*d.y; a[2][2]+=w*d.z*d.z;
        b[0]+=w*d.x*dv; b[1]+=w*d.y*dv; b[2]+=w*d.z*dv;
    }

    double scale=0.0;
    for (int i=0;i<3;++i) scale=std::max(scale,std::abs(a[i][i]));
    q.scale=scale;
    if (!q.finite) {
        if (quality) *quality=q;
        throw std::runtime_error("weighted_least_squares_gradient: non-finite stencil geometry");
    }
    if (!(scale>0.0) || !std::isfinite(scale)) {
        if (quality) *quality=q;
        return {};
    }
    const double tol=128.0*std::numeric_limits<double>::epsilon()*scale;
    double m[3][4] = {};
    for(int i=0;i<3;++i){for(int j=0;j<3;++j)m[i][j]=a[i][j];m[i][3]=b[i];}
    double minp=std::numeric_limits<double>::infinity(), maxp=0.0;
    for(int col=0;col<3;++col){
        int p=col; double pa=std::abs(m[col][col]);
        for(int row=col+1;row<3;++row) if(std::abs(m[row][col])>pa){pa=std::abs(m[row][col]);p=row;}
        if(!(pa>tol)) continue;
        if(p!=col) for(int j=col;j<4;++j) std::swap(m[col][j],m[p][j]);
        ++q.rank;
        minp=std::min(minp,std::abs(m[col][col])); maxp=std::max(maxp,std::abs(m[col][col]));
        for(int row=col+1;row<3;++row){const double f=m[row][col]/m[col][col];for(int j=col;j<4;++j)m[row][j]-=f*m[col][j];}
    }
    q.min_pivot=std::isfinite(minp)?minp:0.0; q.max_pivot=maxp; q.full_rank=q.rank==3;
    q.condition_estimate=q.full_rank && q.min_pivot>0.0
        ? q.max_pivot/q.min_pivot
        : std::numeric_limits<double>::infinity();
    if(quality)*quality=q;
    if(q.rank<3) throw std::runtime_error("weighted_least_squares_gradient: rank-deficient stencil");
    double x[3]={};
    for(int i=2;i>=0;--i){double rhs=m[i][3];for(int j=i+1;j<3;++j)rhs-=m[i][j]*x[j];x[i]=rhs/m[i][i];}
    if(!std::isfinite(x[0])||!std::isfinite(x[1])||!std::isfinite(x[2]))
        throw std::runtime_error("weighted_least_squares_gradient: non-finite gradient");
    return {x[0],x[1],x[2]};
}

} // namespace cfdx::core
