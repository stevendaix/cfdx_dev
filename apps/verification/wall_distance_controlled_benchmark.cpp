#include "cfdx/physics/wall_distance.h"

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace cfdx::physics;

namespace {

void add_triangle(WallSurface& s, std::size_t a, std::size_t b, std::size_t c)
{
    s.triangles.push_back({{a, b, c}});
}

void add_box(WallSurface& s, const WallDistanceVec3& lo, const WallDistanceVec3& hi)
{
    const std::size_t b = s.points.size();
    s.points.insert(s.points.end(), {
        {lo.x, lo.y, lo.z}, {hi.x, lo.y, lo.z},
        {hi.x, hi.y, lo.z}, {lo.x, hi.y, lo.z},
        {lo.x, lo.y, hi.z}, {hi.x, lo.y, hi.z},
        {hi.x, hi.y, hi.z}, {lo.x, hi.y, hi.z}
    });

    const std::array<std::array<std::size_t, 4>, 6> q = {{
        {{0, 3, 2, 1}}, {{4, 5, 6, 7}}, {{0, 1, 5, 4}},
        {{1, 2, 6, 5}}, {{2, 3, 7, 6}}, {{3, 0, 4, 7}}
    }};

    for (const auto& f : q) {
        add_triangle(s, b + f[0], b + f[1], b + f[2]);
        add_triangle(s, b + f[0], b + f[2], b + f[3]);
    }
}

bool inside_box(const WallDistanceVec3& p)
{
    // The surface itself belongs to the solid side. This is essential when a
    // validation grid lands exactly on a wall: wall nodes must not become
    // fluid unknowns with a zero-distance BVH query.
    return p.x >= 0.0 && p.x <= 1.0 &&
           p.y >= 0.0 && p.y <= 1.0 &&
           p.z >= 0.0 && p.z <= 1.0;
}

void add_plane_x(WallSurface& s,double x,double y0,double y1,double z0,double z1)
{
    const std::size_t b=s.points.size();
    s.points.insert(s.points.end(),{
        {x,y0,z0},{x,y1,z0},{x,y1,z1},{x,y0,z1}});
    add_triangle(s,b+0,b+1,b+2);
    add_triangle(s,b+0,b+2,b+3);
}

struct AnalyticCaseRow {
    std::string name;
    std::size_t n{};
    double l2{};
    double linf{};
    double residual{};
    std::size_t iterations{};
    bool converged{};
};

template<class Inside, class Reference>
AnalyticCaseRow run_analytic_poisson_case(const std::string& name,
                                          std::size_t n,
                                          const WallSurface& surface,
                                          Inside inside,
                                          Reference reference)
{
    const double h=1.0/static_cast<double>(n-1);
    const auto grid=make_wall_distance_grid(
        n,n,n,{-1.0,-0.5,-0.5},{h,h,h},inside);
    const auto bvh=WallDistanceBvh(surface);
    const auto audit=audit_poisson_operator(bvh,grid);

    if(audit.symmetry_error>1e-14)
        throw std::runtime_error(name+": Poisson operator is not symmetric");

    std::size_t iterations=0;
    double residual=0.0;
    const auto d=poisson_distance(bvh,grid,1000,1.5,&iterations,&residual);

    double sum2=0.0,ref2=0.0,maxe=0.0,maxref=0.0;
    for(std::size_t id=0;id<grid.points.size();++id) {
        if(grid.solid[id]) continue;
        const double exact=reference(grid.points[id]);
        const double e=std::abs(d[id]-exact);
        sum2+=e*e; ref2+=exact*exact;
        maxe=std::max(maxe,e); maxref=std::max(maxref,std::abs(exact));
    }

    return {
        name,n,
        std::sqrt(sum2/std::max(1e-30,ref2)),
        maxe/std::max(1e-30,maxref),
        residual,iterations,
        std::isfinite(residual) && residual<1e-10
    };
}

void run_analytic_validation()
{
    std::vector<double> half_l2,half_linf,channel_l2,channel_linf;

    // Plane half-space: -phi''=1, phi=0 at x=0, zero flux at x=-1.
    // The continuous solution is phi=(-x)-x^2/2 and Spalding's
    // reconstruction gives d=-x exactly. The node-centred outer Neumann
    // boundary and one-sided wall gradient make the implemented distance
    // reconstruction first-order in Linf and >1st-order in L2; these rates
    // are therefore the validation gate rather than an artificially tiny
    // absolute error tolerance.
    for(const std::size_t n : {9u,17u,33u,65u}) {
        WallSurface plane;
        add_plane_x(plane,0.0,-1.0,1.0,-1.0,1.0);
        const auto row=run_analytic_poisson_case(
            "half_space",n,plane,
            [](const WallDistanceVec3& p){ return p.x>=0.0; },
            [](const WallDistanceVec3& p){ return -p.x; });
        half_l2.push_back(row.l2);
        half_linf.push_back(row.linf);
        std::cout << "analytic,half_space," << row.n << ","
                  << row.l2 << "," << row.linf << ","
                  << row.residual << "," << row.iterations << ","
                  << (row.converged?"true":"false") << "\\n";
        if(!row.converged)
            throw std::runtime_error("half-space Poisson linear solve did not converge");
    }

    // Parallel channel: -phi''=1, phi=0 at x=0 and x=1.
    // phi=x(1-x)/2 and the reconstruction tends to min(x,1-x).
    for(const std::size_t n : {9u,17u,33u,65u}) {
        WallSurface channel;
        add_plane_x(channel,0.0,-1.0,1.0,-1.0,1.0);
        add_plane_x(channel,1.0,-1.0,1.0,-1.0,1.0);

        const double h=1.0/static_cast<double>(n-1);
        const auto grid=make_wall_distance_grid(
            n+2,n,n,{-h,-0.5,-0.5},{h,h,h},
            [](const WallDistanceVec3& p) {
                return p.x<=0.0 || p.x>=1.0;
            });
        const auto bvh=WallDistanceBvh(channel);
        const auto audit=audit_poisson_operator(bvh,grid);
        if(audit.symmetry_error>1e-14)
            throw std::runtime_error("channel Poisson operator is not symmetric");

        std::size_t iterations=0; double residual=0.0;
        const auto d=poisson_distance(bvh,grid,1000,1.5,&iterations,&residual);
        double sum2=0.0,ref2=0.0,maxe=0.0,maxref=0.0;
        for(std::size_t id=0;id<grid.points.size();++id) {
            if(grid.solid[id]) continue;
            const double x=grid.points[id].x;
            const double exact=std::min(x,1.0-x);
            const double e=std::abs(d[id]-exact);
            sum2+=e*e; ref2+=exact*exact;
            maxe=std::max(maxe,e); maxref=std::max(maxref,std::abs(exact));
        }
        const double l2=std::sqrt(sum2/std::max(1e-30,ref2));
        const double linf=maxe/std::max(1e-30,maxref);
        channel_l2.push_back(l2);
        channel_linf.push_back(linf);
        std::cout << "analytic,channel," << n << ","
                  << l2 << "," << linf << ","
                  << residual << "," << iterations << ","
                  << (std::isfinite(residual)&&residual<1e-10?"true":"false") << "\\n";
        if(!(std::isfinite(residual)&&residual<1e-10))
            throw std::runtime_error("channel Poisson linear solve did not converge");
    }

    auto order=[](double e0,double e1) {
        return std::log(e0/std::max(e1,1e-300))/std::log(2.0);
    };
    double half_l2_min=std::numeric_limits<double>::infinity();
    double half_linf_min=std::numeric_limits<double>::infinity();
    double channel_l2_min=std::numeric_limits<double>::infinity();
    double channel_linf_min=std::numeric_limits<double>::infinity();
    for(std::size_t i=1;i<half_l2.size();++i) {
        half_l2_min=std::min(half_l2_min,order(half_l2[i-1],half_l2[i]));
        half_linf_min=std::min(half_linf_min,order(half_linf[i-1],half_linf[i]));
        channel_l2_min=std::min(channel_l2_min,order(channel_l2[i-1],channel_l2[i]));
        channel_linf_min=std::min(channel_linf_min,order(channel_linf[i-1],channel_linf[i]));
    }
    std::cout << "analytic_order,half_space,l2=" << half_l2_min
              << ",linf=" << half_linf_min << "\\n";
    std::cout << "analytic_order,channel,l2=" << channel_l2_min
              << ",linf=" << channel_linf_min << "\\n";

    // The solve itself is second-order in the interior; the current
    // wall-gradient reconstruction is first-order at the first fluid node.
    // The gates below require a genuine decreasing error and a positive
    // asymptotic order without masking failures through relaxed residuals.
    if(half_l2_min<1.1 || half_linf_min<0.8)
        throw std::runtime_error("half-space Poisson convergence order is insufficient");
    if(channel_l2_min<0.8 || channel_linf_min<0.8)
        throw std::runtime_error("channel Poisson convergence order is insufficient");
}

} // namespace

int main(int argc, char** argv)
{
    std::cout << "CFDX wall-distance analytic Poisson validation\n";
    run_analytic_validation();

    const std::string output = argc > 1
        ? argv[1]
        : "wall_distance_controlled_benchmark.csv";
    const auto sizes = parse_sizes(argc, argv);

    std::vector<Row> rows;
    bool poisson_failure = false;

    for (const auto n : sizes) {
        const double lo = -0.5;
        const double hi = 1.5;
        const double h = (hi - lo) / static_cast<double>(n - 1);

        WallSurface surface;
        add_box(surface, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});

        const auto grid = make_wall_distance_grid(
            n, n, n, {lo, lo, lo}, {h, h, h}, inside_box);

        const auto geometric_reference = exact_reference(surface, grid);

        // Independent analytical reference for the same exact planar box.
        std::vector<double> analytic(grid.points.size(), std::numeric_limits<double>::infinity());
        for (std::size_t i = 0; i < grid.points.size(); ++i)
            if (!grid.solid[i])
                analytic[i] = analytic_box_distance(grid.points[i]);

        double reference_mismatch = 0.0;
        for (std::size_t i = 0; i < grid.points.size(); ++i) {
            if (grid.solid[i]) continue;
            reference_mismatch = std::max(
                reference_mismatch,
                std::abs(geometric_reference[i] - analytic[i]));
        }

        if (reference_mismatch > 1e-12) {
            std::cerr << "FAIL analytical/geometric reference mismatch at N="
                      << n << ": " << reference_mismatch << "\n";
            return 2;
        }

        const WallDistanceBvh bvh(surface);
        const auto op = audit_poisson_operator(bvh,grid);
        std::cout << "N=" << n << ", h=" << std::setprecision(12) << h
                  << ", samples=" << grid.points.size() << "\n"
                  << "Poisson operator: fluid=" << op.fluid_nodes
                  << ", fluid_fluid_faces=" << op.fluid_fluid_faces
                  << ", solid_faces=" << op.solid_faces
                  << ", outer_faces=" << op.outer_faces
                  << ", diag=[" << op.min_diagonal << "," << op.max_diagonal << "]"
                  << ", symmetry_error=" << op.symmetry_error << "\n";

        for (const auto method : methods()) {
            const auto t0 = std::chrono::steady_clock::now();
            const auto result = compute_wall_distance(method, surface, grid, 500);
            const auto t1 = std::chrono::steady_clock::now();

            const auto m = compare_wall_distance(
                grid, analytic, result.distance, 2.0 * h);

            double l1 = 0.0;
            double denom = 0.0;
            for (std::size_t i = 0; i < analytic.size(); ++i) {
                if (grid.solid[i] || !std::isfinite(analytic[i]) ||
                    !std::isfinite(result.distance[i]))
                    continue;
                l1 += std::abs(result.distance[i] - analytic[i]);
                denom += std::abs(analytic[i]);
            }
            l1 /= std::max(1e-30, denom);

            const double ms =
                std::chrono::duration<double, std::milli>(t1 - t0).count();

            if (method == WallDistanceMethod::POISSON && !result.converged) {
                poisson_failure = true;
                std::cerr << "FAIL Poisson did not converge at N=" << n
                          << ": residual_inf=" << result.residual_inf
                          << ", iterations=" << result.iterations
                          << ", stopping_reason=" << result.stopping_reason << "\n";
            }

            rows.push_back({
                n,
                result.method,
                l1,
                m.l2_relative,
                m.linf_relative,
                m.near_wall_l2_relative,
                result.residual_inf,
                m.monotonicity_violations,
                result.iterations,
                result.converged,
                result.stopping_reason,
                ms
            });

            std::cout << result.method << ", "
                      << "L1=" << l1
                      << ", L2=" << m.l2_relative
                      << ", Linf=" << m.linf_relative
                      << ", nearL2=" << m.near_wall_l2_relative
                      << ", residual=" << result.residual_inf
                      << ", violations=" << m.monotonicity_violations
                      << ", iterations=" << result.iterations
                      << ", converged=" << (result.converged ? "true" : "false")
                      << ", stopping_reason=" << result.stopping_reason
                      << ", ms=" << ms << "\n";
        }
    }

    std::ofstream csv(output);
    if (!csv) {
        std::cerr << "cannot open output: " << output << "\n";
        return 3;
    }

    csv << "N,h,method,l1_relative,l2_relative,linf_relative,"
           "near_wall_l2_relative,residual_inf,monotonicity_violations,"
           "iterations,converged,stopping_reason,time_ms\n";

    for (const auto& r : rows) {
        const double h = 2.0 / static_cast<double>(r.n - 1);
        csv << r.n << ',' << h << ',' << r.method << ','
            << std::setprecision(16)
            << r.l1 << ',' << r.l2 << ',' << r.linf << ','
            << r.near_l2 << ',' << r.residual << ','
            << r.violations << ',' << r.iterations << ','
            << (r.converged ? "true" : "false") << ',' << r.stopping_reason << ','
            << r.time_ms
            << '\n';
    }

    return poisson_failure ? 4 : 0;
}
