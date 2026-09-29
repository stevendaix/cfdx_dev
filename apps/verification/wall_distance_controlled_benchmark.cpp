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
    return p.x >= 0.0 && p.x <= 1.0 &&
           p.y >= 0.0 && p.y <= 1.0 &&
           p.z >= 0.0 && p.z <= 1.0;
}

WallDistanceVec3 analytic_box_distance_gradient(const WallDistanceVec3& p)
{
    const WallDistanceVec3 q{
        std::clamp(p.x,0.0,1.0),
        std::clamp(p.y,0.0,1.0),
        std::clamp(p.z,0.0,1.0)};
    const WallDistanceVec3 v=p-q;
    const double d=wd_norm(v);
    if(!(d>1e-14)) return {0.0,0.0,0.0};
    return v*(1.0/d);
}

double reconstruct_from_gradient(double phi,const WallDistanceVec3& grad)
{
    const double gn=wd_norm(grad);
    return std::max(0.0,std::sqrt(std::max(0.0,gn*gn+2.0*phi))-gn);
}

double analytic_box_distance(const WallDistanceVec3& p)
{
    const double dx = std::min(std::abs(p.x), std::abs(p.x - 1.0));
    const double dy = std::min(std::abs(p.y), std::abs(p.y - 1.0));
    const double dz = std::min(std::abs(p.z), std::abs(p.z - 1.0));

    // For points outside the box, distance to the box is the Euclidean
    // distance to its closest point. This is exact for the planar box.
    const double ox = p.x < 0.0 ? -p.x : (p.x > 1.0 ? p.x - 1.0 : 0.0);
    const double oy = p.y < 0.0 ? -p.y : (p.y > 1.0 ? p.y - 1.0 : 0.0);
    const double oz = p.z < 0.0 ? -p.z : (p.z > 1.0 ? p.z - 1.0 : 0.0);
    if (ox > 0.0 || oy > 0.0 || oz > 0.0)
        return std::sqrt(ox * ox + oy * oy + oz * oz);

    return std::min({dx, dy, dz});
}

struct Row {
    std::size_t n{};
    std::string method;
    double l1{};
    double l2{};
    double linf{};
    double near_l2{};
    double residual{};
    double violations{};
    std::size_t iterations{};
    bool converged{};
    std::string stopping_reason;
    double time_ms{};
    std::size_t poisson_wall_ray_hits{};
    std::size_t poisson_wall_ray_misses{};
    std::size_t poisson_wall_fallbacks{};
    std::size_t poisson_wall_bad_alignment{};
    double poisson_wall_min_alignment{};
    double poisson_phi_min{};
    double poisson_phi_max{};
    double poisson_grad_min{};
    double poisson_grad_max{};
    double poisson_distance_l2_error{};
    double poisson_distance_linf_error{};
    double poisson_exact_gradient_l2_error{};
    double poisson_exact_gradient_linf_error{};
    double poisson_exact_phi_exact_gradient_l2_error{};
    double poisson_numerical_phi_exact_gradient_l2_error{};
    double poisson_exact_phi_numerical_gradient_l2_error{};
    double poisson_full_reconstruction_l2_error{};
    double poisson_exact_phi_exact_gradient_linf_error{};
    double poisson_numerical_phi_exact_gradient_linf_error{};
    double poisson_exact_phi_numerical_gradient_linf_error{};
    double poisson_full_reconstruction_linf_error{};
};

std::vector<WallDistanceMethod> methods()
{
    return {
        WallDistanceMethod::EXACT_GEOMETRIC,
        WallDistanceMethod::SEARCH_BASED,
        WallDistanceMethod::MESH_WAVE,
        WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
        WallDistanceMethod::POISSON,
        WallDistanceMethod::EIKONAL,
        WallDistanceMethod::HAMILTON_JACOBI,
        WallDistanceMethod::ADVECTION_DIFFUSION,
        WallDistanceMethod::HYBRID_POISSON_EIKONAL
    };
}

std::vector<std::size_t> parse_sizes(int argc, char** argv)
{
    std::vector<std::size_t> sizes;
    for (int i = 2; i < argc; ++i)
        sizes.push_back(static_cast<std::size_t>(std::stoul(argv[i])));

    if (sizes.empty())
        sizes = {16, 32, 64};

    return sizes;
}

} // namespace

int main(int argc, char** argv)
{
    const std::string output = argc > 1
        ? argv[1]
        : "wall_distance_controlled_benchmark.csv";
    const auto sizes = parse_sizes(argc, argv);

    std::vector<Row> rows;
    bool poisson_failure = false;
    bool advection_diffusion_failure = false;

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
        const auto op = audit_poisson_operator(bvh, grid);
        const auto offsets = audit_poisson_wall_offsets(bvh, grid);
        const std::size_t max_id=offsets.max_wall_coefficient_cell;
        const std::size_t max_k=max_id/(grid.nx*grid.ny);
        const std::size_t max_rem=max_id%(grid.nx*grid.ny);
        const std::size_t max_j=max_rem/grid.nx;
        const std::size_t max_i=max_rem%grid.nx;
        const auto max_p=grid.points[max_id];
        std::cout << "N=" << n << ", h=" << std::setprecision(12) << h
                  << ", samples=" << grid.points.size() << "\n"
                  << "Poisson operator: fluid=" << op.fluid_nodes
                  << ", fluid_fluid_faces=" << op.fluid_fluid_faces
                  << ", solid_faces=" << op.solid_faces
                  << ", outer_faces=" << op.outer_faces
                  << ", diag=[" << op.min_diagonal << "," << op.max_diagonal << "]"
                  << ", min_diagonal_dominance=" << op.min_diagonal_dominance
                  << ", symmetry_error=" << op.symmetry_error << "\n"
                  << "Poisson wall offsets: faces=" << offsets.solid_face_count
                  << ", delta=[" << offsets.min_delta << "," << offsets.max_delta << "]"
                  << ", delta_over_h=[" << offsets.min_delta_over_h << "," << offsets.max_delta_over_h << "]"
                  << ", max_wall_coefficient=" << offsets.max_wall_coefficient
                  << ", max_wall_coefficient_cell=" << offsets.max_wall_coefficient_cell
                  << ", max_wall_coefficient_ijk=" << max_i << ":" << max_j << ":" << max_k
                  << ", max_wall_coefficient_point=" << max_p.x << ":" << max_p.y << ":" << max_p.z
                  << ", max_wall_coefficient_delta=" << offsets.max_wall_coefficient_delta
                  << ", max_wall_coefficient_h=" << offsets.max_wall_coefficient_h
                  << ", degenerate=" << offsets.degenerate_count << "\n";

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

            double poisson_exact_gradient_l2=0.0;
            double poisson_exact_gradient_linf=0.0;
            if(method==WallDistanceMethod::POISSON) {
                // Re-run only the elliptic solve and replace the reconstructed
                // gradient by the exact geometric gradient of the box.  This
                // is a controlled error-budget diagnostic: it does not claim
                // that the Poisson phi has an exact analytical solution.
                std::size_t diagnostic_iterations=0;
                double diagnostic_residual=0.0;
                const auto diagnostic_phi=poisson_potential(
                    bvh,grid,500,1.5,&diagnostic_iterations,&diagnostic_residual);
                double e2=0.0,r2=0.0,em=0.0,rm=0.0;
                double exact_phi_exact_grad_e2=0.0, exact_phi_exact_grad_em=0.0;
                double num_phi_exact_grad_e2=0.0, num_phi_exact_grad_em=0.0;
                double exact_phi_num_grad_e2=0.0, exact_phi_num_grad_em=0.0;
                double full_e2=0.0, full_em=0.0;
                for(std::size_t id=0;id<grid.points.size();++id) {
                    if(grid.solid[id]) continue;
                    const auto exact_grad=analytic_box_distance_gradient(grid.points[id]);
                    const double d=analytic[id];
                    // Exact reconstruction input for the Tucker relation:
                    // d = sqrt(|grad d|^2 + 2 phi) - |grad d|.
                    // This is an input manufactured from the geometric
                    // distance; it is not claimed to solve the Poisson PDE
                    // exactly in edge/corner regions.
                    const double exact_phi=d-0.5*d*d;
                    const double d_exact_exact=reconstruct_from_gradient(exact_phi,exact_grad);
                    const double d_num_exact=reconstruct_from_gradient(diagnostic_phi[id],exact_grad);
                    const double gx=poisson_reconstruction_gradient_component(bvh,diagnostic_phi,grid,id,0);
                    const double gy=poisson_reconstruction_gradient_component(bvh,diagnostic_phi,grid,id,1);
                    const double gz=poisson_reconstruction_gradient_component(bvh,diagnostic_phi,grid,id,2);
                    const WallDistanceVec3 numerical_grad{gx,gy,gz};
                    const double d_exact_num=reconstruct_from_gradient(exact_phi,numerical_grad);
                    const double d_full=reconstruct_from_gradient(diagnostic_phi[id],numerical_grad);
                    const double e0=std::abs(d_exact_exact-d);
                    const double e1=std::abs(d_num_exact-d);
                    const double e2n=std::abs(d_exact_num-d);
                    const double e3=std::abs(d_full-d);
                    exact_phi_exact_grad_e2+=e0*e0;
                    exact_phi_exact_grad_em=std::max(exact_phi_exact_grad_em,e0);
                    num_phi_exact_grad_e2+=e1*e1;
                    num_phi_exact_grad_em=std::max(num_phi_exact_grad_em,e1);
                    exact_phi_num_grad_e2+=e2n*e2n;
                    exact_phi_num_grad_em=std::max(exact_phi_num_grad_em,e2n);
                    full_e2+=e3*e3;
                    full_em=std::max(full_em,e3);
                    r2+=d*d;
                    rm=std::max(rm,std::abs(d));
                    e2+=e1*e1;
                    em=std::max(em,e1);
                }
                poisson_exact_gradient_l2=std::sqrt(e2/std::max(1e-30,r2));
                poisson_exact_gradient_linf=em/std::max(1e-30,rm);
                poisson_exact_phi_exact_gradient_l2_error=std::sqrt(exact_phi_exact_grad_e2/std::max(1e-30,r2));
                poisson_numerical_phi_exact_gradient_l2_error=std::sqrt(num_phi_exact_grad_e2/std::max(1e-30,r2));
                poisson_exact_phi_numerical_gradient_l2_error=std::sqrt(exact_phi_num_grad_e2/std::max(1e-30,r2));
                poisson_full_reconstruction_l2_error=std::sqrt(full_e2/std::max(1e-30,r2));
                poisson_exact_phi_exact_gradient_linf_error=exact_phi_exact_grad_em/std::max(1e-30,rm);
                poisson_numerical_phi_exact_gradient_linf_error=num_phi_exact_grad_em/std::max(1e-30,rm);
                poisson_exact_phi_numerical_gradient_linf_error=exact_phi_num_grad_em/std::max(1e-30,rm);
                poisson_full_reconstruction_linf_error=full_em/std::max(1e-30,rm);
                std::cout << "Poisson reconstruction decomposition"
                          << ", exact_phi_exact_grad_L2=" << poisson_exact_phi_exact_gradient_l2_error
                          << ", numerical_phi_exact_grad_L2=" << poisson_numerical_phi_exact_gradient_l2_error
                          << ", exact_phi_numerical_grad_L2=" << poisson_exact_phi_numerical_gradient_l2_error
                          << ", full_L2=" << poisson_full_reconstruction_l2_error
                          << ", exact_phi_exact_grad_Linf=" << poisson_exact_phi_exact_gradient_linf_error
                          << ", numerical_phi_exact_grad_Linf=" << poisson_numerical_phi_exact_gradient_linf_error
                          << ", exact_phi_numerical_grad_Linf=" << poisson_exact_phi_numerical_gradient_linf_error
                          << ", full_Linf=" << poisson_full_reconstruction_linf_error
                          << ", diagnostic_residual=" << diagnostic_residual
                          << ", diagnostic_iterations=" << diagnostic_iterations << "\n";
            }

            if (method == WallDistanceMethod::POISSON && !result.converged) {
                poisson_failure = true;
                std::cerr << "FAIL Poisson did not converge at N=" << n
                          << ": residual_inf=" << result.residual_inf
                          << ", iterations=" << result.iterations
                          << ", stopping_reason=" << result.stopping_reason << "\n";
            }
            if (method == WallDistanceMethod::ADVECTION_DIFFUSION && !result.converged) {
                advection_diffusion_failure = true;
                std::cerr << "FAIL Advection-diffusion did not converge at N=" << n
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
                ms,
                result.poisson_wall_ray_hits,
                result.poisson_wall_ray_misses,
                result.poisson_wall_fallbacks,
                result.poisson_wall_bad_alignment,
                result.poisson_wall_min_alignment,
                result.poisson_phi_min,
                result.poisson_phi_max,
                result.poisson_grad_min,
                result.poisson_grad_max,
                result.poisson_distance_l2_error,
                result.poisson_distance_linf_error,
                poisson_exact_gradient_l2,
                poisson_exact_gradient_linf,
                poisson_exact_phi_exact_gradient_l2_error,
                poisson_numerical_phi_exact_gradient_l2_error,
                poisson_exact_phi_numerical_gradient_l2_error,
                poisson_full_reconstruction_l2_error,
                poisson_exact_phi_exact_gradient_linf_error,
                poisson_numerical_phi_exact_gradient_linf_error,
                poisson_exact_phi_numerical_gradient_linf_error,
                poisson_full_reconstruction_linf_error
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
                      << ", ms=" << ms;
            if(method==WallDistanceMethod::POISSON) {
                std::cout << ", wall_ray_hits=" << result.poisson_wall_ray_hits
                          << ", wall_ray_misses=" << result.poisson_wall_ray_misses
                          << ", wall_fallbacks=" << result.poisson_wall_fallbacks
                          << ", wall_bad_alignment=" << result.poisson_wall_bad_alignment
                          << ", wall_min_alignment=" << result.poisson_wall_min_alignment
                          << ", phi=[" << result.poisson_phi_min << "," << result.poisson_phi_max << "]"
                          << ", grad_phi=[" << result.poisson_grad_min << "," << result.poisson_grad_max << "]"
                          << ", poisson_distance_L2=" << result.poisson_distance_l2_error
                          << ", poisson_distance_Linf=" << result.poisson_distance_linf_error;
            }
            std::cout << "\n";
        }
    }

    std::ofstream csv(output);
    if (!csv) {
        std::cerr << "cannot open output: " << output << "\n";
        return 3;
    }

    csv << "N,h,method,l1_relative,l2_relative,linf_relative,"
           "near_wall_l2_relative,residual_inf,monotonicity_violations,"
           "iterations,converged,stopping_reason,time_ms,"
           "poisson_wall_ray_hits,poisson_wall_ray_misses,poisson_wall_fallbacks,"
           "poisson_wall_bad_alignment,poisson_wall_min_alignment,poisson_phi_min,"
           "poisson_phi_max,poisson_grad_min,poisson_grad_max,poisson_distance_l2_error,"
           "poisson_distance_linf_error,poisson_exact_gradient_l2_error,"
           "poisson_exact_gradient_linf_error,"
           "poisson_exact_phi_exact_gradient_l2_error,poisson_numerical_phi_exact_gradient_l2_error,"
           "poisson_exact_phi_numerical_gradient_l2_error,poisson_full_reconstruction_l2_error,"
           "poisson_exact_phi_exact_gradient_linf_error,poisson_numerical_phi_exact_gradient_linf_error,"
           "poisson_exact_phi_numerical_gradient_linf_error,poisson_full_reconstruction_linf_error\n";

    for (const auto& r : rows) {
        const double h = 2.0 / static_cast<double>(r.n - 1);
        csv << r.n << ',' << h << ',' << r.method << ','
            << std::setprecision(16)
            << r.l1 << ',' << r.l2 << ',' << r.linf << ','
            << r.near_l2 << ',' << r.residual << ','
            << r.violations << ',' << r.iterations << ','
            << (r.converged ? "true" : "false") << ',' << r.stopping_reason << ','
            << r.time_ms << ','
            << r.poisson_wall_ray_hits << ',' << r.poisson_wall_ray_misses << ','
            << r.poisson_wall_fallbacks << ',' << r.poisson_wall_bad_alignment << ','
            << r.poisson_wall_min_alignment << ',' << r.poisson_phi_min << ','
            << r.poisson_phi_max << ',' << r.poisson_grad_min << ','
            << r.poisson_grad_max << ',' << r.poisson_distance_l2_error << ','
            << r.poisson_distance_linf_error << ','
            << r.poisson_exact_gradient_l2_error << ','
            << r.poisson_exact_gradient_linf_error << ','
            << r.poisson_exact_phi_exact_gradient_l2_error << ','
            << r.poisson_numerical_phi_exact_gradient_l2_error << ','
            << r.poisson_exact_phi_numerical_gradient_l2_error << ','
            << r.poisson_full_reconstruction_l2_error << ','
            << r.poisson_exact_phi_exact_gradient_linf_error << ','
            << r.poisson_numerical_phi_exact_gradient_linf_error << ','
            << r.poisson_exact_phi_numerical_gradient_linf_error << ','
            << r.poisson_full_reconstruction_linf_error
            << '\n';
    }

    if(poisson_failure) return 4;
    if(advection_diffusion_failure) return 5;
    return 0;
}
