#include "cfdx/physics/wall_distance.h"

#include <algorithm>
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

// Controlled-box diagnostics intentionally keep the production Poisson operator unchanged.\n
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

WallDistanceVec3 analytic_box_distance_gradient(const WallDistanceVec3& p, bool* differentiable = nullptr)
{
    const WallDistanceVec3 q{
        std::clamp(p.x,0.0,1.0),
        std::clamp(p.y,0.0,1.0),
        std::clamp(p.z,0.0,1.0)};
    const WallDistanceVec3 v=p-q;
    const double outside_distance=wd_norm(v);
    if (outside_distance > 1e-14) {
        if (differentiable) *differentiable = true;
        return v*(1.0/outside_distance);
    }

    const double distances[6] = {p.x, 1.0-p.x, p.y, 1.0-p.y, p.z, 1.0-p.z};
    std::size_t best = 0;
    for (std::size_t i=1; i<6; ++i)
        if (distances[i] < distances[best]) best = i;
    double second = std::numeric_limits<double>::infinity();
    for (std::size_t i=0; i<6; ++i)
        if (i != best) second = std::min(second, distances[i]);
    const double scale = std::max(1.0, std::abs(distances[best]));
    const bool unique = (second - distances[best]) > 1e-12 * scale;
    if (differentiable) *differentiable = unique;
    if (!unique) return {0.0,0.0,0.0};
    switch (best) {
        case 0: return {1.0,0.0,0.0};
        case 1: return {-1.0,0.0,0.0};
        case 2: return {0.0,1.0,0.0};
        case 3: return {0.0,-1.0,0.0};
        case 4: return {0.0,0.0,1.0};
        default: return {0.0,0.0,-1.0};
    }
}

void print_poisson_stencil_microscope(const WallDistanceBvh& bvh,
                                     const std::vector<double>& phi,
                                     const WallDistanceGrid& g,
                                     std::size_t id)
{
    if (id >= g.points.size() || g.solid[id]) return;
    const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
    const auto p=g.points[id];
    double aP=0.0, lhs=0.0;
    std::cout << "POISSON_STENCIL_MICROSCOPE cell=" << id
              << " ijk=" << i << ":" << j << ":" << k
              << " point=" << p.x << ":" << p.y << ":" << p.z
              << " phi=" << phi[id] << " target_b=1\\n";
    auto face=[&](const char* name,std::size_t q,bool exists,double h) {
        if(!exists) {
            std::cout << "  face=" << name << " type=outer coefficient=0 contribution=0\\n";
            return;
        }
        if(g.solid[q]) {
            const auto wd=poisson_wall_offset_diagnostic(bvh,g.points[id],g.points[q],h);
            const double coeff=(wd.delta>0.0 && std::isfinite(wd.delta)) ? 1.0/(h*wd.delta) : 0.0;
            const double contribution=coeff*phi[id];
            aP+=coeff;
            lhs+=contribution;
            std::cout << "  face=" << name << " type=solid_cut neighbour=" << q
                      << " delta=" << wd.delta << " delta_over_h=" << wd.delta/h
                      << " coefficient=" << coeff
                      << " contribution=" << contribution
                      << " ray_hit=" << wd.ray_hit
                      << " fallback=" << wd.fallback
                      << " alignment=" << wd.alignment << "\\n";
        } else {
            const double coeff=1.0/(h*h);
            const double contribution=coeff*(phi[id]-phi[q]);
            aP+=coeff;
            lhs+=contribution;
            std::cout << "  face=" << name << " type=fluid neighbour=" << q
                      << " coefficient=" << coeff
                      << " contribution=" << contribution
                      << " phi_neighbour=" << phi[q] << "\\n";
        }
    };
    face("xm",i>0?g.index(i-1,j,k):0,i>0,g.spacing.x);
    face("xp",i+1<g.nx?g.index(i+1,j,k):0,i+1<g.nx,g.spacing.x);
    face("ym",j>0?g.index(i,j-1,k):0,j>0,g.spacing.y);
    face("yp",j+1<g.ny?g.index(i,j+1,k):0,j+1<g.ny,g.spacing.y);
    face("zm",k>0?g.index(i,j,k-1):0,k>0,g.spacing.z);
    face("zp",k+1<g.nz?g.index(i,j,k+1):0,k+1<g.nz,g.spacing.z);
    std::cout << "  stencil_summary aP=" << aP
              << " A_h_phi=" << lhs
              << " operator_residual=" << std::abs(lhs-1.0) << "\\n";
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
    std::size_t poisson_gradient_nondifferentiable_cells{};
    double poisson_exact_phi_operator_residual_inf{};
    double poisson_exact_phi_cut_residual_inf{};
    double poisson_exact_phi_face_residual_inf{};
    double poisson_exact_phi_edge_residual_inf{};
    double poisson_exact_phi_corner_residual_inf{};
    std::size_t poisson_exact_phi_cut_cells{};
    std::size_t poisson_exact_phi_face_cells{};
    std::size_t poisson_exact_phi_edge_cells{};
    std::size_t poisson_exact_phi_corner_cells{};
    std::size_t poisson_exact_phi_max_residual_cell{};
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
            double poisson_exact_phi_exact_gradient_l2_error=0.0;
            double poisson_numerical_phi_exact_gradient_l2_error=0.0;
            double poisson_exact_phi_numerical_gradient_l2_error=0.0;
            double poisson_full_reconstruction_l2_error=0.0;
            double poisson_exact_phi_exact_gradient_linf_error=0.0;
            double poisson_numerical_phi_exact_gradient_linf_error=0.0;
            double poisson_exact_phi_numerical_gradient_linf_error=0.0;
            double poisson_full_reconstruction_linf_error=0.0;
            std::size_t poisson_gradient_nondifferentiable_cells=0;
            double poisson_exact_phi_operator_residual_inf=0.0;
            double poisson_exact_phi_cut_residual_inf=0.0;
            double poisson_exact_phi_face_residual_inf=0.0;
            double poisson_exact_phi_edge_residual_inf=0.0;
            double poisson_exact_phi_corner_residual_inf=0.0;
            std::size_t poisson_exact_phi_cut_cells=0;
            std::size_t poisson_exact_phi_face_cells=0;
            std::size_t poisson_exact_phi_edge_cells=0;
            std::size_t poisson_exact_phi_corner_cells=0;
            std::size_t poisson_exact_phi_max_residual_cell=0;
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
                std::vector<double> exact_phi_field(grid.points.size(), 0.0);
                for(std::size_t id=0;id<grid.points.size();++id) {
                    if(grid.solid[id]) continue;
                    const double d=analytic[id];
                    exact_phi_field[id]=d-0.5*d*d;
                }
                for(std::size_t id=0;id<grid.points.size();++id) {
                    if(grid.solid[id]) continue;
                    const std::size_t k=id/(grid.nx*grid.ny), rem=id%(grid.nx*grid.ny), j=rem/grid.nx, i=rem%grid.nx;
                    const bool outer = (i==0 || i+1==grid.nx || j==0 || j+1==grid.ny || k==0 || k+1==grid.nz);
                    std::size_t solid_neighbours=0;
                    auto count_solid=[&](std::size_t q,bool exists){ if(exists && grid.solid[q]) ++solid_neighbours; };
                    count_solid(i>0?grid.index(i-1,j,k):0,i>0); count_solid(i+1<grid.nx?grid.index(i+1,j,k):0,i+1<grid.nx);
                    count_solid(j>0?grid.index(i,j-1,k):0,j>0); count_solid(j+1<grid.ny?grid.index(i,j+1,k):0,j+1<grid.ny);
                    count_solid(k>0?grid.index(i,j,k-1):0,k>0); count_solid(k+1<grid.nz?grid.index(i,j,k+1):0,k+1<grid.nz);
                    // The manufactured field is designed for the solid-wall Poisson
                    // equation. It does not satisfy the artificial outer Neumann
                    // boundary, so outer-boundary cells are excluded from the MMS
                    // consistency gate and classified separately from physical walls.
                    if (outer) continue;
                    const double op_residual=std::abs(-poisson_laplacian_at(bvh,exact_phi_field,grid,id)-1.0);
                    if(op_residual>poisson_exact_phi_operator_residual_inf) poisson_exact_phi_max_residual_cell=id;
                    poisson_exact_phi_operator_residual_inf=std::max(poisson_exact_phi_operator_residual_inf,op_residual);
                    if(solid_neighbours>0) ++poisson_exact_phi_cut_cells;
                    if(solid_neighbours==1) { ++poisson_exact_phi_face_cells; poisson_exact_phi_face_residual_inf=std::max(poisson_exact_phi_face_residual_inf,op_residual); }
                    if(solid_neighbours==2) { ++poisson_exact_phi_edge_cells; poisson_exact_phi_edge_residual_inf=std::max(poisson_exact_phi_edge_residual_inf,op_residual); }
                    if(solid_neighbours>=3) { ++poisson_exact_phi_corner_cells; poisson_exact_phi_corner_residual_inf=std::max(poisson_exact_phi_corner_residual_inf,op_residual); }
                }
                print_poisson_stencil_microscope(bvh, exact_phi_field, grid, poisson_exact_phi_max_residual_cell);
                std::cout << "Poisson exact-phi operator audit"
                          << ", max_abs_Aphi_minus_b=" << poisson_exact_phi_operator_residual_inf
                          << ", max_cut_cell=" << poisson_exact_phi_cut_residual_inf
                          << ", max_face_region=" << poisson_exact_phi_face_residual_inf
                          << ", max_edge_region=" << poisson_exact_phi_edge_residual_inf
                          << ", max_corner_region=" << poisson_exact_phi_corner_residual_inf
                          << ", cut_cells=" << poisson_exact_phi_cut_cells
                          << ", face_cells=" << poisson_exact_phi_face_cells
                          << ", edge_cells=" << poisson_exact_phi_edge_cells
                          << ", corner_cells=" << poisson_exact_phi_corner_cells
                          << ", max_cell=" << poisson_exact_phi_max_residual_cell
                          << ", max_point=" << grid.points[poisson_exact_phi_max_residual_cell].x
                          << ":" << grid.points[poisson_exact_phi_max_residual_cell].y
                          << ":" << grid.points[poisson_exact_phi_max_residual_cell].z << "\n";
                for(std::size_t id=0;id<grid.points.size();++id) {
                    if(grid.solid[id]) continue;
                    bool exact_grad_defined = false;
                    const auto exact_grad=analytic_box_distance_gradient(grid.points[id], &exact_grad_defined);
                    if (!exact_grad_defined) {
                        ++poisson_gradient_nondifferentiable_cells;
                        continue;
                    }
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
                poisson_full_reconstruction_linf_error,
                poisson_gradient_nondifferentiable_cells,
                poisson_exact_phi_operator_residual_inf,
                poisson_exact_phi_cut_residual_inf,
                poisson_exact_phi_face_residual_inf,
                poisson_exact_phi_edge_residual_inf,
                poisson_exact_phi_corner_residual_inf,
                poisson_exact_phi_cut_cells,
                poisson_exact_phi_face_cells,
                poisson_exact_phi_edge_cells,
                poisson_exact_phi_corner_cells,
                poisson_exact_phi_max_residual_cell
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
           "poisson_exact_phi_numerical_gradient_linf_error,poisson_full_reconstruction_linf_error,"
           "poisson_gradient_nondifferentiable_cells,poisson_exact_phi_operator_residual_inf,"
           "poisson_exact_phi_cut_residual_inf,poisson_exact_phi_face_residual_inf,"
           "poisson_exact_phi_edge_residual_inf,poisson_exact_phi_corner_residual_inf,"
           "poisson_exact_phi_cut_cells,poisson_exact_phi_face_cells,poisson_exact_phi_edge_cells,"
           "poisson_exact_phi_corner_cells,poisson_exact_phi_max_residual_cell\n";

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
            << r.poisson_full_reconstruction_linf_error << ','
            << r.poisson_gradient_nondifferentiable_cells << ','
            << r.poisson_exact_phi_operator_residual_inf << ','
            << r.poisson_exact_phi_cut_residual_inf << ','
            << r.poisson_exact_phi_face_residual_inf << ','
            << r.poisson_exact_phi_edge_residual_inf << ','
            << r.poisson_exact_phi_corner_residual_inf << ','
            << r.poisson_exact_phi_cut_cells << ','
            << r.poisson_exact_phi_face_cells << ','
            << r.poisson_exact_phi_edge_cells << ','
            << r.poisson_exact_phi_corner_cells << ','
            << r.poisson_exact_phi_max_residual_cell
            << '\n';
    }

    if(poisson_failure) return 4;
    if(advection_diffusion_failure) return 5;
    return 0;
}
