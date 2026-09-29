#include "cfdx/physics/wall_distance.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace cfdx::physics;

namespace {

void add_triangle(WallSurface& s,std::size_t a,std::size_t b,std::size_t c) {
    s.triangles.push_back({{a,b,c}});
}

void add_box(WallSurface& s,const WallDistanceVec3& lo,const WallDistanceVec3& hi) {
    const std::size_t b=s.points.size();
    s.points.insert(s.points.end(),{
        {lo.x,lo.y,lo.z},{hi.x,lo.y,lo.z},{hi.x,hi.y,lo.z},{lo.x,hi.y,lo.z},
        {lo.x,lo.y,hi.z},{hi.x,lo.y,hi.z},{hi.x,hi.y,hi.z},{lo.x,hi.y,hi.z}});
    const std::array<std::array<std::size_t,4>,6> q={{
        {{0,3,2,1}},{{4,5,6,7}},{{0,1,5,4}},{{1,2,6,5}},{{2,3,7,6}},{{3,0,4,7}}}};
    for(const auto& f:q) { add_triangle(s,b+f[0],b+f[1],b+f[2]); add_triangle(s,b+f[0],b+f[2],b+f[3]); }
}

void add_cylinder_x(WallSurface& s,double x0,double x1,double cy,double cz,double r,std::size_t n) {
    const std::size_t b=s.points.size();
    for(std::size_t i=0;i<n;++i) {
        const double a=2.0*3.14159265358979323846*static_cast<double>(i)/static_cast<double>(n);
        s.points.push_back({x0,cy+r*std::cos(a),cz+r*std::sin(a)});
        s.points.push_back({x1,cy+r*std::cos(a),cz+r*std::sin(a)});
    }
    const std::size_t c0=s.points.size(); s.points.push_back({x0,cy,cz});
    const std::size_t c1=s.points.size(); s.points.push_back({x1,cy,cz});
    for(std::size_t i=0;i<n;++i) {
        const std::size_t j=(i+1)%n;
        add_triangle(s,b+2*i,b+2*j,b+2*j+1); add_triangle(s,b+2*i,b+2*j+1,b+2*i+1);
        add_triangle(s,c0,b+2*j,b+2*i);
        add_triangle(s,c1,b+2*i+1,b+2*j+1);
    }
}

WallSurface make_complex_wing_body_tail() {
    WallSurface s;
    // Deliberately non-trivial multi-body-like surface: fuselage, swept-ish
    // wing slabs, tail planes and a vertical fin. The reference is the
    // triangulated surface distance, not an analytic approximation.
    add_cylinder_x(s,0.0,9.0,0.0,0.0,0.72,64);
    add_box(s,{3.0,-3.8,-0.10},{6.4,3.8,0.10});
    add_box(s,{7.0,-1.9,-0.08},{8.7,1.9,0.08});
    add_box(s,{7.2,-0.12,0.0},{8.8,0.12,1.45});
    add_box(s,{-0.15,-1.35,-0.55},{1.9,1.35,0.0});
    return s;
}

bool inside_complex(const WallDistanceVec3& p) {
    // Keep the computational solid exactly consistent with add_cylinder_x():
    // the surface is a finite cylinder x in [0,9] with circular radius 0.72.
    // The previous ellipsoidal test silently changed the solid topology near
    // both cylinder end caps, invalidating the complex-geometry benchmark.
    if(p.x>=0.0 && p.x<=9.0 &&
       p.y*p.y+p.z*p.z<=0.72*0.72) return true;
    if(p.x>=3.0 && p.x<=6.4 && std::abs(p.y)<=3.8 && std::abs(p.z)<=0.10) return true;
    if(p.x>=7.0 && p.x<=8.7 && std::abs(p.y)<=1.9 && std::abs(p.z)<=0.08) return true;
    if(p.x>=7.2 && p.x<=8.8 && std::abs(p.y)<=0.12 && p.z>=0.0 && p.z<=1.45) return true;
    if(p.x>=-0.15 && p.x<=1.9 && std::abs(p.y)<=1.35 && p.z>=-0.55 && p.z<=0.0) return true;
    return false;
}

struct Row { std::string method; double l2,linf,near_l2,violations,residual; std::size_t invalid,iterations; double ms,init_ms,poisson_ms; std::size_t poisson_iterations; double poisson_residual; bool converged; };

} // namespace

int main(int argc,char** argv) {
    std::string output = argc>1 ? argv[1] : "wall_distance_complex_benchmark.csv";
    const WallSurface surface=make_complex_wing_body_tail();
    const WallDistanceGrid grid=make_wall_distance_grid(
        40,30,22,{-1.5,-5.0,-2.2},{0.30,0.34,0.21},inside_complex);
    const auto reference=exact_reference(surface,grid);
    const double h=std::min({grid.spacing.x,grid.spacing.y,grid.spacing.z});
    const WallDistanceBvh poisson_bvh(surface);
    const auto poisson_operator_audit=audit_poisson_operator(poisson_bvh,grid);
    const auto poisson_offset_audit=audit_poisson_wall_offsets(poisson_bvh,grid);
    const std::size_t max_id=poisson_offset_audit.max_wall_coefficient_cell;
    const std::size_t max_k=max_id/(grid.nx*grid.ny);
    const std::size_t max_rem=max_id%(grid.nx*grid.ny);
    const std::size_t max_j=max_rem/grid.nx;
    const std::size_t max_i=max_rem%grid.nx;
    const auto max_p=grid.points[max_id];

    const std::array<WallDistanceMethod,9> methods={{
        WallDistanceMethod::EXACT_GEOMETRIC,
        WallDistanceMethod::SEARCH_BASED,
        WallDistanceMethod::MESH_WAVE,
        WallDistanceMethod::DIRECTIONAL_MESH_WAVE,
        WallDistanceMethod::POISSON,
        WallDistanceMethod::EIKONAL,
        WallDistanceMethod::HAMILTON_JACOBI,
        WallDistanceMethod::ADVECTION_DIFFUSION,
        WallDistanceMethod::HYBRID_POISSON_EIKONAL}};

    std::vector<Row> rows;
    std::cout << "CFDX wall-distance complex-geometry benchmark\n";
    std::cout << "surface_vertices=" << surface.points.size()
              << " surface_triangles=" << surface.triangles.size()
              << " samples=" << grid.points.size() << " h=" << h << "\n";
    std::cout << "poisson_operator,fluid_nodes=" << poisson_operator_audit.fluid_nodes
              << ",fluid_fluid_faces=" << poisson_operator_audit.fluid_fluid_faces
              << ",solid_faces=" << poisson_operator_audit.solid_faces
              << ",outer_faces=" << poisson_operator_audit.outer_faces
              << ",min_diagonal=" << poisson_operator_audit.min_diagonal
              << ",max_diagonal=" << poisson_operator_audit.max_diagonal
              << ",min_diagonal_dominance=" << poisson_operator_audit.min_diagonal_dominance
              << ",symmetry_error=" << poisson_operator_audit.symmetry_error
              << ",min_delta_over_h=" << poisson_offset_audit.min_delta_over_h
              << ",max_delta_over_h=" << poisson_offset_audit.max_delta_over_h
              << ",max_wall_coefficient=" << poisson_offset_audit.max_wall_coefficient
              << ",max_wall_coefficient_cell=" << poisson_offset_audit.max_wall_coefficient_cell
              << ",max_wall_coefficient_ijk=" << max_i << ":" << max_j << ":" << max_k
              << ",max_wall_coefficient_point=" << max_p.x << ":" << max_p.y << ":" << max_p.z
              << ",max_wall_coefficient_delta=" << poisson_offset_audit.max_wall_coefficient_delta
              << ",max_wall_coefficient_h=" << poisson_offset_audit.max_wall_coefficient_h
              << ",degenerate_wall_offsets=" << poisson_offset_audit.degenerate_count << "\n";
    constexpr std::size_t benchmark_iterations=500;
    std::cout << "method,l2_relative,linf_relative,near_wall_l2_relative,monotonicity_violations,invalid,iterations,residual_inf,time_ms,eikonal_init_ms,poisson_stage_ms,poisson_iterations,poisson_residual_inf,converged\n";

    for(const auto method:methods) {
        double init_ms=0.0, poisson_ms=0.0;
        if(method==WallDistanceMethod::HAMILTON_JACOBI) {
            const auto p0=std::chrono::steady_clock::now();
            std::size_t init_iter=0;
            (void)eikonal_fast_sweep(surface,grid,benchmark_iterations,&init_iter);
            const auto p1=std::chrono::steady_clock::now();
            init_ms=std::chrono::duration<double,std::milli>(p1-p0).count();
        } else if(method==WallDistanceMethod::HYBRID_POISSON_EIKONAL) {
            const auto p0=std::chrono::steady_clock::now();
            std::size_t poisson_iter=0; double poisson_residual=0.0;
            (void)poisson_distance(surface,grid,benchmark_iterations,1.5,&poisson_iter,&poisson_residual);
            const auto p1=std::chrono::steady_clock::now();
            poisson_ms=std::chrono::duration<double,std::milli>(p1-p0).count();
        }
        const auto t0=std::chrono::steady_clock::now();
        const auto result=compute_wall_distance(method,surface,grid,benchmark_iterations);
        const auto t1=std::chrono::steady_clock::now();
        const double ms=std::chrono::duration<double,std::milli>(t1-t0).count();
        const auto m=compare_wall_distance(grid,reference,result.distance,2.0*h);
        std::size_t invalid=0; for(std::size_t i=0;i<result.distance.size();++i) if(!result.valid[i]) ++invalid;
        rows.push_back({result.method,m.l2_relative,m.linf_relative,m.near_wall_l2_relative,m.monotonicity_violations,result.residual_inf,invalid,result.iterations,ms,init_ms,poisson_ms,result.auxiliary_iterations,result.auxiliary_residual_inf,result.converged});
        std::cout << result.method << "," << std::setprecision(8)
                  << m.l2_relative << "," << m.linf_relative << ","
                  << m.near_wall_l2_relative << "," << m.monotonicity_violations << ","
                  << invalid << "," << result.iterations << "," << result.residual_inf << ","
                  << ms << "," << result.auxiliary_iterations << "," << result.auxiliary_residual_inf
                  << "," << (result.converged ? "true" : "false") << "\n";
    }

    std::ofstream csv(output);
    if(!csv) throw std::runtime_error("cannot open benchmark output: "+output);
    csv << "method,l2_relative,linf_relative,near_wall_l2_relative,monotonicity_violations,invalid,iterations,residual_inf,time_ms,eikonal_init_ms,poisson_stage_ms,poisson_iterations,poisson_residual_inf,converged\n";
    for(const auto& r:rows)
        csv << r.method << "," << r.l2 << "," << r.linf << "," << r.near_l2 << ","
            << r.violations << "," << r.invalid << "," << r.iterations << "," << r.residual << "," << r.ms
            << "," << r.init_ms << "," << r.poisson_ms << "," << r.poisson_iterations
            << "," << r.poisson_residual << "," << (r.converged ? "true" : "false") << "\n";
    csv.close();

    // The exact method must be an exact self-reference. This is a regression
    // gate for the benchmark infrastructure, not a physical validation claim.
    if(rows.front().l2>1e-14 || rows.front().near_l2>1e-14)
        throw std::runtime_error("exact wall-distance reference is not self-consistent");
    return 0;
}
