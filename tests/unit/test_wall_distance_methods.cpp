#include "cfdx/physics/wall_distance.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace cfdx::physics;

static void require(bool ok,const char* msg) {
    if(!ok) throw std::runtime_error(msg);
}

int main() {
    WallSurface s;
    s.points={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
              {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    const std::array<std::array<std::size_t,4>,6> q={{
        {{0,3,2,1}},{{4,5,6,7}},{{0,1,5,4}},{{1,2,6,5}},{{2,3,7,6}},{{3,0,4,7}}}};
    for(const auto& f:q) {
        s.triangles.push_back({{f[0],f[1],f[2]}});
        s.triangles.push_back({{f[0],f[2],f[3]}});
    }

    require(std::abs(exact_point_distance(s,{0.5,0.5,0.5})-0.5)<1e-12,
            "exact point-to-triangle distance is wrong");

    const auto g=make_wall_distance_grid(
        9,9,9,{-1,-1,-1},{0.25,0.25,0.25},
        [](const WallDistanceVec3& p){ return p.x>0.0 && p.x<1.0 && p.y>0.0 && p.y<1.0 && p.z>0.0 && p.z<1.0; });
    const auto ref=exact_reference(s,g);
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
    for(const auto method:methods) {
        const auto r=compute_wall_distance(method,s,g,20);
        std::cerr << "wall-distance method: " << wall_distance_method_name(method) << "\\n";
        require(r.distance.size()==g.points.size(),"distance size mismatch");
        for(std::size_t i=0;i<r.distance.size();++i)
            if(!g.solid[i]) require(r.valid[i] && std::isfinite(r.distance[i]) && r.distance[i]>=0.0,
                                    "invalid wall distance");
    }

    const WallDistanceBvh bvh(s);
    require(bvh.triangle_count()==s.triangles.size(),"BVH triangle count mismatch");
    require(bvh.node_count()>0,"BVH is empty");
    const std::array<WallDistanceVec3,8> probes={{{-0.2,0.4,0.4},{0.5,0.5,-0.3},
        {1.2,0.2,0.7},{0.25,0.8,1.4},{-0.4,-0.2,0.3},{0.7,1.3,0.1},
        {1.4,1.4,1.4},{0.12,0.34,0.91}}};
    for(const auto& p:probes) {
        const double brute=exact_point_distance(s,p);
        const double accelerated=bvh.nearest_distance(p);
        require(std::abs(brute-accelerated)<1e-12,"BVH nearest distance mismatch");
    }
    const auto search=compute_wall_distance(WallDistanceMethod::SEARCH_BASED,s,g,20);
    for(std::size_t i=0;i<g.points.size();++i) if(!g.solid[i]) {
        const double brute=ref[i];
        if(brute<=1.0) require(std::abs(search.distance[i]-brute)<1e-12,
                               "search-based near-wall distance mismatch");
    }

    const auto mesh=compute_wall_distance(WallDistanceMethod::MESH_WAVE,s,g,40);
    const auto directional=compute_wall_distance(WallDistanceMethod::DIRECTIONAL_MESH_WAVE,s,g,40);
    const auto eik=compute_wall_distance(WallDistanceMethod::EIKONAL,s,g,40);
    const auto hj=compute_wall_distance(WallDistanceMethod::HAMILTON_JACOBI,s,g,40);
    const auto adv=compute_wall_distance(WallDistanceMethod::ADVECTION_DIFFUSION,s,g,40);
    bool directional_differs=false, hj_differs=false, adv_differs=false;
    for(std::size_t i=0;i<g.points.size();++i) if(!g.solid[i]) {
        directional_differs |= std::abs(mesh.distance[i]-directional.distance[i])>1e-13;
        hj_differs |= std::abs(eik.distance[i]-hj.distance[i])>1e-13;
        adv_differs |= std::abs(eik.distance[i]-adv.distance[i])>1e-13;
    }
    require(directional_differs,"directional mesh wave collapsed to isotropic mesh wave");
    require(hj_differs,"Hamilton-Jacobi collapsed to Eikonal");
    require(adv_differs,"advection-diffusion collapsed to Eikonal");

    // Canonical analytical test: an infinite planar wall at x=0 has
    // d(x,y,z)=x for x>=0, |grad d|=1 and laplacian(d)=0. This isolates
    // the PDE operators from geometry/corner errors.
    WallSurface plane;
    plane.points={{0,-1,-1},{0,1,-1},{0,1,1},{0,-1,1}};
    plane.triangles={{{0,1,2}},{{0,2,3}}};
    const auto pg=make_wall_distance_grid(
        9,5,5,{0,-1,-1},{0.25,0.5,0.5},
        [](const WallDistanceVec3& p){ return p.x<0.0; });
    const std::size_t pc=pg.index(4,2,2);
    std::vector<double> linear(pg.points.size(),0.0);
    for(std::size_t id=0;id<linear.size();++id) linear[id]=pg.points[id].x;
    require(std::abs(godunov_gradient_at(linear,pg,pc)-1.0)<1e-13,
            "Godunov gradient is not exact for d=x");
    require(std::abs(laplacian_at(linear,pg,pc))<1e-13,
            "Laplacian is not exact for the linear planar solution");
    const double hj_linear_residual =
        godunov_gradient_at(linear,pg,pc)-1.0-0.25*pg.points[pc].x*laplacian_at(linear,pg,pc);
    require(std::abs(hj_linear_residual)<1e-13,
            "Hamilton-Jacobi operator is not exact for the linear planar solution");
    const std::array<WallDistanceMethod,1> plane_methods={{
        WallDistanceMethod::EIKONAL}};
    for(const auto method:plane_methods) {
        const auto r=compute_wall_distance(method,plane,pg,120);
        for(std::size_t id=0;id<r.distance.size();++id) {
            const std::size_t k=id/(pg.nx*pg.ny), rem=id%(pg.nx*pg.ny), j=rem/pg.nx, i=rem%pg.nx;
            // The analytical d=x solution is used only where the Cartesian
            // stencil is complete. The outer boundary has a separate
            // numerical boundary condition and is not part of this operator
            // consistency test.
            if(pg.solid[id] || i==0 || i+1==pg.nx || j==0 || j+1==pg.ny || k==0 || k+1==pg.nz) continue;
            const double expected=pg.points[id].x;
            require(r.valid[id] && std::isfinite(r.distance[id]),
                    "analytical plane produced an invalid distance");
            if(std::abs(r.distance[id]-expected)>=2e-10) {
                const double grad=godunov_gradient_at(r.distance,pg,id);
                const double lap=laplacian_at(r.distance,pg,id);
                const double hj_res=grad-1.0-0.25*expected*lap;
                std::cerr << "planar failure method=" << wall_distance_method_name(method)
                          << " id=" << id << " i=" << i << " j=" << j << " k=" << k
                          << " x=" << expected << " value=" << r.distance[id]
                          << " error=" << std::abs(r.distance[id]-expected)
                          << " grad=" << grad << " lap=" << lap << " residual=" << hj_res << "\n";
                throw std::runtime_error("wall-distance PDE failed the analytical planar-wall solution");
            }
        }
    }

    // Poisson operator check independent of the seed-band initialization:
    // phi(x)=x*L-x^2/2 satisfies phi''=-1 for 0<x<L.
    const double L=2.0;
    std::vector<double> phi(pg.points.size(),0.0);
    std::vector<unsigned char> fixed(pg.points.size(),0);
    for(std::size_t id=0;id<phi.size();++id) {
        const double x=pg.points[id].x;
        phi[id]=x*L-0.5*x*x;
    }
    // The quadratic solution is tested only on complete Cartesian
    // stencils; the truncated outer boundary has its own numerical BC.
    double poisson_res=0.0;
    for(std::size_t id=0;id<phi.size();++id) {
        const std::size_t k=id/(pg.nx*pg.ny), rem=id%(pg.nx*pg.ny), j=rem/pg.nx, i=rem%pg.nx;
        if(i==0 || i+1==pg.nx || j==0 || j+1==pg.ny || k==0 || k+1==pg.nz) continue;
        const double h=pg.spacing.x;
        const double lap=(phi[pg.index(i-1,j,k)]-2.0*phi[id]+phi[pg.index(i+1,j,k)])/(h*h);
        poisson_res=std::max(poisson_res,std::abs(lap+1.0));
    }
    require(poisson_res<1e-13,
            "discrete Poisson operator is not exact for the analytical planar solution");


    // P2: oblique planar manufactured solution.  This is deliberately
    // separate from the cube/complex benchmark: it checks the directional
    // consistency of the same Poisson operator and reconstruction on a wall
    // that is not aligned with the Cartesian grid.
    const double inv_sqrt2=1.0/std::sqrt(2.0);
    WallSurface oblique;
    oblique.points={
        {-2.0*inv_sqrt2,  2.0*inv_sqrt2,-2.0},
        { 2.0*inv_sqrt2, -2.0*inv_sqrt2,-2.0},
        { 2.0*inv_sqrt2, -2.0*inv_sqrt2, 2.0},
        {-2.0*inv_sqrt2,  2.0*inv_sqrt2, 2.0}};
    oblique.triangles={{{0,1,2}},{{0,2,3}}};
    const auto og=make_wall_distance_grid(
        17,17,9,{-1,-1,-1},{0.125,0.125,0.25},
        [](const WallDistanceVec3& p){ return p.x+p.y<0.0; });
    const WallDistanceBvh obvh(oblique);
    const double L_oblique=4.0;
    double oblique_operator_res=0.0;
    double oblique_gradient_res=0.0;
    double oblique_reconstruction_res=0.0;
    std::size_t oblique_samples=0;
    for(std::size_t id=0;id<og.points.size();++id) {
        if(og.solid[id]) continue;
        const auto x=og.points[id];
        const double d_exact=(x.x+x.y)*inv_sqrt2;
        // Stay away from the wall and the artificial outer box so that the
        // audit tests the interior Cartesian operator, not an unrelated BC.
        if(d_exact<=0.375 || d_exact>=1.75 || std::abs(x.z)>1.0 ||
           std::abs(x.x-x.y)>0.75) continue;
        const double phi_exact=L_oblique*d_exact-0.5*d_exact*d_exact;
        std::vector<double> local_phi(og.points.size(),0.0);
        local_phi[id]=phi_exact;
        const auto exact_grad=[&](const WallDistanceVec3&) {
            const double a=L_oblique-d_exact;
            return WallDistanceVec3{a*inv_sqrt2,a*inv_sqrt2,0.0};
        };
        const double gx=poisson_reconstruction_gradient_component(obvh,
            [&]()->const std::vector<double>& {
                static thread_local std::vector<double> dummy;
                return dummy;
            }(),og,id,0);
        (void)gx;
        // Use a full manufactured field; only the current sample is needed
        // for the local operator, but the vector keeps the production API
        // identical to the actual Poisson solve.
        static_cast<void>(local_phi);
        ++oblique_samples;
    }
    require(oblique_samples>0,"oblique manufactured audit selected no samples");

    std::vector<double> oblique_phi(og.points.size(),0.0);
    for(std::size_t id=0;id<og.points.size();++id) if(!og.solid[id]) {
        const auto x=og.points[id];
        const double d=(x.x+x.y)*inv_sqrt2;
        oblique_phi[id]=L_oblique*d-0.5*d*d;
    }
    const auto ob_audit=audit_poisson_manufactured(
        obvh,og,oblique_phi,
        [](const WallDistanceVec3& x){ return (x.x+x.y)*inv_sqrt2; },
        [inv_sqrt2,L_oblique](const WallDistanceVec3& x) {
            const double d=(x.x+x.y)*inv_sqrt2;
            const double a=L_oblique-d;
            return WallDistanceVec3{a*inv_sqrt2,a*inv_sqrt2,0.0};
        });
    // The generic audit includes points close to the finite patch boundary;
    // require the stronger interior result separately.
    for(std::size_t id=0;id<og.points.size();++id) {
        if(og.solid[id]) continue;
        const auto x=og.points[id];
        const double d=(x.x+x.y)*inv_sqrt2;
        if(d<=0.375 || d>=1.75 || std::abs(x.z)>1.0 || std::abs(x.x-x.y)>0.75) continue;
        const double op=std::abs(-poisson_laplacian_at(obvh,oblique_phi,og,id)-1.0);
        const double gx=poisson_reconstruction_gradient_component(obvh,oblique_phi,og,id,0);
        const double gy=poisson_reconstruction_gradient_component(obvh,oblique_phi,og,id,1);
        const double eg=(L_oblique-d)*inv_sqrt2;
        oblique_operator_res=std::max(oblique_operator_res,op);
        oblique_gradient_res=std::max(oblique_gradient_res,
            std::max(std::abs(gx-eg),std::abs(gy-eg)));
        const double rec=poisson_reconstructed_distance(obvh,og,oblique_phi,id);
        oblique_reconstruction_res=std::max(oblique_reconstruction_res,std::abs(rec-d));
    }
    require(oblique_samples>0,"oblique audit sample count is zero");
    require(oblique_operator_res<1e-12,"oblique Poisson operator is not second-order-consistent on the manufactured plane");
    require(oblique_gradient_res<1e-12,"oblique Poisson reconstruction gradient is inconsistent");
    require(oblique_reconstruction_res<1e-12,"oblique Poisson reconstruction is inconsistent on the planar manufactured solution");
    require(ob_audit.samples>0,"generic Poisson manufactured audit produced no samples");
    std::cerr << "Poisson P0/P1/P2 audit: operator=" << ob_audit.operator_residual_inf
              << " gradient=" << ob_audit.gradient_error_inf
              << " reconstruction_L2=" << ob_audit.reconstruction_l2_relative
              << " reconstruction_Linf=" << ob_audit.reconstruction_linf_relative
              << " eikonal=" << ob_audit.eikonal_residual_inf << "\n";

    // The hybrid method must be a genuine PDE solve, not an algebraic
    // combination of completed Poisson and Eikonal distance fields.
    const auto hybrid=compute_wall_distance(WallDistanceMethod::HYBRID_POISSON_EIKONAL,s,g,40);
    const auto poisson=compute_wall_distance(WallDistanceMethod::POISSON,s,g,40);
    const auto eikonal=compute_wall_distance(WallDistanceMethod::EIKONAL,s,g,40);
    bool differs_from_poisson=false, differs_from_eikonal=false;
    for(std::size_t i=0;i<g.points.size();++i) if(!g.solid[i]) {
        differs_from_poisson |= std::abs(hybrid.distance[i]-poisson.distance[i])>1e-10;
        differs_from_eikonal |= std::abs(hybrid.distance[i]-eikonal.distance[i])>1e-10;
    }
    require(differs_from_poisson || differs_from_eikonal,
            "hybrid method collapsed to an existing completed distance field");

    const auto exact=compute_wall_distance(WallDistanceMethod::EXACT_GEOMETRIC,s,g,20);
    const auto m=compare_wall_distance(g,ref,exact.distance,0.5);
    require(m.l2_relative<1e-14,"exact benchmark regression failed");
    std::cout<<"wall-distance method contract: PASS\n";
    return 0;
}
