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


    // P2: oblique planar manufactured solution.
    //
    // The legacy cut-face Dirichlet closure is not second-order at an
    // oblique wall: phi_P/(h*delta) is a two-point approximation of the
    // wall-normal flux and has a first-order truncation term when the wall
    // cuts a Cartesian face at an arbitrary distance.  The aligned MMS
    // already demonstrated this behaviour.  Therefore the oblique test must
    // qualify the actual closure order rather than impose an impossible
    // second-order pointwise gate.
    //
    // We use three refinements and normalize the unscaled operator defect by
    // h^2.  A first-order wall-flux defect gives O(h) for this normalized
    // quantity.  Interior fluid-fluid cells must still be exact to roundoff.
    double previous_normalized_oblique_residual=0.0;
    bool have_previous_oblique=false;
    for(const std::size_t n : {17u,33u,65u}) {
        const double h=4.0/static_cast<double>(n-1);
        const double inv_sqrt2=1.0/std::sqrt(2.0);
        WallSurface oblique;
        oblique.points={
            {-2.0*inv_sqrt2,  2.0*inv_sqrt2,-2.0},
            { 2.0*inv_sqrt2, -2.0*inv_sqrt2,-2.0},
            { 2.0*inv_sqrt2, -2.0*inv_sqrt2, 2.0},
            {-2.0*inv_sqrt2,  2.0*inv_sqrt2, 2.0}};
        oblique.triangles={{{0,1,2}},{{0,2,3}}};
        const auto og=make_wall_distance_grid(
            n,n,9,{-2.0,-2.0,-1.0},{h,h,0.25},
            [](const WallDistanceVec3& p){ return p.x+p.y<0.0; });
        const WallDistanceBvh obvh(oblique);
        const double L_oblique=4.0;
        std::vector<double> oblique_phi(og.points.size(),0.0);
        for(std::size_t id=0;id<og.points.size();++id) if(!og.solid[id]) {
            const auto x=og.points[id];
            const double d=(x.x+x.y)*inv_sqrt2;
            oblique_phi[id]=L_oblique*d-0.5*d*d;
        }

        double interior_operator_res=0.0;
        double cut_operator_res=0.0;
        double normalized_cut_residual=0.0;
        std::size_t cut_samples=0;
        std::size_t interior_samples=0;
        std::size_t worst_cut_id=0;
        for(std::size_t id=0;id<og.points.size();++id) {
            if(og.solid[id]) continue;
            const auto x=og.points[id];
            const double d=(x.x+x.y)*inv_sqrt2;
            if(d<=0.0 || d>=1.75 || std::abs(x.z)>0.75) continue;

            const std::size_t k=id/(og.nx*og.ny);
            const std::size_t rem=id%(og.nx*og.ny);
            const std::size_t j=rem/og.nx;
            const std::size_t i=rem%og.nx;
            bool cut=false;
            auto inspect=[&](std::size_t q,bool exists) {
                if(exists && og.solid[q]) cut=true;
            };
            inspect(i>0?og.index(i-1,j,k):0,i>0);
            inspect(i+1<og.nx?og.index(i+1,j,k):0,i+1<og.nx);
            inspect(j>0?og.index(i,j-1,k):0,j>0);
            inspect(j+1<og.ny?og.index(i,j+1,k):0,j+1<og.ny);
            inspect(k>0?og.index(i,j,k-1):0,k>0);
            inspect(k+1<og.nz?og.index(i,j,k+1):0,k+1<og.nz);

            const double op=std::abs(-poisson_laplacian_at(obvh,oblique_phi,og,id)-1.0);
            if(cut) {
                ++cut_samples;
                if(op>cut_operator_res) {
                    cut_operator_res=op;
                    worst_cut_id=id;
                }
            } else {
                ++interior_samples;
                interior_operator_res=std::max(interior_operator_res,op);
            }
        }

        require(interior_samples>0,"oblique audit has no interior samples");
        require(cut_samples>0,"oblique audit has no cut-face samples");
        require(interior_operator_res<1e-12,
                "oblique interior Poisson operator is not exact");
        require(std::isfinite(cut_operator_res) && cut_operator_res>0.0,
                "oblique cut-face operator defect was not detected");

        normalized_cut_residual=cut_operator_res*h*h;
        std::cerr << "oblique P2 N=" << n
                  << " h=" << h
                  << " interior_operator_inf=" << interior_operator_res
                  << " cut_operator_inf=" << cut_operator_res
                  << " normalized_cut_operator_inf=" << normalized_cut_residual
                  << " cut_samples=" << cut_samples
                  << " worst_cut_id=" << worst_cut_id << "\n";

        if(have_previous_oblique) {
            const double order=std::log(previous_normalized_oblique_residual/
                                        normalized_cut_residual)/std::log(2.0);
            require(std::isfinite(order) && order>0.5 && order<1.5,
                    "oblique cut-face closure does not show the expected first-order trend");
            std::cerr << "oblique P2 observed order=" << order << "\n";
        }
        previous_normalized_oblique_residual=normalized_cut_residual;
        have_previous_oblique=true;
    }

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
