#include "cfdx/physics/wall_distance.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace cfdx::physics;

namespace {

struct Metrics {
    double phi_l2{};
    double phi_linf{};
    double operator_exact_inf{};
    double residual_inf{};
    double grad_exact_l2{};
    double grad_exact_linf{};
    double recon_exact_l2{};
    double recon_phi_exact_grad_l2{};
    double recon_phi_num_grad_exact_l2{};
    double recon_full_l2{};
    double recon_exact_linf{};
    double recon_phi_exact_grad_linf{};
    double recon_phi_num_grad_exact_linf{};
    double recon_full_linf{};
    std::size_t iterations{};
    bool converged{};
};

double rel_l2(const std::vector<double>& a,const std::vector<double>& b,
              const std::vector<unsigned char>& solid) {
    double e2=0.0,r2=0.0;
    for(std::size_t i=0;i<a.size();++i) {
        if(solid[i] || !std::isfinite(a[i]) || !std::isfinite(b[i])) continue;
        const double e=a[i]-b[i];
        e2+=e*e; r2+=b[i]*b[i];
    }
    return std::sqrt(e2/std::max(1e-30,r2));
}

double rel_linf(const std::vector<double>& a,const std::vector<double>& b,
                const std::vector<unsigned char>& solid) {
    double em=0.0,rm=0.0;
    for(std::size_t i=0;i<a.size();++i) {
        if(solid[i] || !std::isfinite(a[i]) || !std::isfinite(b[i])) continue;
        em=std::max(em,std::abs(a[i]-b[i]));
        rm=std::max(rm,std::abs(b[i]));
    }
    return em/std::max(1e-30,rm);
}

double relative_distance_error(const std::vector<double>& d,
                               const std::vector<double>& ref,
                               const std::vector<unsigned char>& solid,
                               bool linf) {
    double e2=0.0,r2=0.0,em=0.0,rm=0.0;
    for(std::size_t i=0;i<d.size();++i) {
        if(solid[i] || !std::isfinite(d[i]) || !std::isfinite(ref[i])) continue;
        const double e=std::abs(d[i]-ref[i]);
        e2+=e*e; r2+=ref[i]*ref[i]; em=std::max(em,e); rm=std::max(rm,std::abs(ref[i]));
    }
    return linf ? em/std::max(1e-30,rm) : std::sqrt(e2/std::max(1e-30,r2));
}

double reconstruct(double phi,double gx,double gy,double gz) {
    const double g=std::sqrt(gx*gx+gy*gy+gz*gz);
    return std::max(0.0,std::sqrt(std::max(0.0,g*g+2.0*phi))-g);
}

} // namespace

int main(int argc,char** argv) {
    const std::string output=argc>1 ? argv[1] : "wall_distance_poisson_decomposition.csv";
    const std::vector<std::size_t> sizes={16,32,64,128};
    std::ofstream csv(output);
    if(!csv) throw std::runtime_error("cannot open output: "+output);

    csv << "N,h,phi_l2_relative,phi_linf_relative,operator_exact_inf,residual_inf,"
           "grad_exact_l2_relative,grad_exact_linf_relative,"
           "recon_exact_l2_relative,recon_phi_num_grad_exact_l2_relative,"
           "recon_phi_exact_grad_num_l2_relative,recon_full_l2_relative,"
           "recon_exact_linf_relative,recon_phi_num_grad_exact_linf_relative,"
           "recon_phi_exact_grad_num_linf_relative,recon_full_linf_relative,"
           "iterations,converged\n";

    for(const auto n:sizes) {
        const double L=1.0;
        const double h=L/static_cast<double>(n);
        WallSurface surface;
        surface.points={{-2.0,0.0,-2.0},{2.0,0.0,-2.0},{2.0,0.0,2.0},{-2.0,0.0,2.0}};
        surface.triangles={{{0,1,2}},{{0,2,3}}};

        const auto grid=make_wall_distance_grid(
            n,n+1,n,{-2.0,-0.5*h,-2.0},
            {4.0/static_cast<double>(n-1),h,4.0/static_cast<double>(n-1)},
            [](const WallDistanceVec3& p){ return p.y<=0.0; });
        const WallDistanceBvh bvh(surface);

        std::size_t iterations=0;
        double residual=0.0;
        const auto phi=poisson_potential(bvh,grid,5000,1.0,&iterations,&residual);

        std::vector<double> phi_exact(phi.size(),std::numeric_limits<double>::infinity());
        std::vector<double> d_exact(phi.size(),std::numeric_limits<double>::infinity());
        std::vector<double> d_exact_inputs(phi.size(),std::numeric_limits<double>::infinity());
        std::vector<double> d_num_phi_exact_grad(phi.size(),std::numeric_limits<double>::infinity());
        std::vector<double> d_exact_phi_num_grad(phi.size(),std::numeric_limits<double>::infinity());
        std::vector<double> d_full(phi.size(),std::numeric_limits<double>::infinity());

        Metrics m;
        for(std::size_t id=0;id<grid.points.size();++id) {
            if(grid.solid[id]) continue;
            const double y=grid.points[id].y;
            phi_exact[id]=L*y-0.5*y*y;
            d_exact[id]=y;

            const double exact_gx=0.0, exact_gy=L-y, exact_gz=0.0;
            const double gex=std::sqrt(exact_gx*exact_gx+exact_gy*exact_gy+exact_gz*exact_gz);
            d_exact_inputs[id]=reconstruct(phi_exact[id],exact_gx,exact_gy,exact_gz);

            const double gx_phi=poisson_reconstruction_gradient_component(
                bvh,phi,g,id,0);
            const double gy_phi=poisson_reconstruction_gradient_component(
                bvh,phi,g,id,1);
            const double gz_phi=poisson_reconstruction_gradient_component(
                bvh,phi,g,id,2);

            const double gx_exact=poisson_reconstruction_gradient_component(
                bvh,phi_exact,g,id,0);
            const double gy_exact=poisson_reconstruction_gradient_component(
                bvh,phi_exact,g,id,1);
            const double gz_exact=poisson_reconstruction_gradient_component(
                bvh,phi_exact,g,id,2);

            d_num_phi_exact_grad[id]=reconstruct(phi[id],exact_gx,exact_gy,exact_gz);
            d_exact_phi_num_grad[id]=reconstruct(phi_exact[id],gx_phi,gy_phi,gz_phi);
            d_full[id]=reconstruct(phi[id],gx_phi,gy_phi,gz_phi);

            const double ge=std::sqrt(gx_exact*gx_exact+gy_exact*gy_exact+gz_exact*gz_exact);
            const double ephi=std::abs(phi[id]-phi_exact[id]);
            m.operator_exact_inf=std::max(
                m.operator_exact_inf,
                std::abs(-poisson_laplacian_at(bvh,phi_exact,grid,id)-1.0));
            m.grad_exact_l2=std::max(m.grad_exact_l2,std::abs(ge-gex));
            m.grad_exact_linf=std::max(m.grad_exact_linf,
                std::max({std::abs(gx_exact-exact_gx),
                          std::abs(gy_exact-exact_gy),
                          std::abs(gz_exact-exact_gz)}));
            m.phi_l2=std::max(m.phi_l2,ephi);
            m.phi_linf=std::max(m.phi_linf,ephi);
        }

        m.phi_l2=rel_l2(phi,phi_exact,grid.solid);
        m.phi_linf=rel_linf(phi,phi_exact,grid.solid);
        m.residual_inf=poisson_residual_inf(bvh,phi,grid,grid.solid);
        m.grad_exact_l2=0.0;
        m.grad_exact_linf=0.0;

        double ge2=0.0, gr2=0.0, gem=0.0, grm=0.0;
        for(std::size_t id=0;id<grid.points.size();++id) {
            if(grid.solid[id]) continue;
            const double y=grid.points[id].y;
            const double gx=poisson_reconstruction_gradient_component(bvh,phi_exact,grid,id,0);
            const double gy=poisson_reconstruction_gradient_component(bvh,phi_exact,grid,id,1);
            const double gz=poisson_reconstruction_gradient_component(bvh,phi_exact,grid,id,2);
            const double ex=0.0, ey=L-y, ez=0.0;
            const double e2=(gx-ex)*(gx-ex)+(gy-ey)*(gy-ey)+(gz-ez)*(gz-ez);
            ge2+=e2; gr2+=ey*ey; gem=std::max(gem,std::sqrt(e2)); grm=std::max(grm,std::abs(ey));
        }
        m.grad_exact_l2=std::sqrt(ge2/std::max(1e-30,gr2));
        m.grad_exact_linf=gem/std::max(1e-30,grm);

        m.recon_exact_l2=relative_distance_error(d_exact_inputs,d_exact,grid.solid,false);
        m.recon_phi_num_grad_exact_l2=relative_distance_error(d_num_phi_exact_grad,d_exact,grid.solid,false);
        m.recon_phi_exact_grad_num_l2=relative_distance_error(d_exact_phi_num_grad,d_exact,grid.solid,false);
        m.recon_full_l2=relative_distance_error(d_full,d_exact,grid.solid,false);
        m.recon_exact_linf=relative_distance_error(d_exact_inputs,d_exact,grid.solid,true);
        m.recon_phi_num_grad_exact_linf=relative_distance_error(d_num_phi_exact_grad,d_exact,grid.solid,true);
        m.recon_phi_exact_grad_num_linf=relative_distance_error(d_exact_phi_num_grad,d_exact,grid.solid,true);
        m.recon_full_linf=relative_distance_error(d_full,d_exact,grid.solid,true);
        m.iterations=iterations;
        m.converged=residual<1e-8;

        csv << n << ',' << std::setprecision(16) << h << ','
            << m.phi_l2 << ',' << m.phi_linf << ',' << m.operator_exact_inf << ','
            << m.residual_inf << ',' << m.grad_exact_l2 << ',' << m.grad_exact_linf << ','
            << m.recon_exact_l2 << ',' << m.recon_phi_num_grad_exact_l2 << ','
            << m.recon_phi_exact_grad_num_l2 << ',' << m.recon_full_l2 << ','
            << m.recon_exact_linf << ',' << m.recon_phi_num_grad_exact_linf << ','
            << m.recon_phi_exact_grad_num_linf << ',' << m.recon_full_linf << ','
            << m.iterations << ',' << (m.converged ? "true" : "false") << '\n';

        std::cout << "N=" << n
                  << " phi_L2=" << m.phi_l2
                  << " operator_exact_inf=" << m.operator_exact_inf
                  << " residual=" << m.residual_inf
                  << " grad_exact_L2=" << m.grad_exact_l2
                  << " recon_exact=" << m.recon_exact_l2
                  << " recon_phi_num_grad_exact=" << m.recon_phi_num_grad_exact_l2
                  << " recon_phi_exact_grad_num=" << m.recon_phi_exact_grad_num_l2
                  << " recon_full=" << m.recon_full_l2
                  << " iterations=" << m.iterations
                  << " converged=" << (m.converged ? "true" : "false") << '\n';

        if(!m.converged) {
            std::cerr << "FAIL Poisson decomposition solve did not converge at N=" << n << '\n';
            return 4;
        }
    }

    return 0;
}
