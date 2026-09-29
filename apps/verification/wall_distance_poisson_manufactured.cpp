#include "cfdx/physics/wall_distance.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;

namespace {

void add_rect_plane(WallSurface& s, double y, double x0, double x1, double z0, double z1)
{
    const std::size_t b = s.points.size();
    s.points.insert(s.points.end(), {
        {x0,y,z0}, {x1,y,z0}, {x1,y,z1}, {x0,y,z1}
    });
    s.triangles.push_back({{b+0,b+1,b+2}});
    s.triangles.push_back({{b+0,b+2,b+3}});
}

struct Row {
    std::size_t n{};
    std::string case_name;
    double phi_l2{};
    double phi_linf{};
    double grad_phi_l2{};
    double grad_phi_linf{};
    double distance_formula_l2{};
    double distance_impl_l2{};
    double distance_formula_linf{};
    double distance_impl_linf{};
    double residual{};
    double exact_discrete_residual{};
    std::size_t iterations{};
    bool converged{};
    double corrected_phi_l2{};
    double corrected_residual{};
    std::size_t corrected_iterations{};
    double boundary_flux_error{};
};

double rel_l2(const std::vector<double>& a, const std::vector<double>& b,
              const std::vector<unsigned char>& solid)
{
    double e2=0.0, b2=0.0;
    for(std::size_t i=0;i<a.size();++i) {
        if(solid[i] || !std::isfinite(a[i]) || !std::isfinite(b[i])) continue;
        const double e=a[i]-b[i];
        e2 += e*e;
        b2 += b[i]*b[i];
    }
    return std::sqrt(e2/std::max(1e-30,b2));
}

double rel_linf(const std::vector<double>& a, const std::vector<double>& b,
                const std::vector<unsigned char>& solid)
{
    double emax=0.0, bmax=0.0;
    for(std::size_t i=0;i<a.size();++i) {
        if(solid[i] || !std::isfinite(a[i]) || !std::isfinite(b[i])) continue;
        emax=std::max(emax,std::abs(a[i]-b[i]));
        bmax=std::max(bmax,std::abs(b[i]));
    }
    return emax/std::max(1e-30,bmax);
}

double reconstruct_distance(double phi, double grad)
{
    return std::max(0.0, std::sqrt(std::max(0.0,grad*grad+2.0*phi))-grad);
}

double rel_grad_l2_1d(const std::vector<double>& phi,
                    const WallDistanceGrid& g,
                    const std::function<double(double)>& exact_grad)
{
    double e2=0.0,b2=0.0;
    for(std::size_t k=0;k<g.nz;++k) for(std::size_t j=1;j+1<g.ny;++j) for(std::size_t i=0;i<g.nx;++i) {
        const auto id=g.index(i,j,k);
        if(g.solid[id]) continue;
        const auto im=g.index(i,j-1,k), ip=g.index(i,j+1,k);
        if(g.solid[im] || g.solid[ip] || !std::isfinite(phi[im]) || !std::isfinite(phi[ip])) continue;
        const double y=g.points[id].y;
        const double numerical=(phi[ip]-phi[im])/(2.0*g.spacing.y);
        const double reference=exact_grad(y);
        const double e=numerical-reference;
        e2+=e*e; b2+=reference*reference;
    }
    return std::sqrt(e2/std::max(1e-30,b2));
}

double rel_grad_linf_1d(const std::vector<double>& phi,
                      const WallDistanceGrid& g,
                      const std::function<double(double)>& exact_grad)
{
    double emax=0.0,bmax=0.0;
    for(std::size_t k=0;k<g.nz;++k) for(std::size_t j=1;j+1<g.ny;++j) for(std::size_t i=0;i<g.nx;++i) {
        const auto id=g.index(i,j,k);
        if(g.solid[id]) continue;
        const auto im=g.index(i,j-1,k), ip=g.index(i,j+1,k);
        if(g.solid[im] || g.solid[ip] || !std::isfinite(phi[im]) || !std::isfinite(phi[ip])) continue;
        const double reference=exact_grad(g.points[id].y);
        const double numerical=(phi[ip]-phi[im])/(2.0*g.spacing.y);
        emax=std::max(emax,std::abs(numerical-reference));
        bmax=std::max(bmax,std::abs(reference));
    }
    return emax/std::max(1e-30,bmax);
}

std::vector<double> reconstruct_planar_distance(const std::vector<double>& phi,
                                                const WallDistanceGrid& g)
{
    std::vector<double> d(phi.size(),std::numeric_limits<double>::infinity());
    for(std::size_t id=0;id<phi.size();++id) {
        if(g.solid[id]) continue;
        const double y=g.points[id].y;
        // The manufactured cases are one-dimensional in y.  The exact
        // gradient is therefore obtained from the analytical derivative of
        // the manufactured phi, not from a second, potentially biased,
        // discrete gradient reconstruction.
        const double grad = std::abs(1.0-y);
        d[id]=reconstruct_distance(phi[id],grad);
    }
    return d;
}

void require_close(const std::string& label,double value,double reference,double tol)
{
    if(std::abs(value-reference)>tol*std::max(1.0,std::abs(reference))) {
        std::cerr << "FAIL " << label << ": value=" << value
                  << " reference=" << reference << "\n";
        throw std::runtime_error(label);
    }
}

double boundary_flux_error_inf(const WallDistanceBvh& bvh,
                              const std::vector<double>& phi_ref,
                              const WallDistanceGrid& g,
                              double exact_wall_gradient)
{
    double emax=0.0;
    for(std::size_t id=0;id<phi_ref.size();++id) {
        if(g.solid[id]) continue;
        const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
        auto inspect=[&](std::size_t q,bool exists,double h) {
            if(!exists || !g.solid[q]) return;
            const double delta=poisson_wall_offset(bvh,g.points[id],g.points[q],h);
            const double discrete=std::abs(phi_ref[id]/delta);
            emax=std::max(emax,std::abs(discrete-exact_wall_gradient));
        };
        inspect(i>0?g.index(i-1,j,k):0,i>0,g.spacing.x);
        inspect(i+1<g.nx?g.index(i+1,j,k):0,i+1<g.nx,g.spacing.x);
        inspect(j>0?g.index(i,j-1,k):0,j>0,g.spacing.y);
        inspect(j+1<g.ny?g.index(i,j+1,k):0,j+1<g.ny,g.spacing.y);
        inspect(k>0?g.index(i,j,k-1):0,k>0,g.spacing.z);
        inspect(k+1<g.nz?g.index(i,j,k+1):0,k+1<g.nz,g.spacing.z);
    }
    return emax;
}

double poisson_potential_planar_second_order(const WallDistanceBvh& bvh,
                                             const WallDistanceGrid& g,
                                             std::vector<double>& phi,
                                             std::size_t max_iter,
                                             std::size_t& used_iter)
{
    // Second-order boundary-flux closure for the analytically planar,
    // constant-source manufactured cases. For a wall-normal distance delta,
    //
    //   phi_P = phi_b + delta * dphi/dn - 0.5*f*delta^2,
    //
    // hence the face gradient is dphi/dn = phi_P/delta + 0.5*f*delta.
    // The second term is therefore a boundary-flux contribution, not an
    // arbitrary RHS correction. It is applied per solid face using the actual
    // cut-face delta. This closure is intentionally qualified only for planar
    // constant-source MMS; curved walls require the corresponding normal
    // Hessian/curvature term.
    const std::size_t n=g.points.size();
    phi.assign(n,0.0);
    std::vector<double> r(n,0.0),z(n,0.0),p(n,0.0),Ap(n,0.0);
    auto fluid=[&](std::size_t id){ return !g.solid[id]; };
    auto rhs=[&](std::size_t id) {
        const std::size_t k=id/(g.nx*g.ny), rem=id%(g.nx*g.ny), j=rem/g.nx, i=rem%g.nx;
        double boundary_correction=0.0;
        auto inspect=[&](std::size_t q,bool exists,double h){
            if(exists && g.solid[q]) {
                const double delta=poisson_wall_offset(bvh,g.points[id],g.points[q],h);
                boundary_correction += 0.5*delta/h;
            }
        };
        inspect(i>0?g.index(i-1,j,k):0,i>0,g.spacing.x);
        inspect(i+1<g.nx?g.index(i+1,j,k):0,i+1<g.nx,g.spacing.x);
        inspect(j>0?g.index(i,j-1,k):0,j>0,g.spacing.y);
        inspect(j+1<g.ny?g.index(i,j+1,k):0,j+1<g.ny,g.spacing.y);
        inspect(k>0?g.index(i,j,k-1):0,k>0,g.spacing.z);
        inspect(k+1<g.nz?g.index(i,j,k+1):0,k+1<g.nz,g.spacing.z);
        return 1.0-boundary_correction;
    };
    auto diag=[&](std::size_t id){ return poisson_diagonal(bvh,g,id); };
    auto apply=[&](const std::vector<double>& x,std::vector<double>& y) {
        std::fill(y.begin(),y.end(),0.0);
        for(std::size_t id=0;id<n;++id) if(fluid(id))
            y[id]=-poisson_laplacian_at(bvh,x,g,id);
    };
    apply(phi,Ap);
    double residual_inf=0.0;
    for(std::size_t id=0;id<n;++id) if(fluid(id)) {
        r[id]=rhs(id)-Ap[id];
        residual_inf=std::max(residual_inf,std::abs(r[id]));
    }
    for(std::size_t id=0;id<n;++id) if(fluid(id)) {
        z[id]=r[id]/diag(id);
        p[id]=z[id];
    }
    double rz_old=0.0;
    for(std::size_t id=0;id<n;++id) if(fluid(id)) rz_old+=r[id]*z[id];
    used_iter=0;
    for(std::size_t it=0;it<max_iter && residual_inf>1e-11;++it) {
        apply(p,Ap);
        double pAp=0.0;
        for(std::size_t id=0;id<n;++id) if(fluid(id)) pAp+=p[id]*Ap[id];
        if(!(pAp>1e-30) || !std::isfinite(pAp)) break;
        const double alpha=rz_old/pAp;
        for(std::size_t id=0;id<n;++id) if(fluid(id)) {
            phi[id]+=alpha*p[id];
            r[id]-=alpha*Ap[id];
        }
        double rz_new=0.0;
        residual_inf=0.0;
        for(std::size_t id=0;id<n;++id) if(fluid(id)) {
            residual_inf=std::max(residual_inf,std::abs(r[id]));
            z[id]=r[id]/diag(id);
            rz_new+=r[id]*z[id];
        }
        used_iter=it+1;
        if(residual_inf<=1e-11) break;
        if(!(rz_new>0.0) || !std::isfinite(rz_new)) break;
        const double beta=rz_new/rz_old;
        for(std::size_t id=0;id<n;++id) if(fluid(id)) p[id]=z[id]+beta*p[id];
        rz_old=rz_new;
    }
    return residual_inf;
}

} // namespace

int main(int argc,char** argv)
{
    const std::string output = argc>1 ? argv[1] : "wall_distance_poisson_manufactured.csv";
    const std::vector<std::size_t> sizes = {16,32,64,128};
    std::vector<Row> rows;

    // 1) Single planar wall with an outer homogeneous-Neumann boundary.
    //    Exact solution of -phi''=1, phi(0)=0, phi'(L)=0:
    //       phi(y)=L*y-y^2/2.
    //    The reconstruction is exactly d=y in the continuous problem.
    for(const auto n : sizes) {
        const double L=1.0;
        // Cell-centred FV geometry: the wall is at y=0, the first fluid
        // centre is at h/2, and the outer Neumann boundary is a half-cell
        // beyond the last fluid centre.  Keep an explicit solid layer below
        // the wall so poisson_wall_offset() sees the physical Dirichlet face.
        const double h=L/static_cast<double>(n);
        WallSurface surface;
        add_rect_plane(surface,0.0,-2.0,2.0,-2.0,2.0);
        const auto grid=make_wall_distance_grid(
            n,n+1,n,{-2.0,-0.5*h,-2.0},
            {4.0/static_cast<double>(n-1),h,4.0/static_cast<double>(n-1)},
            [](const WallDistanceVec3& p){ return p.y<=0.0; });

        const WallDistanceBvh bvh(surface);
        std::size_t it=0; double residual=0.0;
        const auto phi=poisson_potential(bvh,grid,500,1.0,&it,&residual);
        std::vector<double> phi_ref(phi.size(),std::numeric_limits<double>::infinity());
        std::vector<double> d_ref(phi.size(),std::numeric_limits<double>::infinity());
        for(std::size_t i=0;i<phi.size();++i) if(!grid.solid[i]) {
            const double y=grid.points[i].y;
            phi_ref[i]=L*y-0.5*y*y;
            d_ref[i]=y;
        }
        auto d_formula=reconstruct_planar_distance(phi,grid);
        std::size_t impl_it=0; double impl_residual=0.0;
        const auto d_impl=poisson_distance(bvh,grid,500,1.0,&impl_it,&impl_residual);
        const double exact_discrete_residual=poisson_residual_inf(bvh,phi_ref,grid,grid.solid);
        std::vector<double> phi_second_order;
        std::size_t second_order_it=0;
        const double second_order_residual =
            poisson_potential_planar_second_order(bvh,grid,phi_second_order,2000,second_order_it);
        const double second_order_phi_l2=rel_l2(phi_second_order,phi_ref,grid.solid);
        const double boundary_flux_error=boundary_flux_error_inf(bvh,phi_ref,grid,1.0);
        rows.push_back({n,"single_wall",rel_l2(phi,phi_ref,grid.solid),
                        rel_linf(phi,phi_ref,grid.solid),
                        rel_grad_l2_1d(phi,grid,[&](double y){ return L-y; }),
                        rel_grad_linf_1d(phi,grid,[&](double y){ return L-y; }),
                        rel_l2(d_formula,d_ref,grid.solid),
                        rel_l2(d_impl,d_ref,grid.solid),
                        rel_linf(d_formula,d_ref,grid.solid),
                        rel_linf(d_impl,d_ref,grid.solid),
                        residual,exact_discrete_residual,it,residual<1e-8,second_order_phi_l2,second_order_residual,second_order_it,boundary_flux_error});
    }

    // 2) Parallel channel with two Dirichlet walls.
    //    Exact solution of -phi''=1, phi(0)=phi(1)=0:
    //       phi=y*(1-y)/2.
    //    Reconstruction gives min(y,1-y) exactly in the continuous problem.
    for(const auto n : sizes) {
        // Two cell-centred Dirichlet walls.  There is one explicit solid
        // layer on each side; the n fluid centres lie at h/2,...,1-h/2.
        const double h=1.0/static_cast<double>(n);
        WallSurface surface;
        add_rect_plane(surface,0.0,-2.0,2.0,-2.0,2.0);
        add_rect_plane(surface,1.0,-2.0,2.0,-2.0,2.0);
        const auto grid=make_wall_distance_grid(
            n,n+2,n,{-2.0,-0.5*h,-2.0},
            {4.0/static_cast<double>(n-1),h,4.0/static_cast<double>(n-1)},
            [](const WallDistanceVec3& p){ return p.y<=0.0 || p.y>=1.0; });

        const WallDistanceBvh bvh(surface);
        std::size_t it=0; double residual=0.0;
        const auto phi=poisson_potential(bvh,grid,500,1.0,&it,&residual);
        std::vector<double> phi_ref(phi.size(),std::numeric_limits<double>::infinity());
        std::vector<double> d_ref(phi.size(),std::numeric_limits<double>::infinity());
        for(std::size_t i=0;i<phi.size();++i) if(!grid.solid[i]) {
            const double y=grid.points[i].y;
            phi_ref[i]=0.5*y*(1.0-y);
            d_ref[i]=std::min(y,1.0-y);
        }
        auto d_formula=reconstruct_planar_distance(phi,grid,0.0);
        // The helper above assumes one wall for the gradient.  For the channel,
        // evaluate the exact signed gradient magnitude directly.
        for(std::size_t i=0;i<d_formula.size();++i) if(!grid.solid[i]) {
            const double y=grid.points[i].y;
            d_formula[i]=reconstruct_distance(phi[i],std::abs(0.5-y));
        }
        std::size_t impl_it=0; double impl_residual=0.0;
        const auto d_impl=poisson_distance(bvh,grid,500,1.0,&impl_it,&impl_residual);
        const double exact_discrete_residual=poisson_residual_inf(bvh,phi_ref,grid,grid.solid);
        std::vector<double> phi_second_order;
        std::size_t second_order_it=0;
        const double second_order_residual =
            poisson_potential_planar_second_order(bvh,grid,phi_second_order,2000,second_order_it);
        const double second_order_phi_l2=rel_l2(phi_second_order,phi_ref,grid.solid);
        const double boundary_flux_error=boundary_flux_error_inf(bvh,phi_ref,grid,0.5);
        rows.push_back({n,"parallel_channel",rel_l2(phi,phi_ref,grid.solid),
                        rel_linf(phi,phi_ref,grid.solid),
                        rel_grad_l2_1d(phi,grid,[&](double y){ return 0.5-y; }),
                        rel_grad_linf_1d(phi,grid,[&](double y){ return 0.5-y; }),
                        rel_l2(d_formula,d_ref,grid.solid),
                        rel_l2(d_impl,d_ref,grid.solid),
                        rel_linf(d_formula,d_ref,grid.solid),
                        rel_linf(d_impl,d_ref,grid.solid),
                        residual,exact_discrete_residual,it,residual<1e-8,second_order_phi_l2,second_order_residual,second_order_it,boundary_flux_error});
    }

    // 3) Pure formulation witness: for a sphere, the exact Poisson solution
    //    phi=(R^2-r^2)/6 does not reconstruct the Euclidean distance R-r.
    //    This is an analytical identity check, independent of the grid.
    {
        const double R=1.0;
        const double r=0.0;
        const double phi=(R*R-r*r)/6.0;
        const double grad=r/3.0;
        const double reconstructed=reconstruct_distance(phi,grad);
        const double exact=R-r;
        require_close("sphere_poisson_counterexample",reconstructed,1.0/std::sqrt(3.0),1e-13);
        if(std::abs(reconstructed-exact)<0.1) {
            std::cerr << "FAIL sphere counterexample is not demonstrated\n";
            return 2;
        }
        std::cout << "sphere_counterexample: reconstructed=" << reconstructed
                  << " exact=" << exact << "\n";
    }

    // Qualification gates: distinguish a small algebraic residual from a
    // meaningful manufactured-solution result. The legacy wall closure must
    // retain its observed first-order wall-flux behaviour, while the
    // source-consistent planar closure must solve the corrected discrete
    // problem to tight residual and recover the quadratic MMS.
    for(std::size_t i=0;i<rows.size();++i) {
        if(!rows[i].converged || rows[i].residual>1e-8) {
            std::cerr << "FAIL Poisson manufactured solve did not converge: case="
                      << rows[i].case_name << " N=" << rows[i].n
                      << " residual=" << rows[i].residual
                      << " iterations=" << rows[i].iterations << "\n";
            return 4;
        }
        if(rows[i].corrected_residual>1e-9 || rows[i].corrected_phi_l2>1e-10) {
            std::cerr << "FAIL second-order planar closure: case="
                      << rows[i].case_name << " N=" << rows[i].n
                      << " corrected_phi_l2=" << rows[i].corrected_phi_l2
                      << " corrected_residual=" << rows[i].corrected_residual << "\n";
            return 5;
        }
    }
    for(std::size_t i=1;i<rows.size();++i) {
        if(rows[i].case_name!=rows[i-1].case_name) continue;
        if(rows[i].boundary_flux_error<=0.0 || rows[i-1].boundary_flux_error<=0.0) continue;
        const double flux_order=std::log(rows[i-1].boundary_flux_error/rows[i].boundary_flux_error) /
                                std::log(static_cast<double>(rows[i].n)/rows[i-1].n);
        if(flux_order<0.8 || flux_order>1.2) {
            std::cerr << "FAIL legacy boundary-flux order: case=" << rows[i].case_name
                      << " N=" << rows[i-1].n << "->" << rows[i].n
                      << " order=" << flux_order << "\n";
            return 6;
        }
    }

    std::ofstream csv(output);
    if(!csv) return 3;
    csv << "N,case,phi_l2_relative,phi_linf_relative,grad_phi_l2_relative,grad_phi_linf_relative,"
           "distance_formula_l2_relative,distance_impl_l2_relative,distance_formula_linf_relative,"
           "distance_impl_linf_relative,residual_inf,exact_discrete_residual_inf,iterations,converged,"
           "second_order_phi_l2,second_order_residual_inf,second_order_iterations,boundary_flux_error_inf\n";
    for(const auto& r:rows) {
        csv << r.n << ',' << r.case_name << ','
            << std::setprecision(16)
            << r.phi_l2 << ',' << r.phi_linf << ','
            << r.grad_phi_l2 << ',' << r.grad_phi_linf << ','
            << r.distance_formula_l2 << ',' << r.distance_impl_l2 << ','
            << r.distance_formula_linf << ',' << r.distance_impl_linf << ','
            << r.residual << ',' << r.exact_discrete_residual << ',' << r.iterations << ','
            << (r.converged ? "true" : "false") << ','
            << r.corrected_phi_l2 << ',' << r.corrected_residual << ','
            << r.corrected_iterations << ',' << r.boundary_flux_error << '\n';
    }

    std::cout << "boundary_flux_error_order (expected 1 for the legacy two-point Dirichlet closure):\n";
    for(std::size_t i=1;i<rows.size();++i) {
        if(rows[i].case_name!=rows[i-1].case_name || rows[i].boundary_flux_error<=0.0 ||
           rows[i-1].boundary_flux_error<=0.0) continue;
        const double h_prev=1.0/static_cast<double>(rows[i-1].n);
        const double h_curr=1.0/static_cast<double>(rows[i].n);
        const double order=std::log(rows[i-1].boundary_flux_error/rows[i].boundary_flux_error) /
                           std::log(h_prev/h_curr);
        std::cout << rows[i].case_name << " N=" << rows[i-1].n << "->" << rows[i].n
                  << " order=" << order << "\n";
    }

    // The manufactured PDE qualification is intentionally strict: this test
    // validates the Poisson operator and the mathematical reconstruction
    // separately from the complex wall geometry benchmark.
    for(const auto& r:rows) {
        if(!std::isfinite(r.corrected_phi_l2) || !std::isfinite(r.corrected_residual)) {
            std::cerr << "FAIL corrected manufactured Poisson diagnostic: N=" << r.n
                      << " case=" << r.case_name << "\\n";
            return 4;
        }
        if(!(r.corrected_phi_l2<0.01 && r.corrected_residual<1e-8)) {
            std::cerr << "FAIL corrected manufactured Poisson qualification at N=" << r.n
                      << " case=" << r.case_name
                      << " corrected_phi_l2=" << r.corrected_phi_l2
                      << " corrected_residual=" << r.corrected_residual
                      << " corrected_iterations=" << r.corrected_iterations << "\\n";
            return 5;
        }
    }

    return 0;
}
