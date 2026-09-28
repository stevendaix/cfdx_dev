#include "cfdx/physics/wall_distance.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
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
    double distance_formula_l2{};
    double distance_impl_l2{};
    double distance_formula_linf{};
    double distance_impl_linf{};
    double residual{};
    std::size_t iterations{};
    bool converged{};
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

std::vector<double> reconstruct_planar_distance(const std::vector<double>& phi,
                                                const WallDistanceGrid& g,
                                                double wall_y)
{
    std::vector<double> d(phi.size(),std::numeric_limits<double>::infinity());
    for(std::size_t id=0;id<phi.size();++id) {
        if(g.solid[id]) continue;
        const double y=g.points[id].y;
        // The manufactured cases are one-dimensional in y.  The exact
        // gradient is therefore obtained from the analytical derivative of
        // the manufactured phi, not from a second, potentially biased,
        // discrete gradient reconstruction.
        const double grad = std::abs(wall_y-y);
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
        const double lo=-2.0/(static_cast<double>(n)-1.0);
        const double h=(L-lo)/(static_cast<double>(n)-1.0);
        WallSurface surface;
        add_rect_plane(surface,0.0,-2.0,2.0,-2.0,2.0);
        const auto grid=make_wall_distance_grid(
            n,n,n,{lo,lo,lo},{h,h,h},
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
        const auto d_formula=reconstruct_planar_distance(phi,grid,0.0);
        std::size_t impl_it=0; double impl_residual=0.0;
        const auto d_impl=poisson_distance(bvh,grid,500,1.0,&impl_it,&impl_residual);
        rows.push_back({n,"single_wall",rel_l2(phi,phi_ref,grid.solid),
                        rel_linf(phi,phi_ref,grid.solid),
                        rel_l2(d_formula,d_ref,grid.solid),
                        rel_l2(d_impl,d_ref,grid.solid),
                        rel_linf(d_formula,d_ref,grid.solid),
                        rel_linf(d_impl,d_ref,grid.solid),
                        residual,it,residual<1e-8});
    }

    // 2) Parallel channel with two Dirichlet walls.
    //    Exact solution of -phi''=1, phi(0)=phi(1)=0:
    //       phi=y*(1-y)/2.
    //    Reconstruction gives min(y,1-y) exactly in the continuous problem.
    for(const auto n : sizes) {
        const double lo=-2.0/(static_cast<double>(n)-1.0);
        const double hi=1.0-lo;
        const double h=(hi-lo)/(static_cast<double>(n)-1.0);
        WallSurface surface;
        add_rect_plane(surface,0.0,-2.0,2.0,-2.0,2.0);
        add_rect_plane(surface,1.0,-2.0,2.0,-2.0,2.0);
        const auto grid=make_wall_distance_grid(
            n,n,n,{lo,lo,lo},{h,h,h},
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
        const auto d_formula=reconstruct_planar_distance(phi,grid,0.0);
        // The helper above assumes one wall for the gradient.  For the channel,
        // evaluate the exact signed gradient magnitude directly.
        for(std::size_t i=0;i<d_formula.size();++i) if(!grid.solid[i]) {
            const double y=grid.points[i].y;
            d_formula[i]=reconstruct_distance(phi[i],std::abs(0.5-y));
        }
        std::size_t impl_it=0; double impl_residual=0.0;
        const auto d_impl=poisson_distance(bvh,grid,500,1.0,&impl_it,&impl_residual);
        rows.push_back({n,"parallel_channel",rel_l2(phi,phi_ref,grid.solid),
                        rel_linf(phi,phi_ref,grid.solid),
                        rel_l2(d_formula,d_ref,grid.solid),
                        rel_l2(d_impl,d_ref,grid.solid),
                        rel_linf(d_formula,d_ref,grid.solid),
                        rel_linf(d_impl,d_ref,grid.solid),
                        residual,it,residual<1e-8});
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

    std::ofstream csv(output);
    if(!csv) return 3;
    csv << "N,case,phi_l2_relative,phi_linf_relative,distance_l2_relative,"
           "distance_linf_relative,residual_inf,iterations,converged\n";
    for(const auto& r:rows) {
        csv << r.n << ',' << r.case_name << ','
            << std::setprecision(16)
            << r.phi_l2 << ',' << r.phi_linf << ','
            << r.distance_formula_l2 << ',' << r.distance_impl_l2 << ','
            << r.distance_formula_linf << ',' << r.distance_impl_linf << ','
            << r.residual << ',' << r.iterations << ','
            << (r.converged ? "true" : "false") << '\n';
    }

    // The manufactured PDE qualification is intentionally strict: this test
    // validates the Poisson operator and the mathematical reconstruction
    // separately from the complex wall geometry benchmark.
    for(const auto& r:rows) {
        if(!r.converged) return 4;
        if(!(r.phi_l2<0.05 && r.distance_formula_l2<0.05 && r.distance_impl_l2<0.05)) {
            std::cerr << "FAIL manufactured Poisson qualification at N=" << r.n
                      << " case=" << r.case_name
                      << " phi_l2=" << r.phi_l2
                      << " formula_distance_l2=" << r.distance_formula_l2
                      << " implementation_distance_l2=" << r.distance_impl_l2 << "\n";
            return 5;
        }
    }
    return 0;
}
