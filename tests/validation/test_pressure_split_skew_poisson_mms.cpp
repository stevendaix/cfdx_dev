// Independent solution-level benchmark for experimental pressure face splits.
// Solve -Delta(phi)=-4 on a sheared square with exact Dirichlet data.
// The gradient in the deferred term is supplied analytically (MMS): this
// qualifies the face split and CG solution, NOT the production pressure
// reconstruction or Rhie-Chow coupling.
#include "cfdx/physics/pressure_flux_split.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {
void require(bool cond, const char* message) {
    if (!cond) throw std::runtime_error(message);
}
Vec3 center(std::size_t i, std::size_t j, int n, double shear) {
    const double h=1.0/static_cast<double>(n);
    const double x=(static_cast<double>(i)+0.5)*h;
    const double y=(static_cast<double>(j)+0.5)*h;
    return {x+shear*y,y,0.0};
}
double exact(const Vec3& x, bool affine) {
    return affine ? 0.3+2.0*x.x-3.0*x.y : x.x*x.x+x.y*x.y;
}
Vec3 gradient(const Vec3& x, bool affine) {
    return affine ? Vec3{2.0,-3.0,0.0} : Vec3{2.0*x.x,2.0*x.y,0.0};
}
struct Metrics {double error=0.0; double true_residual=0.0; std::size_t iterations=0;};
Metrics solve_grid(int n, double shear, PressureFluxSplitScheme scheme, bool affine) {
    const std::size_t count=static_cast<std::size_t>(n*n);
    const double h=1.0/static_cast<double>(n), volume=h*h;
    SparseMatrix A(count,count);
    Vector b(count,0.0), p(count,0.0);
    auto id=[n](int i,int j){return static_cast<std::size_t>(j*n+i);};
    const Vec3 response{1.0,1.0,1.0};
    for(int j=0;j<n;++j) for(int i=0;i<n;++i) {
        const std::size_t row=id(i,j);
        const Vec3 xc=center(static_cast<std::size_t>(i),static_cast<std::size_t>(j),n,shear);
        b(row)=affine ? 0.0 : -4.0*volume;
        // Cell-face-area vectors oriented OUTWARDS and boundary midpoints.
        const std::array<std::array<int,2>,4> dirs{{{{-1,0}},{{1,0}},{{0,-1}},{{0,1}}}};
        for(const auto& dir:dirs) {
            const int ni=i+dir[0], nj=j+dir[1];
            const bool internal=ni>=0 && ni<n && nj>=0 && nj<n;
            const Vec3 S = dir[0]!=0
                ? Vec3{static_cast<double>(dir[0])*h,-static_cast<double>(dir[0])*shear*h,0.0}
                : Vec3{0.0,static_cast<double>(dir[1])*h,0.0};
            const Vec3 face{xc.x+0.5*h*(dir[0]+shear*dir[1]),
                            xc.y+0.5*h*dir[1],0.0};
            const Vec3 other = internal
                ? center(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj),n,shear)
                : face;
            const Vec3 d = other-xc;
            const auto split=split_pressure_flux(S,d,response,gradient(face,affine),scheme);
            A.push_back(row,row,split.implicit_coefficient);
            b(row)+=split.deferred_flux;
            if(internal) A.push_back(row,id(ni,nj),-split.implicit_coefficient);
            else b(row)+=split.implicit_coefficient*exact(face,affine);
        }
    }
    A.finalize();
    const auto r=solve_cg(A,b,p,4000,1e-11);
    require(r.status==SolverStatus::CONVERGED,"CG did not converge on skew Poisson");
    const auto Ap=A.matvec(p);
    double r2=0.0,b2=0.0,error2=0.0;
    for (int j=0;j<n;++j) for(int i=0;i<n;++i) {
        const auto row=id(i,j);
        const double residual=Ap[row]-b(row);
        r2+=residual*residual; b2+=b(row)*b(row);
        const double diff=p(row)-exact(center(static_cast<std::size_t>(i), static_cast<std::size_t>(j), n, shear),affine);
        error2+=diff*diff;
    }
    const double relres=std::sqrt(r2)/std::max(std::sqrt(b2),1e-300);
    require(relres<1e-9,"independent true residual above gate");
    return {std::sqrt(error2/static_cast<double>(count)),relres,r.iterations};
}
}
int main() {
    try {
        std::cout<<std::setprecision(12);
        for (const double shear : {0.0,1.0,3.732050807568877}) {
            for (const auto scheme : {PressureFluxSplitScheme::MINIMUM_CORRECTION,
                                      PressureFluxSplitScheme::OVER_RELAXED}) {
                const auto affine=solve_grid(8,shear,scheme,true);
                require(affine.error<1e-9,"affine manufactured solution not exact");
                const auto coarse=solve_grid(8,shear,scheme,false);
                const auto medium=solve_grid(16,shear,scheme,false);
                const auto fine=solve_grid(32,shear,scheme,false);
                std::cout<<"POISSON_SKEW_MMS shear="<<shear
                         <<" scheme="<<(scheme==PressureFluxSplitScheme::OVER_RELAXED?"over":"minimum")
                         <<" e8="<<coarse.error<<" e16="<<medium.error
                         <<" e32="<<fine.error<<" affine="<<affine.error
                         <<" relres32="<<fine.true_residual
                         <<" cg32="<<fine.iterations<<'\n';
                require(fine.error<medium.error && medium.error<coarse.error,
                        "quadratic solution error does not decrease on refinement");
            }
        }
        std::cout<<"PRESSURE_SPLIT_SKEW_POISSON_MMS: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr<<"PRESSURE_SPLIT_SKEW_POISSON_MMS: FAIL "<<e.what()<<'\n';
        return 1;
    }
}
