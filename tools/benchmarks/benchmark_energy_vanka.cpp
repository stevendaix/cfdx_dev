// Opt-in comparative benchmark. Never promotes a solver by timing alone:
// every row certifies the original FP64 ||A*x-b|| / ||b|| and solution error.
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/vanka_patch_preconditioner.h"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <tuple>
using namespace cfdx::core;
namespace {
SparseMatrix elliptic(std::size_t side,double h,double v,double oblique) {
    SparseMatrix A(side*side,side*side);
    for(std::size_t y=0;y<side;++y) for(std::size_t x=0;x<side;++x) {
        const auto i=y*side+x;
        double diag=0.5;
        auto edge=[&](std::size_t j,double w){A.push_back(i,j,-w);diag+=w;};
        if(x) edge(i-1,h); if(x+1<side) edge(i+1,h);
        if(y) edge(i-side,v); if(y+1<side) edge(i+side,v);
        if(x&&y) edge(i-side-1,oblique);
        if(x+1<side&&y+1<side) edge(i+side+1,oblique);
        if(x+1<side&&y) edge(i-side+1,oblique);
        if(x&&y+1<side) edge(i+side-1,oblique);
        A.push_back(i,i,diag);
    }
    A.finalize(); return A;
}
SparseMatrix coupled(std::size_t nc,double anisotropy) {
    SparseMatrix A(4*nc,4*nc);
    for(std::size_t i=0;i<nc;++i) {
        const auto nv=3*nc;
        const double scale=i%2==0?anisotropy:1.0;
        A.push_back(i,i,5.0*scale);
        if(i) A.push_back(i,i-1,-0.25);
        if(i+1<nc) A.push_back(i,i+1,-0.25);
        A.push_back(i,nv+i,1.0);
        if(i+1<nc) A.push_back(i,nv+i+1,-1.0);
        A.push_back(nc+i,nc+i,4.0*scale);
        A.push_back(2*nc+i,2*nc+i,4.0*scale);
        if(i==0) A.push_back(nv,nv,1.0);
        else {
            A.push_back(nv+i,i,-1.0);
            A.push_back(nv+i,i-1,1.0);
            A.push_back(nv+i,nv+i,0.1);
        }
    }
    A.finalize(); return A;
}
Vector exact(std::size_t n){
    Vector x(n,0.0);
    for(std::size_t i=0;i<n;++i)
        x(i)=1.0+0.6*std::sin(0.011*static_cast<double>(i))+
             0.25*std::cos(0.23*static_cast<double>(i));
    return x;
}
double error(const Vector& x,const Vector& y){
    double num=0,den=0;
    for(std::size_t i=0;i<x.size();++i){
        num=std::hypot(num,x(i)-y(i));den=std::hypot(den,y(i));
    }
    return num/den;
}
void row(const std::string& kind,const std::string& matrix,const std::string& method,
         std::size_t n,std::size_t iterations,double residual,double relative_error,
         double elapsed_ms,bool converged){
    std::cout<<kind<<','<<matrix<<','<<method<<','<<n<<','<<iterations<<','
             <<std::setprecision(17)<<residual<<','<<relative_error<<','
             <<elapsed_ms<<','<<(converged?"PASS":"FAIL")<<'\n';
    if(!converged) throw std::runtime_error("uncertified benchmark: "+matrix+" "+method);
}
}
int main(){
    constexpr double tolerance=1e-10;
    std::cout<<"family,matrix,method,n,iterations,relative_true_residual,relative_solution_error,wall_ms,status\n";
    for(auto [name,h,v,o] : std::vector<std::tuple<std::string,double,double,double>>{
            {"x_1e3",1000,1,0.01},{"y_1e3",1,1000,0.01},
            {"oblique_1e3",1,1,1000}}){
        for(std::size_t side : {16u,32u}) {
            const auto A=elliptic(side,h,v,o);
            const auto truth=exact(A.n_rows());
            Vector b(A.n_rows(),0.0); A.matvec_into(truth,b);
            for(const auto preconditioner : {
                    PreconditionerModel::NativeAMG,
                    PreconditionerModel::DirectionalSemiCoarseningAMG,
                    PreconditionerModel::EnergyMinimizingAMG}){
                LinearSolverRequest request;
                request.preconditioner=preconditioner; request.gmres_restart=50;
                Vector x(A.n_rows(),0.0);
                const auto t0=std::chrono::steady_clock::now();
                const auto result=solve_linear_system(A,b,x,
                    LinearProblemKind::PressurePoisson,request,400,tolerance);
                const auto t1=std::chrono::steady_clock::now();
                const double ms=std::chrono::duration<double,std::milli>(t1-t0).count();
                const double relative=physical_relative_residual(A,x,b);
                const double relerror=error(x,truth);
                const bool ok=result.result.status==SolverStatus::CONVERGED &&
                    relative<=tolerance && relerror<=2e-7;
                row("elliptic",name,to_string(preconditioner),A.n_rows(),
                    result.result.iterations,relative,relerror,ms,ok);
            }
        }
    }
    for(double strength : {1.0,1000.0}) {
        const auto A=coupled(30,strength);
        const auto truth=exact(A.n_rows());Vector b(A.n_rows(),0.0); A.matvec_into(truth,b);
        LinearOperator op;
        op.size=A.n_rows();
        op.apply=[&](const Vector& in,Vector& out){
            if(out.size()!=A.n_rows()) out.resize(A.n_rows());
            A.matvec_into(in,out);
        };
        op.apply_fp64=op.apply;
        for(bool use_vanka : {false,true}) {
            VankaPressurePatchPreconditioner vanka(30);
            Vector x(A.n_rows(),0.0);
            const auto start=std::chrono::steady_clock::now();
            if(use_vanka&&!vanka.setup(A))
                throw std::runtime_error("Vanka setup: "+vanka.last_error());
            constexpr double coupled_tolerance=1e-12;
            const auto result=solve_fgmres(op,b,x,50,400,coupled_tolerance,
                use_vanka?static_cast<Preconditioner*>(&vanka):nullptr);
            const auto finish=std::chrono::steady_clock::now();
            const double residual=physical_relative_residual(A,x,b);
            const double relerror=error(x,truth);
            const bool ok=result.status==SolverStatus::CONVERGED &&
                residual<=coupled_tolerance && relerror<=1e-8;
            row("coupled",strength>1?"stretched":"normal",
                use_vanka?"vanka_pressure_patch":"none",A.n_rows(),
                result.iterations,residual,relerror,
                std::chrono::duration<double,std::milli>(finish-start).count(),ok);
        }
    }
    return 0;
}