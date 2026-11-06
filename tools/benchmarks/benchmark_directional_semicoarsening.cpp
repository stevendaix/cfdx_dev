// Opt-in deterministic benchmark of entire FGMRES solve, not only one V-cycle.
// Every timing is admitted only if A*x=b and manufactured solution error are
// independently certified. Includes hierarchy build, iteration and total.
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/hypre_amg.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
namespace {
SparseMatrix make_matrix(std::size_t side,double ax,double ay,double ad) {
    SparseMatrix A(side*side,side*side);
    for (std::size_t y=0;y<side;++y)
        for (std::size_t x=0;x<side;++x) {
            const std::size_t i=y*side+x;
            double diag=1.0;
            const auto link=[&](std::size_t j,double a) {
                A.push_back(i,j,-a); diag+=a;
            };
            if (x) link(i-1,ax);
            if (x+1<side) link(i+1,ax);
            if (y) link(i-side,ay);
            if (y+1<side) link(i+side,ay);
            if (x && y) link(i-side-1,ad);
            if (x+1<side && y+1<side) link(i+side+1,ad);
            if (x+1<side && y) link(i-side+1,ad);
            if (x && y+1<side) link(i+side-1,ad);
            A.push_back(i,i,diag);
        }
    A.finalize();
    return A;
}

double elapsed_ms(const std::chrono::steady_clock::time_point& begin,
                  const std::chrono::steady_clock::time_point& end) {
    return std::chrono::duration<double,std::milli>(end-begin).count();
}
double median(std::vector<double> a) {
    std::sort(a.begin(),a.end());
    return a[a.size()/2];
}

struct Scenario { const char* name; double ax, ay, ad; };
struct Solver { const char* name; int id; };
std::unique_ptr<Preconditioner> make_pc(int id) {
    if (id==0) return std::make_unique<NativeBoomerAMGPreconditioner>(false,false);
    if (id==1) return std::make_unique<NativeAnisotropicSmoothedAggregationAMGPreconditioner>();
    if (id==2) return std::make_unique<NativeDirectionalSemiCoarseningAMGPreconditioner>();
    throw std::invalid_argument("bad preconditioner id");
}
}
int main(int argc,char** argv) {
    try {
        const std::size_t side=argc>1?std::stoul(argv[1]):32;
        const int trials=argc>2?std::stoi(argv[2]):3;
        const std::string selected=argc>3?argv[3]:"all";
        if (selected!="all" && selected!="anisotropic" &&
            selected!="isotropic")
            throw std::invalid_argument("scenario selector must be all, anisotropic, or isotropic");
        if (side<8 || side>128 || trials<1 || trials>20)
            throw std::invalid_argument("usage: benchmark_directional_semicoarsening side(8..128) trials(1..20)");
        const Scenario cases[]={{"strong_x_1e3",1000.,1.,0.01},
                                {"strong_y_1e3",1.,1000.,0.01},
                                {"oblique_1e3",0.5,0.5,1000.},
                                {"strong_x_1e5",1e5,1.,0.01},
                                {"isotropic",1.,1.,0.01}};
        const Solver solvers[]={{"native_boomer",0},{"anisotropic_SA_lines",1},
                                {"directional_semicoarsen_line_GS",2}};
        std::cout << "case,n,solver,setup_median_ms,fgmres_median_ms,total_median_ms,iterations,physical_residual,relative_solution_error,line_coverage,directional_pairs,operator_complexity\n";
        for (const auto& sc : cases) {
            if (selected=="anisotropic" && std::string(sc.name)=="isotropic") continue;
            if (selected=="isotropic" && std::string(sc.name)!="isotropic") continue;
            const auto A=make_matrix(side,sc.ax,sc.ay,sc.ad);
            Vector exact(A.n_rows(),0.0), b(A.n_rows(),0.0);
            for (std::size_t i=0;i<exact.size();++i)
                exact(i)=1.0+std::sin(0.018*static_cast<double>(i))+
                         0.2*std::cos(0.13*static_cast<double>(i));
            A.matvec_into(exact,b);
            for (const auto& sv : solvers) {
                std::vector<double> prep,solver,total;
                int iterations=-1;
                double residual=0.0, error=0.0,complexity=0.0;
                std::size_t coverage=0,pairs=0;
                for(int t=0;t<trials;++t) {
                    auto pc=make_pc(sv.id);
                    const auto t0=std::chrono::steady_clock::now();
                    if (!pc->setup(A)) {
                        const auto* native=dynamic_cast<NativeAMGPreconditioner*>(pc.get());
                        throw std::runtime_error(std::string("setup failure for ")+
                            sv.name+" "+sc.name+" " + (native?native->last_error():""));
                    }
                    const auto t1=std::chrono::steady_clock::now();
                    Vector x(A.n_rows(),0.0);
                    const auto result=solve_fgmres(A,b,x,50,300,1e-10,pc.get(),
                                                   PrecisionPolicy{},false);
                    const auto t2=std::chrono::steady_clock::now();
                    Vector Ax(A.n_rows(),0.0);
                    A.matvec_into(x,Ax);
                    double rn=0.0,bn=0.0,en=0.0,xn=0.0;
                    for (std::size_t i=0;i<b.size();++i) {
                        rn=std::hypot(rn,b(i)-Ax(i)); bn=std::hypot(bn,b(i));
                        en=std::hypot(en,x(i)-exact(i));
                        xn=std::hypot(xn,exact(i));
                    }
                    residual=rn/std::max(bn,1e-300);
                    error=en/std::max(xn,1e-300);
                    if (result.status!=SolverStatus::CONVERGED ||
                        !(residual<=1e-10) || !(error<2e-6)) {
                        std::cerr << "FAIL case=" << sc.name << " solver=" << sv.name
                            << " trial=" << t << " status=" << static_cast<int>(result.status)
                            << " iterations=" << result.iterations << " residual=" << residual
                            << " solution_error=" << error << '\n';
                        return 1; // fail-closed: no speedup from an invalid solve.
                    }
                    if (auto* a=dynamic_cast<NativeAMGPreconditioner*>(pc.get())) {
                        coverage=a->finest_line_coverage();
                        pairs=a->finest_directional_pairs();
                        complexity=a->operator_complexity();
                    }
                    iterations=static_cast<int>(result.iterations);
                    prep.push_back(elapsed_ms(t0,t1));
                    solver.push_back(elapsed_ms(t1,t2));
                    total.push_back(elapsed_ms(t0,t2));
                }
                std::cout << sc.name << ',' << A.n_rows() << ',' << sv.name << ','
                    << std::setprecision(8) << median(prep) << ',' << median(solver)
                    << ',' << median(total) << ',' << iterations << ',' << residual
                    << ',' << error << ',' << coverage << ',' << pairs << ','
                    << complexity << '\n';
            }
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "DIRECTIONAL_BENCH_FAIL: " << e.what() << '\n';
        return 2;
    }
}
