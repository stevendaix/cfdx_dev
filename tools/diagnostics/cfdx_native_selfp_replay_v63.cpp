@@ -0,0 +1,90 @@
// Reproduce the seventh Schur option on captured REAL coupled 4N MatrixMarket.
// Read-only original A,b,x0; compute the TRUE residual on the SAME A.
#include "cfdx/core/linalg/native_selfp_schur_v63.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
using namespace cfdx::core;
static SparseMatrix read_matrix(const char* path) {
    std::ifstream in(path);
    std::string line;
    if(!std::getline(in,line) ||
       line.find("%%MatrixMarket matrix coordinate real general")!=0)
        throw std::runtime_error("requires real general coordinate MatrixMarket");
    do {if(!std::getline(in,line))throw std::runtime_error("missing shape");}
    while(line.empty() || line[0]=='%');
    std::istringstream meta(line);
    std::size_t nr=0,nc=0,nnz=0;
    if(!(meta>>nr>>nc>>nnz) || !nr || nr!=nc)
        throw std::runtime_error("invalid saddle matrix dimensions");
    SparseMatrix A(nr,nc);
    for(std::size_t k=0;k<nnz;++k) {
        std::size_t i=0,j=0;double value=0;
        if(!(in>>i>>j>>value) || !i || !j || i>nr || j>nc)
            throw std::runtime_error("invalid matrix value");
        A.push_back(i-1,j-1,value);
    }
    A.finalize();
    return A;
}
static Vector read_vector(const char* path,std::size_t dimension) {
    std::ifstream in(path);
    std::string line;
    if(!std::getline(in,line) ||
       line.find("%%MatrixMarket matrix array real general")!=0)
        throw std::runtime_error("requires real general array MatrixMarket");
    do {if(!std::getline(in,line))throw std::runtime_error("missing vector dimensions");}
    while(line.empty() || line[0]=='%');
    std::istringstream meta(line);
    std::size_t nr=0,nc=0;
    if(!(meta>>nr>>nc) || nr!=dimension || nc!=1)
        throw std::runtime_error("vector dimensions mismatch");
    Vector result(dimension);
    for(std::size_t i=0;i<dimension;++i)
        if(!(in>>result(i)) || !std::isfinite(result(i)))
            throw std::runtime_error("invalid vector entry");
    return result;
}
int main(int argc,char** argv) {
    if(argc!=4) {
        std::cerr<<"usage: cfdx_native_selfp_replay_v63 A.mtx b.mtx x0.mtx\n";
        return 2;
    }
    try {
        const auto A=read_matrix(argv[1]);
        const auto b=read_vector(argv[2],A.n_rows());
        auto x=read_vector(argv[3],A.n_rows());
        if(A.n_rows()%4)return 2;
        NativeSelfpILUSchurV63 pc(A.n_rows()/4);
        if(!pc.setup(A)) {
            std::cerr<<"V63_SELFP_REPLAY setup_failed="<<pc.last_error()<<"\n";
            return 3;
        }
        KrylovControls controls;
        controls.restart_min=30;controls.restart_max=200;
        controls.adaptive_restart=false;
        const auto result=solve_fgmres(A,b,x,200,500,1e-9,&pc,controls,{},false);
        Vector ax(A.n_rows());A.matvec_into(x,ax);
        long double r2=0.,b2=0.;
        for(std::size_t i=0;i<b.size();++i) {
            const long double r=static_cast<long double>(b(i))-ax(i);
            r2+=r*r;b2+=static_cast<long double>(b(i))*b(i);
        }
        const auto relative=std::sqrt(static_cast<double>(r2/b2));
        std::cout<<std::scientific<<std::setprecision(14)
                 <<"V63_SELFP_REPLAY n="<<A.n_rows()
                 <<" nnz="<<A.nnz()
                 <<" schur_nnz="<<pc.schur_matrix().nnz()
                 <<" iterations="<<result.iterations
                 <<" true_relative="<<relative
                 <<" claimed_relative="<<result.residual_relative
                 <<" status="<<static_cast<int>(result.status)
                 <<" physical_qualified=0\n";
        if(!std::isfinite(relative))return 4;
        return relative<=1e-9 ? 0:5;
    }catch(const std::exception& error){std::cerr<<error.what()<<"\n";return 6;}
}
