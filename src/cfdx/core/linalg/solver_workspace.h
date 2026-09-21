#pragma once
#include "cfdx/core/linalg/vector.h"
#include <cstddef>
#include <vector>
namespace cfdx::core {
struct KrylovWorkspace{Vector r,r_hat,p,v,z,s,z_s,t,ap;void resize(std::size_t n){r.resize(n,0.0);r_hat.resize(n,0.0);p.resize(n,0.0);v.resize(n,0.0);z.resize(n,0.0);s.resize(n,0.0);z_s.resize(n,0.0);t.resize(n,0.0);ap.resize(n,0.0);}};
struct GmresWorkspace{std::size_t n=0,restart=0;Vector r,w,z;std::vector<double> basis,zbasis,h,cs,sn,g,y;
void resize(std::size_t n_,std::size_t m_){n=n_;restart=m_;r.resize(n,0.0);w.resize(n,0.0);z.resize(n,0.0);basis.assign((m_+1)*n_,0.0);zbasis.assign(m_*n_,0.0);h.assign((m_+1)*m_,0.0);cs.assign(m_,0.0);sn.assign(m_,0.0);g.assign(m_+1,0.0);y.assign(m_,0.0);}
double*v(std::size_t j){return basis.data()+j*n;}const double*v(std::size_t j)const{return basis.data()+j*n;}
double*zv(std::size_t j){return zbasis.data()+j*n;}const double*zv(std::size_t j)const{return zbasis.data()+j*n;}
};
}