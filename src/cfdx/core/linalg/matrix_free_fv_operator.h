#pragma once
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/mesh/mesh.h"
#include <cmath>
#include <stdexcept>
#include <vector>
namespace cfdx::core {
class FvDiffusionOperator final: public LinearOperatorBase {
public:
    FvDiffusionOperator(const Mesh& mesh,const std::vector<Vec3>& centres,const std::vector<Vec3>& Sf,double gamma):mesh_(mesh),cc_(centres),Sf_(Sf),gamma_(gamma){
        if(cc_.size()!=mesh_.n_cells()||Sf_.size()!=mesh_.n_faces()||gamma_<0.0)throw std::invalid_argument("invalid FvDiffusionOperator geometry");
    }
    std::size_t rows()const noexcept override{return mesh_.n_cells();}
    std::size_t cols()const noexcept override{return mesh_.n_cells();}
    void apply(const Vector& x,Vector& y)const override{
        const std::size_t n=mesh_.n_cells(); if(x.size()!=n)throw std::invalid_argument("FvDiffusionOperator size mismatch"); if(y.size()!=n)y.resize(n,0.0);
        for(std::size_t i=0;i<n;++i)y(i)=0.0; const auto& own=mesh_.ownership();
        for(std::size_t f=0;f<mesh_.n_faces();++f){std::size_t o=own.owner(f);auto nr=own.neighbour(f);double area=Sf_[f].mag();if(!(area>0.0))continue;
            if(nr>=0){std::size_t nb=static_cast<std::size_t>(nr);double d=(cc_[nb]-cc_[o]).mag();if(!(d>0.0))throw std::runtime_error("degenerate diffusion face");double a=gamma_*area/d;double q=a*(x(o)-x(nb));y(o)+=q;y(nb)-=q;}
        }
    }
private: const Mesh& mesh_; const std::vector<Vec3>& cc_; const std::vector<Vec3>& Sf_; double gamma_;
};
}