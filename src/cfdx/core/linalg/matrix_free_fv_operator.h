#pragma once
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/mesh/mesh.h"
#include <cmath>
#include <stdexcept>
#include <vector>
namespace cfdx::core {
class FvDiffusionOperator final:public LinearOperatorBase{
public:
 FvDiffusionOperator(const Mesh& mesh,const std::vector<Vec3>& centres,const std::vector<Vec3>& Sf,double gamma):mesh_(mesh),cc_(centres),Sf_(Sf),gamma_(gamma){if(cc_.size()!=mesh_.n_cells()||Sf_.size()!=mesh_.n_faces()||gamma_<0.0)throw std::invalid_argument("invalid FvDiffusionOperator geometry");}
 std::size_t rows()const noexcept override{return mesh_.n_cells();}
 std::size_t cols()const noexcept override{return mesh_.n_cells();}
 bool has_diagonal()const noexcept override{return true;}
 void diagonal(Vector& d)const override{d.resize(mesh_.n_cells());d.fill(0.0);const auto& own=mesh_.ownership();for(std::size_t f=0;f<mesh_.n_faces();++f){const std::size_t o=own.owner(f);const auto nr=own.neighbour(f);if(nr<0)continue;const std::size_t nb=static_cast<std::size_t>(nr);const double dist=(cc_[nb]-cc_[o]).mag();if(!(dist>0.0))throw std::runtime_error("degenerate diffusion face");const double a=gamma_*Sf_[f].mag()/dist;d(o)+=a;d(nb)+=a;}for(std::size_t i=0;i<d.size();++i)if(std::abs(d(i))<1e-30)d(i)=1.0;}
 void apply(const Vector& x,Vector& y)const override{const std::size_t n=mesh_.n_cells();if(x.size()!=n)throw std::invalid_argument("FvDiffusionOperator size mismatch");if(y.size()!=n)y.resize(n);y.fill(0.0);const auto& own=mesh_.ownership();for(std::size_t f=0;f<mesh_.n_faces();++f){const std::size_t o=own.owner(f);const auto nr=own.neighbour(f);const double area=Sf_[f].mag();if(!(area>0.0))continue;if(nr>=0){const std::size_t nb=static_cast<std::size_t>(nr);const double d=(cc_[nb]-cc_[o]).mag();if(!(d>0.0))throw std::runtime_error("degenerate diffusion face");const double a=gamma_*area/d;const double q=a*(x(o)-x(nb));y(o)+=q;y(nb)-=q;}}}
private:const Mesh& mesh_;const std::vector<Vec3>& cc_;const std::vector<Vec3>& Sf_;double gamma_;
};
} // namespace cfdx::core
