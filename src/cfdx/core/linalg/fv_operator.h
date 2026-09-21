#pragma once
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/mesh/mesh.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
namespace cfdx::core {
class FvConvectionDiffusionOperator final:public LinearOperatorBase{
public:
 FvConvectionDiffusionOperator(const Mesh& mesh,const std::vector<Vec3>& centres,const std::vector<Vec3>& Sf,const std::vector<double>& mass_flux,double diffusivity,double convection_factor=1.0)
 :mesh_(mesh),centres_(centres),Sf_(Sf),phi_(mass_flux),gamma_(diffusivity),convection_factor_(convection_factor){
  if(centres_.size()!=mesh_.n_cells()||Sf_.size()!=mesh_.n_faces()||phi_.size()!=mesh_.n_faces()||gamma_<0.0) throw std::invalid_argument("invalid fused FVM operator geometry");
 }
 std::size_t rows()const noexcept override{return mesh_.n_cells();}
 std::size_t cols()const noexcept override{return mesh_.n_cells();}
 bool has_diagonal()const noexcept override{return true;}
 void diagonal(Vector& d)const override{
  const std::size_t n=mesh_.n_cells(); d.resize(n); d.fill(0.0); const auto& own=mesh_.ownership();
  for(std::size_t f=0;f<mesh_.n_faces();++f){const std::size_t o=own.owner(f);const auto nr=own.neighbour(f);if(nr<0)continue;const std::size_t nb=static_cast<std::size_t>(nr);const double dist=(centres_[nb]-centres_[o]).mag();if(!(dist>0.0))throw std::runtime_error("degenerate FVM face");const double diff=gamma_*Sf_[f].mag()/dist;const double conv=convection_factor_*std::abs(phi_[f]);d(o)+=diff+0.5*conv;d(nb)+=diff+0.5*conv;}
  for(std::size_t i=0;i<n;++i)if(std::abs(d(i))<1e-30)d(i)=1.0;
 }
 void apply(const Vector& x,Vector& y)const override{
  const std::size_t n=mesh_.n_cells();if(x.size()!=n)throw std::invalid_argument("FVM operator input mismatch");if(y.size()!=n)y.resize(n);y.fill(0.0);const auto& own=mesh_.ownership();
  for(std::size_t f=0;f<mesh_.n_faces();++f){const std::size_t o=own.owner(f);const auto nr=own.neighbour(f);if(nr<0)continue;const std::size_t nb=static_cast<std::size_t>(nr);const double dist=(centres_[nb]-centres_[o]).mag();if(!(dist>0.0))throw std::runtime_error("degenerate FVM face");const double diff=gamma_*Sf_[f].mag()/dist;const double flux=convection_factor_*phi_[f];const double q=(diff+0.5*flux)*x(o)+(-diff+0.5*flux)*x(nb);y(o)+=q;y(nb)-=q;}
 }
private:const Mesh& mesh_;const std::vector<Vec3>& centres_;const std::vector<Vec3>& Sf_;const std::vector<double>& phi_;double gamma_,convection_factor_;
};
} // namespace cfdx::core
