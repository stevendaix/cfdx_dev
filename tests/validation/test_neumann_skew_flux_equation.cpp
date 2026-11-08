#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "common/test_harness.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <limits>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {
Mesh grid(std::size_t nx, std::size_t ny, std::size_t nz, double shear = 0.0,
          bool two_dimensional = false, bool normalized = false) {
    Mesh mesh;
    const auto v = [=](std::size_t i, std::size_t j, std::size_t k) {
        return (k * (ny + 1) + j) * (nx + 1) + i;
    };
    mesh.points().resize((nx+1)*(ny+1)*(nz+1));
    for (std::size_t k=0;k<=nz;++k)
        for (std::size_t j=0;j<=ny;++j)
            for (std::size_t i=0;i<=nx;++i)
                mesh.points().set(v(i,j,k),
                    static_cast<double>(i) / (normalized ? static_cast<double>(nx) : 1.0) +
                        shear * static_cast<double>(j) / (normalized ? static_cast<double>(ny) : 1.0),
                    static_cast<double>(j) / (normalized ? static_cast<double>(ny) : 1.0),
                    static_cast<double>(k) / (normalized ? static_cast<double>(nz) : 1.0));
    std::map<std::array<std::size_t,4>,std::size_t> seen;
    std::vector<std::vector<std::size_t>> cf(nx*ny*nz);
    const auto cindex = [=](std::size_t i,std::size_t j,std::size_t k) {
        return (k*ny+j)*nx+i;
    };
    for (std::size_t k=0;k<nz;++k)
        for (std::size_t j=0;j<ny;++j)
            for (std::size_t i=0;i<nx;++i) {
                const std::size_t cell=cindex(i,j,k);
                const auto add=[&](std::array<std::size_t,4> vertices) {
                    auto key=vertices;
                    std::sort(key.begin(),key.end());
                    auto it=seen.find(key);
                    if (it==seen.end()) {
                        const std::size_t id=mesh.faces().n_faces();
                        mesh.faces().push_face({vertices.begin(),vertices.end()});
                        mesh.ownership().resize(mesh.n_faces());
                        mesh.ownership().set_owner(id,cell);
                        mesh.ownership().set_neighbour(id,FaceOwnership::BOUNDARY);
                        seen.emplace(key,id);
                        cf[cell].push_back(id);
                    } else {
                        mesh.ownership().set_neighbour(it->second,
                            static_cast<std::int64_t>(cell));
                        cf[cell].push_back(it->second);
                    }
                };
                const auto a=v(i,j,k),b=v(i+1,j,k),c=v(i+1,j+1,k),d=v(i,j+1,k);
                const auto A=v(i,j,k+1),B=v(i+1,j,k+1),C=v(i+1,j+1,k+1),D=v(i,j+1,k+1);
                add({a,d,c,b}); add({A,B,C,D});
                add({a,b,B,A}); add({d,D,C,c});
                add({a,A,D,d}); add({b,c,C,B});
            }
    for (const auto& faces:cf) mesh.cells().push_cell(faces);
    mesh.boundary().add_patch("wall", PatchType::WALL);
    if (two_dimensional) mesh.boundary().add_patch("frontBack",PatchType::EMPTY);
    for (std::size_t f=0; f<mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f)>=0) continue;
        const auto offset=mesh.faces().offsets_data()[f];
        const auto end=mesh.faces().offsets_data()[f+1];
        const auto z_first=mesh.faces().vertices_data()[offset]/((nx+1)*(ny+1));
        // An EMPTY face is an entire extrusion cap, not a side face whose
        // first two vertices happen to lie at the same z.  The old test
        // helper silently discarded wall faces and invalidated the gauge
        // invariance test on sheared grids.
        bool on_empty=two_dimensional;
        for (auto vi=offset; vi<end; ++vi)
            on_empty = on_empty &&
                mesh.faces().vertices_data()[vi]/((nx+1)*(ny+1))==z_first;
        mesh.boundary().patch(on_empty?1:0).face_ids.push_back(f);
    }
    return mesh;
}


} // namespace
#include <iostream>
int main() {
  try {
    const double gamma=0.8;
    const Vec3 exact_gradient{1.25,-0.8,0.3};
    for (double shear : {0.0, 0.45, 0.75}) {
      auto mesh=grid(2,2,2,shear);
      const auto geo=build_fv_geometry(mesh);
      Field<double,Location::CELL> scalar(mesh.n_cells(),"T","K",1);
      Field<double,Location::FACE> mass_flux(mesh.n_faces(),"mdot","kg/s",1);
      Field<double,Location::CELL> source(mesh.n_cells(),"Su","W/m3",1);
      for (std::size_t c=0;c<mesh.n_cells();++c)
        scalar(c)=4+exact_gradient.dot(geo.cell_centres[c]);
      mass_flux.fill(0.0);
      source.fill(0.0);
      ScalarBoundaryFaceConditions conditions(mesh.n_faces());
      std::size_t face_test=mesh.n_faces();
      for (std::size_t f=0;f<mesh.n_faces();++f) {
        if (mesh.ownership().neighbour(f)>=0) continue;
        const auto cf=geo.face_centres[f];
        conditions.specified[f]=1;
        conditions.conditions[f]={ScalarBoundaryFaceCondition::Type::FIXED_VALUE,
                                  4+exact_gradient.dot(cf)};
        if (face_test==mesh.n_faces() && cf.y>1.99) face_test=f;
      }
      if(face_test==mesh.n_faces())
        throw std::runtime_error("missing prescribed Neumann face");
      const auto sf=geo.face_area_vectors[face_test];
      const double gn=exact_gradient.dot(sf)/sf.mag();
      conditions.conditions[face_test]={ScalarBoundaryFaceCondition::Type::FIXED_GRADIENT,gn};
      const auto reconstructed=reconstruct_scalar_transport_flux(
          mesh,geo,mass_flux,scalar,gamma,{},true,ConvectionScheme::UPWIND,&conditions);
      const auto correction=compute_nonorthogonal_diffusion_face_fluxes(
          mesh,geo,scalar,gamma,{},nullptr,nullptr,&conditions);
      const double expected=-gamma*sf.mag()*gn;
      if (std::abs(reconstructed(face_test)-expected)>1e-11 ||
          std::abs(correction[face_test])>1e-12)
        throw std::runtime_error("Neumann normal flux was changed by non-orthogonal correction");
      // Independently check the integrated physical face divergence against
      // A(phi)-b, not just against the same split-flux helper.
      const auto eq=assemble_scalar_equation(
          mesh,geo,mass_flux,gamma,source,source,{},true,nullptr,
          nullptr,nullptr,nullptr,ConvectionScheme::UPWIND,&scalar,
          &conditions,true);
      std::vector<double> divergence(mesh.n_cells(),0.0);
      for (std::size_t f=0;f<mesh.n_faces();++f) {
        const auto o=mesh.ownership().owner(f);
        divergence[o]+=reconstructed(f);
        const auto n=mesh.ownership().neighbour(f);
        if (n>=0) divergence[static_cast<std::size_t>(n)]-=reconstructed(f);
      }
      double row_error=0.0;
      for (std::size_t row=0;row<mesh.n_cells();++row) {
        double residual=-eq.rhs(row);
        for (auto k=eq.matrix.row_offsets_data()[row];
                  k<eq.matrix.row_offsets_data()[row+1];++k)
          residual+=eq.matrix.values_data()[k] * scalar(eq.matrix.columns_data()[k]);
        row_error=std::max(row_error,std::abs(residual-divergence[row]));
      }
      if(row_error>2e-11)
        throw std::runtime_error("assembled Neumann equation disagrees with physical face divergence");
      std::cout<<"NEUMANN_MMS shear="<<shear<<" flux_expected="<<expected
               <<" flux_actual="<<reconstructed(face_test)<<" equation_flux_error="<<row_error<<"\n";
      // Nonconstant boundary gamma must be the face gamma in gamma*S*gn.
      std::vector<double> wall_gamma(mesh.n_faces(),gamma);
      wall_gamma[face_test]=1.7;
      const auto reconstructed_variable=reconstruct_scalar_transport_flux(
          mesh,geo,mass_flux,scalar,gamma,{},true,ConvectionScheme::UPWIND,
          &conditions,nullptr,nullptr,&wall_gamma);
      if (std::abs(reconstructed_variable(face_test) + 1.7*sf.mag()*gn)>1e-11)
        throw std::runtime_error("Neumann did not apply the prescribed face diffusion coefficient");
    }
    // Patch-level Neumann values are also a complete normal flux, including
    // on skew faces when their gradient is not that of the cell field.
    auto mesh=grid(2,2,2,0.67);
    const auto geo=build_fv_geometry(mesh);
    Field<double,Location::CELL> scalar(mesh.n_cells(),"T","K",1);
    Field<double,Location::FACE> mass_flux(mesh.n_faces(),"mdot","kg/s",1);
    for (std::size_t c=0;c<mesh.n_cells();++c)
      scalar(c)=4+exact_gradient.dot(geo.cell_centres[c]);
    mass_flux.fill(0.0);
    const ScalarBoundaryConditions patch_bc{
      {"wall", {ScalarBoundaryType::FIXED_GRADIENT,0.0,0.4}}};
    const auto flux=reconstruct_scalar_transport_flux(
        mesh,geo,mass_flux,scalar,gamma,patch_bc);
    const auto correction=compute_nonorthogonal_diffusion_face_fluxes(
        mesh,geo,scalar,gamma,patch_bc);
    for (std::size_t f=0;f<mesh.n_faces();++f) {
      if (mesh.ownership().neighbour(f)>=0) continue;
      const double expected=-gamma*geo.face_area_vectors[f].mag()*0.4;
      if(std::abs(flux(f)-expected)>1e-11 || std::abs(correction[f])>1e-12)
         throw std::runtime_error("patch Neumann normal flux is not enforced exactly");
    }
    // A one-cell skewed diffusion problem with only prescribed normal fluxes
    // needs no gradient field for deferred tangential diffusion.  The extra
    // diagonal pins the otherwise singular pure-Neumann operator.
    {
      auto one_cell=grid(1,1,1,0.57);
      const auto geom=build_fv_geometry(one_cell);
      Field<double,Location::FACE> phi(one_cell.n_faces(),"flux","kg/s",1);
      Field<double,Location::CELL> zero(one_cell.n_cells(),"source","W/m3",1);
      phi.fill(0.0);
      zero.fill(0.0);
      const std::vector<double> pin{1.0};
      const auto equation=assemble_scalar_equation(
          one_cell,geom,phi,gamma,zero,zero,patch_bc,true,nullptr,
          &pin,nullptr,nullptr,ConvectionScheme::UPWIND,nullptr,nullptr,true);
      double flux_rhs=0.0;
      for(std::size_t f=0;f<one_cell.n_faces();++f)
        flux_rhs+=gamma*geom.face_area_vectors[f].mag()*0.4;
      if(std::abs(equation.rhs(0)-flux_rhs)>1e-12)
        throw std::runtime_error("skew all-Neumann assembly did not preserve specified flux");
    }
    {
      auto mesh=grid(2,2,2,0.65);
      const auto geom=build_fv_geometry(mesh);
      Field<double,Location::CELL> val(mesh.n_cells(),"T","K",1);
      Field<double,Location::FACE> mass_flux(mesh.n_faces(),"mdot","kg/s",1);
      const Vec3 gradient{1.25,-0.8,0.3};
      for(std::size_t c=0;c<mesh.n_cells();++c) val(c)=4.0+gradient.dot(geom.cell_centres[c]);
      mass_flux.fill(0.0);
      ScalarBoundaryFaceConditions bc(mesh.n_faces());
      std::vector<double> wall_gamma(mesh.n_faces(),0.95);
      for(std::size_t f=0;f<mesh.n_faces();++f) {
        if(mesh.ownership().neighbour(f)>=0) continue;
        bc.specified[f]=1;
        bc.conditions[f]={ScalarBoundaryFaceCondition::Type::FIXED_VALUE,
                          4.0+gradient.dot(geom.face_centres[f])};
      }
      const auto flux=reconstruct_scalar_transport_flux(
         mesh,geom,mass_flux,val,0.0,{},true,ConvectionScheme::UPWIND,
         &bc,nullptr,nullptr,&wall_gamma);
      double worst=0.0;
      for(std::size_t f=0;f<mesh.n_faces();++f) {
        if(mesh.ownership().neighbour(f)>=0) continue;
        const double want=-0.95*geom.face_area_vectors[f].dot(gradient);
        const double err=std::abs(flux(f)-want);
        if(err>worst) worst=err;
      }
      if(worst>2e-11)
        throw std::runtime_error("boundary-face diffusion was omitted from skew correction");
      Field<double,Location::CELL> src(mesh.n_cells(),"source","W/m3",1);
      src.fill(0.0);
      const auto equation=assemble_scalar_equation(
          mesh,geom,mass_flux,0.0,src,src,{},true,nullptr,
          nullptr,nullptr,nullptr,ConvectionScheme::UPWIND,&val,
          &bc,true,&wall_gamma);
      std::vector<double> divergence(mesh.n_cells(),0.0);
      for(std::size_t f=0;f<mesh.n_faces();++f) {
        const auto o=mesh.ownership().owner(f);
        divergence[o]+=flux(f);
        const auto nbr=mesh.ownership().neighbour(f);
        if(nbr>=0) divergence[static_cast<std::size_t>(nbr)]-=flux(f);
      }
      double balance_error=0.0;
      for(std::size_t row=0;row<mesh.n_cells();++row) {
        double residual=-equation.rhs(row);
        for(auto k=equation.matrix.row_offsets_data()[row];
                 k<equation.matrix.row_offsets_data()[row+1];++k)
          residual+=equation.matrix.values_data()[k]*val(equation.matrix.columns_data()[k]);
        balance_error=std::max(balance_error,std::abs(residual-divergence[row]));
      }
      if(balance_error>2e-11)
        throw std::runtime_error("boundary-face-only diffusion matrix and physical flux disagree");
      std::cout<<"BOUNDARY_FACE_ONLY_DIRICHLET_MMS max_flux_error="<<worst
               <<" equation_flux_error="<<balance_error<<"\n";
    }
    std::cout<<"NEUMANN_SKEW_FLUX_EQUATION PASS\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"NEUMANN_SKEW_FLUX_EQUATION FAIL: "<<error.what()<<"\n";
    return 1;
  }
}
