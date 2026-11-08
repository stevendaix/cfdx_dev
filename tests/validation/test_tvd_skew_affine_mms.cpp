// Multidimensional TVD affine-skew MMS and matrix/face-flux identity.
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
#include <stdexcept>

int main() {
    auto m = grid(3,3,2,0.35);
    const auto v=[](std::size_t i,std::size_t j,std::size_t k) {
        return (k*4+j)*4+i;
    };
    for (std::size_t k=0; k<=2; ++k) {
        const std::size_t id=v(1,1,k);
        m.points().set(id, m.points().x_data()[id],
                       m.points().y_data()[id]+0.65,
                       m.points().z_data()[id]);
    }
    const auto g=build_fv_geometry(m);
    Field<double,Location::CELL> phi(m.n_cells(),"affine_y","1",1);
    Field<double,Location::CELL> zero(m.n_cells(),"source","1",1);
    Field<double,Location::FACE> flux(m.n_faces(),"flux","kg/s",1);
    zero.fill(0.0);
    for (std::size_t c=0;c<m.n_cells();++c) phi(c)=g.cell_centres[c].y;
    for (std::size_t f=0;f<m.n_faces();++f)
        flux(f)=m.ownership().neighbour(f)>=0 ? 1.0 : 0.0;
    ScalarBoundaryFaceConditions bcs(m.n_faces());
    for (std::size_t f=0;f<m.n_faces();++f)
        if (m.ownership().neighbour(f)<0) {
            bcs.specified[f]=1;
            bcs.conditions[f]={ScalarBoundaryFaceCondition::Type::FIXED_VALUE,
                               g.face_centres[f].y};
        }
    const auto grad=compute_boundary_aware_wls_gradient(phi,m,g,{},nullptr,&bcs);
    double max_gradient_error=0.0;
    for (std::size_t c=0;c<m.n_cells();++c) {
        max_gradient_error=std::max(max_gradient_error,
            std::abs(grad(c,0))+std::abs(grad(c,1)-1)+std::abs(grad(c,2)));
    }
    if (max_gradient_error>1e-11)
        throw std::runtime_error("manufactured gradient is not affine exact");

    for (const auto scheme : {ConvectionScheme::TVD,
                              ConvectionScheme::TVD_VANLEER,
                              ConvectionScheme::TVD_VANALBADA}) {
        const auto transported=reconstruct_scalar_transport_flux(
            m,g,flux,phi,0.0,{},true,scheme,&bcs);
        double face_error=0.0;
        for (std::size_t f=0;f<m.n_faces();++f)
            if (m.ownership().neighbour(f)>=0)
                face_error=std::max(face_error,
                    std::abs(transported(f)-g.face_centres[f].y));
        if (face_error>1e-11)
            throw std::runtime_error("TVD lost affine exactness on a conforming skew mesh");

        const std::vector<double> extra_diag(m.n_cells(),10.0);
        const auto eq=assemble_scalar_equation(
            m,g,flux,0.0,zero,zero,{},true,nullptr,
            &extra_diag,nullptr,nullptr,scheme,&phi,&bcs,false);
        std::vector<double> divergence(m.n_cells(),0.0);
        std::vector<double> div_f(m.n_cells(),0.0);
        for (std::size_t f=0;f<m.n_faces();++f) {
            const auto owner=m.ownership().owner(f);
            divergence[owner]+=transported(f);
            div_f[owner]+=flux(f);
            const auto nb=m.ownership().neighbour(f);
            if (nb>=0) {
                divergence[static_cast<std::size_t>(nb)]-=transported(f);
                div_f[static_cast<std::size_t>(nb)]-=flux(f);
            }
        }
        double equation_error=0.0;
        for (std::size_t row=0;row<m.n_cells();++row) {
            double assembled=-eq.rhs(row);
            for (auto k=eq.matrix.row_offsets_data()[row];
                 k<eq.matrix.row_offsets_data()[row+1];++k)
                assembled+=eq.matrix.values_data()[k]
                    *phi(eq.matrix.columns_data()[k]);
            const double expected=divergence[row]-phi(row)*div_f[row]
                +extra_diag[row]*phi(row);
            equation_error=std::max(equation_error,std::abs(assembled-expected));
        }
        if (equation_error>2e-11)
            throw std::runtime_error("TVD assembly and reconstructed face flux disagree");
        std::cout<<"TVD_SKEW_AFFINE scheme="<<static_cast<int>(scheme)
                 <<" face_Linf="<<face_error
                 <<" matrix_flux_Linf="<<equation_error<<"\n";
    }
    // No artificially enlarged extrema: even with a skew face, the result
    // must remain within the neighbouring one-ring cell/boundary envelope.
    Field<double,Location::CELL> step(m.n_cells(),"step","1",1);
    for (std::size_t c=0;c<m.n_cells();++c)
        step(c)=g.cell_centres[c].y>=1.5 ? 1.0 : 0.0;
    for (std::size_t f=0;f<m.n_faces();++f)
        if (m.ownership().neighbour(f)<0)
            bcs.conditions[f].value=g.face_centres[f].y>=1.5 ? 1.0 : 0.0;
    for (const auto scheme : {ConvectionScheme::TVD,
                              ConvectionScheme::TVD_VANLEER,
                              ConvectionScheme::TVD_VANALBADA}) {
        const auto limited=reconstruct_scalar_transport_flux(
            m,g,flux,step,0.0,{},true,scheme,&bcs);
        for (std::size_t f=0;f<m.n_faces();++f)
            if (m.ownership().neighbour(f)>=0 &&
                (!(limited(f)>=-1e-12) || !(limited(f)<=1+1e-12)))
                throw std::runtime_error("TVD introduced a new local extremum");
    }
    std::cout << "TVD_SKEW_AFFINE_MMS: PASS\n";
}
