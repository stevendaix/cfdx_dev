// Multidimensional TVD affine-skew MMS and matrix/face-flux identity.
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/momentum_same_state_audit.h"
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
    // Relaxed true linear residual and unrelaxed frozen equation must be
    // distinguished, even when a solver would report CONVERGED on A_rel.
    {
        ScalarEquation equation;
        equation.matrix=SparseMatrix(2,2);
        equation.matrix.push_back(0,0,4.0);
        equation.matrix.push_back(0,1,-1.0);
        equation.matrix.push_back(1,0,-1.0);
        equation.matrix.push_back(1,1,6.0);
        equation.matrix.finalize();
        equation.diagonal={4.0,6.0};
        equation.rhs=Vector(2);
        equation.rhs(0)=4.0;
        equation.rhs(1)=6.0;
        Vector old(2),current(2);
        old(0)=2.0;old(1)=1.5;
        current(0)=2.0;current(1)=1.5;
        const double alpha=0.5;
        // Apply the SAME equation-relaxation algebra as the production solver.
        const auto* row=equation.matrix.row_offsets_data();
        const auto* cols=equation.matrix.columns_data();
        auto* vals=equation.matrix.values_data();
        for (std::size_t i=0;i<2;++i) {
            const double added=equation.diagonal[i]*(1.0/alpha-1.0);
            for (auto k=row[i];k<row[i+1];++k)
                if(cols[k]==i)vals[k]+=added;
            equation.diagonal[i]+=added;
            equation.rhs(i)+=added*old(i);
        }
        // Construct the EXACT solution of the relaxed equation by 2x2 LU.
        const double d=8.0*12.0-1.0;
        current(0)=(12.0*equation.rhs(0)+equation.rhs(1))/d;
        current(1)=(equation.rhs(0)+8.0*equation.rhs(1))/d;
        const auto certificate=audit_momentum_linear_relaxation(equation,current,old,alpha);
        if (certificate.relaxed_true_relative>1e-14 ||
            !(certificate.original_operator_relative>1e-2) ||
            !(certificate.relaxation_displacement_l2>0.0))
            throw std::runtime_error("relaxed linear convergence hidden original-equation defect");
        auto bad=equation;
        bad.rhs(0)+=0.05;
        const auto corrupted=audit_momentum_linear_relaxation(bad,current,old,alpha);
        if (!(corrupted.relaxed_true_relative>1e-4))
            throw std::runtime_error("linear true residual oracle missed corrupted RHS");
        std::cout << "RELAXATION_TRUE_LINEAR_CERTIFICATE relaxed="
                  << certificate.relaxed_true_relative
                  << " original=" << certificate.original_operator_relative << "\n";
    }
    // At the exact SAME flux/pressure/predictor state, TVD-UPWIND must
    // equal the signed deferred convection and the Rhie-Chow decomposition
    // must be facewise closed. This is a diagnostic identity, not a proof
    // that either momentum or pressure has converged.
    {
        Field<double,Location::FACE> predictor(m.n_faces(),"predictor","kg/s",1);
        Field<double,Location::FACE> relaxed(m.n_faces(),"relaxed","kg/s",1);
        Field<double,Location::FACE> full(m.n_faces(),"full","kg/s",1);
        Field<double,Location::FACE> authoritative(m.n_faces(),"authoritative","kg/s",1);
        for (std::size_t face=0;face<m.n_faces();++face) {
            predictor(face)=flux(face);
            relaxed(face)=predictor(face)-0.03*static_cast<double>(face%5);
            full(face)=relaxed(face)+0.01*static_cast<double>(face%3);
            authoritative(face)=full(face)+0.002*static_cast<double>(face%2);
        }
        Field<double,Location::FACE> pressure_operator(m.n_faces(),"rc_ref","kg/s",1);
        for (std::size_t f=0;f<m.n_faces();++f)
            pressure_operator(f)=predictor(f)-full(f);
        const auto trace=audit_momentum_same_state_faces(
            m,g,authoritative,phi,momentum_bcs,ConvectionScheme::TVD,1.0,
            predictor,relaxed,full,&pressure_operator);
        const auto tvd=reconstruct_scalar_transport_flux(
            m,g,authoritative,phi,0.0,momentum_bcs,true,ConvectionScheme::TVD);
        const auto upwind=reconstruct_scalar_transport_flux(
            m,g,authoritative,phi,0.0,momentum_bcs,true,ConvectionScheme::UPWIND);
        for (const auto& face:trace.faces) {
            if (std::abs(face.deferred_tvd_flux-
                (tvd(face.face)-upwind(face.face)))>2e-13 ||
                std::abs(face.decomposition_error)>1e-14 ||
                (face.neighbour>=0 && std::abs(face.pressure_operator_error)>1e-14))
                throw std::runtime_error("same-state TVD/RC face audit failed");
        }
        if (trace.faces.size()!=m.n_faces() ||
            !(trace.max_deferred_tvd_flux>1e-5))
            throw std::runtime_error("same-state trace missed nonzero TVD correction");
        std::cout << "SAME_STATE_FACE_TVD_RC: PASS max_deferred="
                  << trace.max_deferred_tvd_flux << "\n";
    }
    // The Picard assembly uses U_old and phi_old, not the final corrected
    // velocity/flux. Test partial continuation and nonzero deferred flux.
    for (double continuation : {0.0,0.5,1.0}) {
        const auto trace=audit_momentum_picard_tvd_faces(
            m,g,flux,phi,momentum_bcs,ConvectionScheme::TVD,continuation);
        const auto expected=reconstruct_scalar_transport_flux(
            m,g,flux,phi,0.0,momentum_bcs,true,ConvectionScheme::TVD,
            nullptr,nullptr,nullptr,nullptr,continuation);
        if (trace.faces.size()!=m.n_faces() ||
            std::abs(trace.high_order_fraction-continuation)>1e-14)
            throw std::runtime_error("Picard TVD trace size or continuation mismatch");
        for (const auto& face:trace.faces)
            if (std::abs(face.selected_flux-expected(face.face))>2e-13 ||
                (continuation==0.0 && std::abs(face.deferred_tvd_flux)>2e-13))
                throw std::runtime_error("Picard TVD did not match frozen operator");
    }
    std::cout << "PICARD_DEFERRED_TVD: PASS\n";
    std::cout << "TVD_SKEW_AFFINE_MMS: PASS\n";
}
