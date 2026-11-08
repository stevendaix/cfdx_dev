// Multidimensional TVD affine-skew MMS and matrix/face-flux identity.
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/nonlinear_multisecant_acceleration.h"
#include "cfdx/physics/momentum_face_budget.h"
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

    // Pairwise TVD limiters are bounded but NOT affine exact on skew faces.
    // The multidimensional MUSCL correction is the affine-exact alternative;
    // never mistake a test of the original TVD method for evidence that it
    // also achieves that stronger property. All schemes must remain
    // assembly/face-flux consistent regardless of their truncation error.
    for (const auto scheme : {ConvectionScheme::TVD,
                              ConvectionScheme::TVD_VANLEER,
                              ConvectionScheme::TVD_VANALBADA,
                              ConvectionScheme::MUSCL_CELL_LIMITED}) {
        const auto transported=reconstruct_scalar_transport_flux(
            m,g,flux,phi,0.0,{},true,scheme,&bcs);
        double face_error=0.0;
        for (std::size_t f=0;f<m.n_faces();++f)
            if (m.ownership().neighbour(f)>=0)
                face_error=std::max(face_error,
                    std::abs(transported(f)-g.face_centres[f].y));
        if (scheme == ConvectionScheme::MUSCL_CELL_LIMITED && face_error>1e-11)
            throw std::runtime_error("multidimensional MUSCL lost affine exactness on a skew mesh");

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
    // Directional Jacobian-vector product of the COMPLETE nonlinear TVD
    // momentum/transport discretization (limiter gradients, all face fluxes,
    // bounded convection term), here at frozen pressure and mass flux.
    // Compare the production A(u)u-b(u) to an INDEPENDENT face-divergence
    // evaluator at the two perturbed states; a frozen UPWIND derivative must
    // NOT pass this test merely because the CSR matrix is solved accurately.
    {
        const std::size_t n=m.n_cells();
        std::vector<double> state(n), direction(n);
        for (std::size_t c=0;c<n;++c) {
            state[c]=0.42+0.21*std::sin(0.73*static_cast<double>(c))+
                     0.017*static_cast<double>(c);
            direction[c]=0.15*std::cos(0.61*static_cast<double>(c))+
                         0.018*static_cast<double>(c);
        }
        const std::vector<double> diag(n,4.0);
        const auto matrix_residual=[&](const std::vector<double>& values,
                                      ConvectionScheme scheme) {
            Field<double,Location::CELL> u(n,"directional_u","",1);
            for(std::size_t c=0;c<n;++c) u(c)=values[c];
            const auto eq=assemble_scalar_equation(
                m,g,flux,0.0,zero,zero,{},true,nullptr,
                &diag,nullptr,nullptr,scheme,&u,&bcs,false);
            std::vector<double> res(n);
            for(std::size_t row=0;row<n;++row) {
                double value=-eq.rhs(row);
                for(auto k=eq.matrix.row_offsets_data()[row];
                    k<eq.matrix.row_offsets_data()[row+1];++k)
                    value+=eq.matrix.values_data()[k]*
                        values[eq.matrix.columns_data()[k]];
                res[row]=value;
            }
            return res;
        };
        const auto face_residual=[&](const std::vector<double>& values,
                                    ConvectionScheme scheme) {
            Field<double,Location::CELL> u(n,"directional_u","",1);
            for(std::size_t c=0;c<n;++c) u(c)=values[c];
            const auto transported=reconstruct_scalar_transport_flux(
                m,g,flux,u,0.0,{},true,scheme,&bcs);
            std::vector<double> r(n,0.0), div_mass(n,0.0);
            for(std::size_t f=0;f<m.n_faces();++f) {
                const auto o=m.ownership().owner(f);
                r[o]+=transported(f);
                div_mass[o]+=flux(f);
                const auto nb=m.ownership().neighbour(f);
                if(nb>=0) {
                    const auto j=static_cast<std::size_t>(nb);
                    r[j]-=transported(f);
                    div_mass[j]-=flux(f);
                }
            }
            for(std::size_t c=0;c<n;++c)
                r[c]+=diag[c]*values[c]-values[c]*div_mass[c];
            return r;
        };
        for(auto scheme : {ConvectionScheme::TVD,
                            ConvectionScheme::TVD_VANLEER,
                            ConvectionScheme::TVD_VANALBADA}) {
            const auto matrix_jv=nonlinear_directional_jvp(
                [&](const auto& u){return matrix_residual(u,scheme);},
                state,direction);
            const auto face_jv=nonlinear_directional_jvp(
                [&](const auto& u){return face_residual(u,scheme);},
                state,direction);
            double error=0.0;
            for(std::size_t c=0;c<n;++c)
                error=std::max(error,std::abs(matrix_jv[c]-face_jv[c]));
            if(error>2e-8)
                throw std::runtime_error(
                    "full TVD direction Jv disagrees with independent face flux");
            const auto upwind=nonlinear_directional_jvp(
                [&](const auto& u){return matrix_residual(u,ConvectionScheme::UPWIND);},
                state,direction);
            double difference=0.0;
            for(std::size_t c=0;c<n;++c)
                difference=std::max(difference,std::abs(matrix_jv[c]-upwind[c]));
            if (difference<1e-5)
                throw std::runtime_error(
                    "TVD directional test is vacuous: equal to frozen UPWIND");
            std::cout<<"TVD_DIRECTIONAL_JV scheme="<<static_cast<int>(scheme)
                     <<" CSR_vs_faces="<<error
                     <<" nonlinear_vs_upwind="<<difference<<"\n";
        }
        // Directional Newton demonstration on the ACTUAL nonlinearly
        // reassembled TVD finite-volume operator (frozen mass flux and
        // pressure). This is deliberately not a claim of coupled NACA Newton.
        // Columns J*e_j come from directional reassemblies, NOT an UPWIND
        // matrix. A pivoted small dense solve is used only on this 18-cell
        // manufactured problem; production would need matrix-free GMRES.
        const auto residual=[&](const std::vector<double>& u) {
            return face_residual(u,ConvectionScheme::TVD);
        };
        const auto initial=residual(state);
        const auto infinity_norm=[](const auto& values) {
            double value=0;
            for(double x:values)value=std::max(value,std::abs(x));
            return value;
        };
        const double start=infinity_norm(initial);
        std::vector<std::vector<double>> J(n,std::vector<double>(n,0.0));
        for(std::size_t column=0;column<n;++column) {
            std::vector<double> direction_e(n,0.0);
            direction_e[column]=1.0;
            const auto product=nonlinear_directional_jvp(
                residual,state,direction_e);
            for(std::size_t row=0;row<n;++row) J[row][column]=product[row];
        }
        std::vector<double> delta(n);
        for(std::size_t i=0;i<n;++i)delta[i]=-initial[i];
        for(std::size_t j=0;j<n;++j) {
            std::size_t pivot=j;
            for(std::size_t i=j+1;i<n;++i)
                if(std::abs(J[i][j])>std::abs(J[pivot][j]))pivot=i;
            if(std::abs(J[pivot][j])<1e-12)
                throw std::runtime_error("directional TVD Jacobian is singular in MMS");
            std::swap(J[j],J[pivot]);std::swap(delta[j],delta[pivot]);
            for(std::size_t i=j+1;i<n;++i) {
                const double ratio=J[i][j]/J[j][j];
                for(std::size_t k=j;k<n;++k)
                    J[i][k]-=ratio*J[j][k];
                delta[i]-=ratio*delta[j];
            }
        }
        for(std::size_t rev=0;rev<n;++rev) {
            const std::size_t i=n-1-rev;
            for(std::size_t j=i+1;j<n;++j)delta[i]-=J[i][j]*delta[j];
            delta[i]/=J[i][i];
        }
        double after=std::numeric_limits<double>::infinity();
        double accepted_step=0.0;
        for(double alpha : {1.0,0.5,0.25,0.125,0.0625,0.03125}) {
            auto candidate=state;
            for(std::size_t i=0;i<n;++i)candidate[i]+=alpha*delta[i];
            const double candidate_norm=infinity_norm(residual(candidate));
            if(std::isfinite(candidate_norm) &&
               candidate_norm<start*(1.0-1e-4*alpha)) {
                accepted_step=alpha;
                after=candidate_norm;
                break;
            }
        }
        if(!(after<start) || accepted_step<=0.0)
            throw std::runtime_error("directional TVD Newton failed to decrease true residual");
        std::cout<<"TVD_DIRECTIONAL_NEWTON_MMS initial="<<start
                 <<" final="<<after<<" step="<<accepted_step<<"\n";
    }
    // NACA-forensic operator identity: TWO independently assembled paths:
    //   (1) signed face residual C + D + P - B - u div(phi) + dt term
    //   (2) production CSR Picard row A(U) U - b(U).
    // Exercise TVD on skew geometry, diffusion, pressure, force and
    // non-solenoidal authoritative mass flux. No convergence tolerances change.
    Field<double,Location::CELL> gradp(m.n_cells(),"gradp_y","Pa/m",1);
    Field<double,Location::CELL> body(m.n_cells(),"body_y","N/m3",1);
    Field<double,Location::CELL> source(m.n_cells(),"source_momentum","N/m3",1);
    std::vector<double> mu(m.n_cells());
    for (std::size_t c=0;c<m.n_cells();++c) {
        gradp(c)=0.27+0.02*static_cast<double>(c);
        body(c)=0.04;
        source(c)=body(c)-gradp(c);
        mu[c]=0.04+0.002*static_cast<double>(c);
    }
    ScalarBoundaryConditions momentum_bcs;
    momentum_bcs["wall"]={ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
    for (auto scheme : {ConvectionScheme::TVD, ConvectionScheme::TVD_VANLEER,
                        ConvectionScheme::TVD_VANALBADA}) {
        const auto equation=assemble_scalar_equation(
            m,g,flux,0.0,source,zero,momentum_bcs,true,nullptr,
            nullptr,nullptr,&mu,scheme,&phi,nullptr,true);
        auto budget=diagnose_momentum_face_budget(
            m,g,flux,phi,gradp,body,mu,momentum_bcs,
            scheme,1.0,true);
        reconcile_momentum_budget_with_matrix(budget,equation,phi);
        if (!momentum_matrix_face_closure_consistent(budget))
            throw std::runtime_error("TVD accepted face budget does not match CSR operator");
        if (budget.matrix_closure_linf>3e-10)
            throw std::runtime_error("TVD matrix/face closure has a significant error");
        std::cout << "MOMENTUM_MATRIX_FACE_CLOSURE scheme=" << static_cast<int>(scheme)
                  << " linf=" << budget.matrix_closure_linf
                  << " relative=" << budget.matrix_closure_relative_linf << "\n";
        // A corrupted CELL away from the worst one MUST fail even when
        // the reported global worst physical residual remains unchanged.
        auto altered = budget;
        const std::size_t altered_cell=altered.worst_cell==0 ? 1 : 0;
        altered.cells[altered_cell].bounded_operator_residual += 0.02;
        reconcile_momentum_budget_with_matrix(altered,equation,phi);
        if (momentum_matrix_face_closure_consistent(altered))
            throw std::runtime_error("cellwise oracle missed displaced residual defect");
    }
    // A physical BE time derivative belongs in the face budget and in A*U-b.
    // A pseudo time diagonal must not be silently inserted in either physical
    // balance. This test inserts one real BE step and checks the sign.
    {
        std::vector<double> mass_diag(m.n_cells());
        std::vector<double> mass_rhs(m.n_cells());
        Field<double,Location::CELL> neg_accum(m.n_cells(),"minus_dUdt","N",1);
        for (std::size_t c=0;c<m.n_cells();++c) {
            mass_diag[c]=1.5*g.cell_volumes[c]/0.25;
            mass_rhs[c]=mass_diag[c]*(phi(c)-0.03);
            neg_accum(c)=-(mass_diag[c]*phi(c)-mass_rhs[c]);
        }
        const auto eq=assemble_scalar_equation(
            m,g,flux,0.0,source,zero,momentum_bcs,true,nullptr,
            &mass_diag,&mass_rhs,&mu,ConvectionScheme::TVD,&phi,nullptr,true);
        auto budget=diagnose_momentum_face_budget(
            m,g,flux,phi,gradp,body,mu,momentum_bcs,
            ConvectionScheme::TVD,1.0,true,nullptr,nullptr,nullptr,&neg_accum);
        reconcile_momentum_budget_with_matrix(budget,eq,phi);
        if (!momentum_matrix_face_closure_consistent(budget))
            throw std::runtime_error("time derivative signs disagree with A*U-b");
    }
    std::cout << "MOMENTUM_MATRIX_FACE_CLOSURE: PASS\n";
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
    // Boundary convective flux is zero in this MMS. The scalar boundary
    // map is intentionally empty; test face values are supplied through
    // bcs above, while the same-state oracle uses this patch-level map.
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
