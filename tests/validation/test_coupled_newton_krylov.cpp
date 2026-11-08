#include "cfdx/physics/coupled_newton_krylov.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::core;

static void require(bool ok,const char* label) {
    if (!ok) throw std::runtime_error(label);
}

// Two-cell manufactured system: 3 momentum components, p/gauge,
// continuity, and three independent Rhie-Chow face constitutive equations.
// An exact Jacobian here necessarily differentiates U, p, phi, and H/A.
static std::vector<double> two_cell_coupled(const std::vector<double>& x) {
    require(x.size()==11,"full coupled state dimension");
    const double u0=x[0],u1=x[1],v0=x[2],v1=x[3];
    const double w0=x[4],w1=x[5],p0=x[6],p1=x[7];
    const double f0=x[8],f1=x[9],f2=x[10];
    // Reassemble a nonlinear H/A predictor (not frozen during Jv).
    const double h0=0.85*u0+0.12*v0*v0+0.04*f1*u1;
    const double h1=0.9*u1-0.07*w1*w1-0.03*f1*u0;
    const double rc0=-h0;
    const double rc1=0.5*(h0+h1)-0.3*(p1-p0);
    const double rc2=h1;
    return {
        1.8*u0-0.3*u1+0.1*f1*u0+(p1-p0),
        1.5*u1-0.2*u0+0.2*f1*u1-(p1-p0),
        1.1*v0-0.08*v1+0.1*f1*v0+0.2*(p1-p0),
        1.3*v1-0.09*v0-0.1*f1*v1-0.2*(p1-p0),
        1.2*w0-0.1*w1+0.05*f1*w0,
        1.4*w1-0.08*w0+0.03*f1*w1,
        p0,                      // pressure-gauge row in place of redundant continuity
        f2-f1,                  // true cell-1 mass balance
        f0-rc0, f1-rc1, f2-rc2  // constitutive constraints on all three faces
    };
}

static Mesh skew_mesh() {
    Mesh m;
    constexpr std::size_t nx=2,ny=2,nz=1;
    const auto v=[](std::size_t i,std::size_t j,std::size_t k) {
        return (k*3+j)*3+i;
    };
    m.points().resize(3*3*2);
    for(std::size_t k=0;k<=nz;++k)
        for(std::size_t j=0;j<=ny;++j)
            for(std::size_t i=0;i<=nx;++i)
                m.points().set(v(i,j,k),double(i)+0.35*double(j),double(j),double(k));
    std::map<std::array<std::size_t,4>,std::size_t> seen;
    for(std::size_t j=0;j<ny;++j)
        for(std::size_t i=0;i<nx;++i) {
            const auto a=v(i,j,0),b=v(i+1,j,0),c=v(i+1,j+1,0),d=v(i,j+1,0);
            const auto A=v(i,j,1),B=v(i+1,j,1),C=v(i+1,j+1,1),D=v(i,j+1,1);
            std::vector<std::size_t> faces;
            const auto add=[&](std::array<std::size_t,4> vertices) {
                auto key=vertices;std::sort(key.begin(),key.end());
                const auto it=seen.find(key);
                if(it==seen.end()) {
                    const auto f=m.n_faces();
                    m.faces().push_face({vertices.begin(),vertices.end()});
                    m.ownership().resize(m.n_faces());
                    m.ownership().set_owner(f,m.n_cells());
                    m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
                    seen.emplace(key,f);faces.push_back(f);
                } else {
                    m.ownership().set_neighbour(it->second,static_cast<std::int64_t>(m.n_cells()));
                    faces.push_back(it->second);
                }
            };
            add({a,d,c,b});add({A,B,C,D});add({a,b,B,A});add({d,D,C,c});
            add({a,A,D,d});add({b,c,C,B});
            m.cells().push_cell(faces);
        }
    m.boundary().add_patch("wall",PatchType::WALL);
    m.boundary().add_patch("frontBack",PatchType::EMPTY);
    for(std::size_t f=0;f<m.n_faces();++f) {
        if(m.ownership().neighbour(f)>=0) continue;
        const auto start=m.faces().offsets_data()[f];
        const auto end=m.faces().offsets_data()[f+1];
        const auto zfirst=m.faces().vertices_data()[start]/9;
        bool empty=true;
        for(auto k=start;k<end;++k)
            empty=empty && (m.faces().vertices_data()[k]/9==zfirst);
        m.boundary().patch(empty?1:0).face_ids.push_back(f);
    }
    return m;
}

int main() {
    // Analytic Jv including the pressure gauge, longitudinal/transverse
    // momentum, mass balance, and Rhie-Chow H/A response on every face.
    const std::vector<double> x{0.17,0.25,-0.11,0.20,0.04,-0.08,0.0,0.06,
                                -0.16,0.21,0.19};
    const std::vector<double> v{0.2,-0.12,0.31,0.07,0.04,-0.16,0.0,-0.25,
                                0.08,-0.11,0.03};
    const auto jv=coupled_newton_directional_jvp(two_cell_coupled,x,v);
    const auto check_jv=coupled_newton_directional_jvp(two_cell_coupled,x,v,3e-6);
    double agreement=0.0;
    for(std::size_t k=0;k<x.size();++k)
        agreement=std::max(agreement,std::abs(jv[k]-check_jv[k]));
    require(agreement<1e-8,"coupled Jv not stable to the finite-difference step");
    auto p_only=std::vector<double>(x.size(),0.0);
    p_only[7]=1.0;
    const auto jp=coupled_newton_directional_jvp(two_cell_coupled,x,p_only);
    require(std::abs(jp[0])>0.1 && std::abs(jp[9])>0.1,
            "pressure does not enter momentum and Rhie-Chow constraint");
    auto velocity_only=std::vector<double>(x.size(),0.0);
    velocity_only[0]=1.0;
    const auto ju=coupled_newton_directional_jvp(two_cell_coupled,x,velocity_only);
    require(std::abs(ju[9])>0.1,"H/A velocity derivative was frozen");
    auto phi_only=std::vector<double>(x.size(),0.0);
    phi_only[9]=1.0;
    const auto jf=coupled_newton_directional_jvp(two_cell_coupled,x,phi_only);
    require(std::abs(jf[0])>1e-3 && std::abs(jf[7])>0.1,
            "convection / continuity face-flux derivative was frozen");
    CoupledNewtonKrylovControls controls;
    controls.krylov_dimension=11;
    controls.linear_relative_tolerance=1e-5;
    const auto step=coupled_newton_krylov_step(two_cell_coupled,x,
                                              std::vector<double>(11,1.0),controls);
    require(step.accepted,"fully coupled Newton-Krylov step rejected on manufactured system");
    require(step.physical_merit_after<0.1*step.physical_merit_before,
            "full coupled Newton step does not reduce true residual");
    require(step.linear_true_relative_residual<1e-5,
            "GMRES returned an invalid true directional residual");
    const auto rejected=coupled_newton_krylov_step(two_cell_coupled,x,
                                                   std::vector<double>(11,1.0),controls,
                                                   [](const auto&){return false;});
    require(!rejected.accepted && rejected.state==x,
            "rejected Newton step mutated accepted state");

    // SECOND oracle: the *actual* multidimensional finite-volume MUSCL
    // assembly, pressure gradient, H/A and Rhie-Chow mass flux, all reassembled
    // at every (U,p,phi) perturbation on a skewed grid. This prevents a frozen
    // transport/pressure/flux derivative from passing the coupled-Jv check.
    auto mesh=skew_mesh();
    const auto geo=build_fv_geometry(mesh);
    const auto nc=mesh.n_cells(),nf=mesh.n_faces();
    VelocityBoundaryConditions ubc{
        {"wall",{VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}}},
        {"frontBack",{VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}}}};
    ScalarBoundaryConditions pbc{
        {"wall",{ScalarBoundaryType::ZERO_GRADIENT,0,0}},
        {"frontBack",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}};
    const std::array<ScalarBoundaryConditions,3> ubcs{{
        {{"wall",{ScalarBoundaryType::FIXED_VALUE,0,0}},
         {"frontBack",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}},
        {{"wall",{ScalarBoundaryType::FIXED_VALUE,0,0}},
         {"frontBack",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}},
        {{"wall",{ScalarBoundaryType::FIXED_VALUE,0,0}},
         {"frontBack",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}}}};
    const std::vector<double> mu(nc,0.1);
    const CoupledResidual finite_volume_full=[&](const std::vector<double>& q) {
        require(q.size()==4*nc+nf,"wrong FV coupled vector dimensions");
        Field<double,Location::CELL> U(nc,"u","m/s",3);
        Field<double,Location::CELL> p(nc,"p","Pa",1);
        Field<double,Location::FACE> phi(nf,"phi","kg/s",1);
        Field<double,Location::CELL> zeros(nc,"zeros","N/m3",1);
        zeros.fill(0.0);
        for(std::size_t d=0;d<3;++d)
            for(std::size_t c=0;c<nc;++c)
                U.component_data(d)[c]=q[d*nc+c];
        for(std::size_t c=0;c<nc;++c) p(c)=q[3*nc+c];
        for(std::size_t f=0;f<nf;++f) phi(f)=q[4*nc+f];
        const auto gp=gauss_gradient_with_boundary(p,mesh,geo,pbc,false,false);
        std::array<ScalarEquation,3> eq;
        std::array<std::vector<double>,3> inverse,ha;
        std::vector<double> residual(q.size(),0.0);
        for(std::size_t d=0;d<3;++d) {
            Field<double,Location::CELL> ui(nc,"ui","m/s",1);
            for(std::size_t c=0;c<nc;++c) ui(c)=q[d*nc+c];
            eq[d]=assemble_momentum_component(mesh,geo,phi,gp,zeros,mu,ubcs[d],d,
                       true,ConvectionScheme::MUSCL_CELL_LIMITED,&ui);
            inverse[d].resize(nc);
            ha[d].resize(nc);
            for(std::size_t c=0;c<nc;++c) {
                const double diag=eq[d].diagonal[c];
                require(diag>0.0,"FV Jacobian invalid momentum diagonal");
                inverse[d][c]=1.0/diag;
                double h=eq[d].rhs(c)+gp.component_data(d)[c]*geo.cell_volumes[c];
                for(auto k=eq[d].matrix.row_offsets_data()[c];
                    k<eq[d].matrix.row_offsets_data()[c+1];++k) {
                    const auto col=eq[d].matrix.columns_data()[k];
                    if(col!=c) h-=eq[d].matrix.values_data()[k]*ui(col);
                    residual[d*nc+c]+=eq[d].matrix.values_data()[k]*ui(col);
                }
                ha[d][c]=h*inverse[d][c];
                residual[d*nc+c]-=eq[d].rhs(c);
            }
        }
        Field<double,Location::CELL> H(nc,"HbyA","m/s",3);
        for(std::size_t d=0;d<3;++d)
            for(std::size_t c=0;c<nc;++c) H.component_data(d)[c]=ha[d][c];
        const auto rc=make_rhie_chow_mass_flux(mesh,geo,H,p,inverse,1.0,ubc,pbc,
                        RhieChowVelocityContract::HBYA_PREDICTOR,false,false);
        for(std::size_t f=0;f<nf;++f) {
            const auto o=mesh.ownership().owner(f);
            residual[3*nc+o]+=phi(f);
            const auto neighbour=mesh.ownership().neighbour(f);
            if(neighbour>=0)
                residual[3*nc+static_cast<std::size_t>(neighbour)]-=phi(f);
            residual[4*nc+f]=phi(f)-rc(f);
        }
        // Pure-Neumann pressure: one redundant continuity row becomes gauge.
        residual[3*nc]=p(0);
        return residual;
    };
    std::vector<double> state(4*nc+nf,0.0),direction(4*nc+nf,0.0);
    for(std::size_t i=0;i<state.size();++i) {
        state[i]=0.02*std::sin(0.81*double(i)+0.5);
        direction[i]=0.04*std::cos(0.51*double(i)+0.2);
    }
    state[3*nc]=0.0;
    const auto fv_jv=coupled_newton_directional_jvp(finite_volume_full,state,direction,5e-7);
    const auto fv_jv2=coupled_newton_directional_jvp(finite_volume_full,state,direction,1e-6);
    double error=0.0;
    for(std::size_t k=0;k<state.size();++k)
        error=std::max(error,std::abs(fv_jv[k]-fv_jv2[k]));
    require(error<2e-7,"complete finite-volume coupled Jv has a step-dependent defect");
    auto dp=std::vector<double>(state.size(),0.0);
    dp[3*nc+1]=1.0;
    const auto pressure_response=coupled_newton_directional_jvp(finite_volume_full,state,dp);
    double momentum_pressure_response=0.0,rc_pressure_response=0.0;
    for(std::size_t i=0;i<3*nc;++i)
        momentum_pressure_response=std::max(momentum_pressure_response,std::abs(pressure_response[i]));
    for(std::size_t i=4*nc;i<state.size();++i)
        rc_pressure_response=std::max(rc_pressure_response,std::abs(pressure_response[i]));
    require(momentum_pressure_response>1e-5 && rc_pressure_response>1e-5,
            "real FV coupled Jv dropped pressure-gradient or Rhie-Chow derivative");
    std::cout << "COUPLED_JFNK_PASS toy_merit="<<step.physical_merit_before
              <<" -> "<<step.physical_merit_after
              <<" GMRES_true="<<step.linear_true_relative_residual
              <<" FV_Jv_step_error="<<error
              <<" FV_pressure_momentum="<<momentum_pressure_response
              <<" FV_pressure_RC="<<rc_pressure_response<<'\n';
}
