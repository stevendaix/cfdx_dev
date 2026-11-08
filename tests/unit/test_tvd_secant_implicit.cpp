#include "cfdx/physics/tvd_secant_implicit_correction.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
static void require(bool pass, const char* reason) {
    if (!pass) throw std::runtime_error(reason);
}
static double residual(const ScalarEquation& e, std::size_t cell,
                       const std::vector<double>& u) {
    const auto* offsets=e.matrix.row_offsets_data();
    const auto* col=e.matrix.columns_data();
    const auto* vals=e.matrix.values_data();
    double r=-e.rhs(cell);
    for (auto k=offsets[cell]; k<offsets[cell+1];++k)
        r+=vals[k]*u[col[k]];
    return r;
}
int main() {
    Mesh mesh;
    mesh.points().resize(4);
    mesh.faces().push_face({0,1,2,3});
    mesh.ownership().resize(1);
    mesh.ownership().set_owner(0,0);
    mesh.ownership().set_neighbour(0,1);
    mesh.cells().push_cell({0});
    mesh.cells().push_cell({0});
    ScalarEquation eq;
    eq.matrix=SparseMatrix(2,2);
    eq.matrix.push_back(0,0,10.0);
    eq.matrix.push_back(0,1,-1.0);
    eq.matrix.push_back(1,0,-1.0);
    eq.matrix.push_back(1,1,10.0);
    eq.matrix.finalize();
    eq.rhs=Vector(2,0.0);
    eq.rhs(0)=-3.0; eq.rhs(1)=3.0;
    eq.diagonal={10.0,10.0};
    Field<double,Location::FACE> face(1,"mflux","kg/s",0.0);
    face(0)=4.0;
    const std::vector<double> u{2.0,4.0}, old{1.0,5.0};
    const std::vector<double> raw{3.0}, prev{1.0};
    const auto percell=tvd_deferred_source_by_cell(mesh,raw);
    require(percell.size()==2 && percell[0]==-3.0 && percell[1]==3.0,
            "Deferred face source must be exactly conservative");
    const double before0=residual(eq,0,u);
    const double before1=residual(eq,1,u);
    TvdSecantImplicitControls controls;
    controls.enabled=true;
    controls.spectral_radius_cap=0.5;
    const auto out=add_tvd_secant_implicit_diagonal(
        mesh,eq,raw,prev,u,old,face,controls);
    require(out.regularized_cells==2, "secant should stabilize both cells");
    require(std::abs(out.max_diagonal_added-2.0)<1e-14,
            "secant diagonal unexpected or cap ignored");
    require(std::abs(eq.diagonal[0]-12.0)<1e-14 &&
            std::abs(eq.diagonal[1]-12.0)<1e-14,
            "ScalarEquation diagonal must match assembled CSR");
    require(std::abs(residual(eq,0,u)-before0)<1e-14 &&
            std::abs(residual(eq,1,u)-before1)<1e-14,
            "secant modification must cancel exactly at frozen input");
    require(std::abs(eq.rhs(0)-1.0)<1e-14 &&
            std::abs(eq.rhs(1)-11.0)<1e-14,
            "secant RHS shift wrong sign");
    // Positive dB/dU is antidiffusive and must NOT become a negative diagonal.
    ScalarEquation negative;
    negative.matrix=SparseMatrix(2,2);
    negative.matrix.push_back(0,0,10.0);
    negative.matrix.push_back(1,1,10.0);
    negative.matrix.finalize();
    negative.rhs=Vector(2,0.0);
    negative.diagonal={10.0,10.0};
    const auto skip=add_tvd_secant_implicit_diagonal(
        mesh,negative,prev,raw,u,old,face,controls);
    require(skip.regularized_cells==0,
            "positive derivative must not induce a negative stabilizer");
    // At the fixed point, tiny changes never cause 1/epsilon secant growth.
    const std::vector<double> almost{2.0+1e-14,4.0};
    ScalarEquation tiny=negative;
    const auto zero=add_tvd_secant_implicit_diagonal(
        mesh,tiny,raw,prev,almost,u,face,controls);
    require(zero.regularized_cells==0,"small delta needs a no-secant guard");
    bool threw=false;
    try {controls.spectral_radius_cap=-0.01; validate_tvd_secant_implicit(controls);}
    catch(const std::invalid_argument&) {threw=true;}
    require(threw,"invalid spectral cap must be rejected");
    threw=false;
    try {controls.spectral_radius_cap=0.5;
         (void)add_tvd_secant_implicit_diagonal(mesh,tiny,{3.},{1.},{2.},{1.},face,controls);}
    catch(const std::invalid_argument&) {threw=true;}
    require(threw,"inconsistent source dimensions must be rejected");
    std::cout << "PASS exact frozen-state residual, conservative TVD source, "
                 "diagonal secant stabilization, guards\n";
}

