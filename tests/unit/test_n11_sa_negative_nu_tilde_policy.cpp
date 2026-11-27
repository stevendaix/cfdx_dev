// N11 F1 (#480): route the existing NegativeNuTildePolicy knob through the
// Spalart-Allmaras transport path. CLAMP_ZERO is the historical default and
// must stay unchanged; REJECT is the fail-closed qualification mode.
#include "cfdx/physics/turbulence_transport.h"
#include "cfdx/physics/turbulence_solver.h"

#include "common/test_harness.h"

#include <cmath>
#include <limits>
#include <stdexcept>

using namespace cfdx::testing;

namespace {

cfdx::physics::TurbulenceTransportControls sa_controls(
    cfdx::physics::NegativeNuTildePolicy policy)
{
    cfdx::physics::TurbulenceTransportControls c;
    c.model = cfdx::physics::TurbulenceModel::SPALART_ALLMARAS;
    c.molecular_viscosity = 1.0e-5;
    c.negative_nu_tilde_policy = policy;
    return c;
}

// Independent transcription of fiche 09, section 2 (SA):
// nu_t = nu_tilde * fv1, fv1 = chi^3/(chi^3 + cv1^3), chi = nu_tilde/nu.
double sa_reference_formula(double nu_tilde, double molecular_nu)
{
    const double chi = nu_tilde / molecular_nu;
    const double chi3 = chi * chi * chi;
    const double cv1 = 7.1;
    return nu_tilde * chi3 / (chi3 + cv1 * cv1 * cv1);
}

void default_policy_is_clamp_zero()
{
    cfdx::physics::TurbulenceTransportControls c;
    EXPECT_TRUE(c.negative_nu_tilde_policy ==
                cfdx::physics::NegativeNuTildePolicy::CLAMP_ZERO);
}

void clamp_zero_clamps_negative_nu_tilde()
{
    const auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::CLAMP_ZERO);
    const double nu = cfdx::physics::turbulence_nu_t(0.0, -1.0e-4, 0.0, 1.0, c);
    EXPECT_TRUE(std::isfinite(nu));
    EXPECT_NEAR(nu, 0.0, 1e-15);
}

void clamp_zero_keeps_positive_formula()
{
    const auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::CLAMP_ZERO);
    const double nu = cfdx::physics::turbulence_nu_t(0.0, 2.0e-5, 0.0, 1.0, c);
    EXPECT_NEAR(nu, sa_reference_formula(2.0e-5, 1.0e-5), 1e-18);
}

void reject_throws_domain_error_on_negative()
{
    const auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::REJECT);
    EXPECT_THROW_WITH(
        cfdx::physics::turbulence_nu_t(0.0, -1.0e-4, 0.0, 1.0, c),
        std::domain_error, "Spalart-Allmaras negative nu_tilde");
}

void reject_throws_invalid_argument_on_nonfinite()
{
    const auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::REJECT);
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW_WITH(
        cfdx::physics::turbulence_nu_t(0.0, nan_value, 0.0, 1.0, c),
        std::invalid_argument, "Spalart-Allmaras invalid inputs");
}

void reject_keeps_positive_path_unchanged()
{
    const auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::REJECT);
    const double nu = cfdx::physics::turbulence_nu_t(0.0, 2.0e-5, 0.0, 1.0, c);
    EXPECT_NEAR(nu, sa_reference_formula(2.0e-5, 1.0e-5), 1e-18);
}

cfdx::core::Mesh make_sa_transport_test_mesh()
{
    using namespace cfdx::core;
    Mesh mesh;
    mesh.points().resize(8);
    const double points[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i=0;i<8;++i)
        mesh.points().set(i,points[i][0],points[i][1],points[i][2]);
    const std::vector<std::vector<FaceIndex>> faces = {
        {0,3,2,1},{4,5,6,7},{0,1,5,4},
        {3,7,6,2},{0,4,7,3},{1,2,6,5}
    };
    for (const auto& face : faces) mesh.faces().push_face(face);
    mesh.ownership().resize(faces.size());
    for (std::size_t face=0;face<faces.size();++face) {
        mesh.ownership().set_owner(face,0);
        mesh.ownership().set_neighbour(face,FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0,1,2,3,4,5});
    Patch walls;
    walls.name="walls";
    walls.type=PatchType::WALL;
    walls.face_ids={0,1,2,3,4,5};
    mesh.boundary().add_patch(walls);
    return mesh;
}

void transport_solver_enforces_reject_policy_on_input_state()
{
    using namespace cfdx::core;
    using namespace cfdx::physics;
    const auto mesh = make_sa_transport_test_mesh();
    const auto geometry = build_fv_geometry(mesh);
    Field<double,Location::FACE> flux(mesh.n_faces(),"phi","kg/s",1);
    flux.fill(0.0);
    Field<double,Location::CELL> nu_tilde(mesh.n_cells(),"nu_tilde","m2/s",1);
    Field<double,Location::CELL> vorticity(mesh.n_cells(),"vorticity","1/s",1);
    Field<double,Location::CELL> wall_distance(mesh.n_cells(),"wall_distance","m",1);
    nu_tilde(0)=-1.0e-4;
    vorticity(0)=0.0;
    wall_distance(0)=0.5;
    auto controls=sa_controls(NegativeNuTildePolicy::REJECT);
    EXPECT_THROW_WITH(
        solve_spalart_allmaras_transport(
            mesh,geometry,flux,nu_tilde,vorticity,wall_distance,controls,{},1,1e-8),
        std::domain_error,"negative input rejected by policy");
    EXPECT_NEAR(nu_tilde(0),-1.0e-4,0.0);
}

void molecular_viscosity_guard_intact()
{
    auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::CLAMP_ZERO);
    c.molecular_viscosity = 0.0;
    EXPECT_THROW_WITH(
        cfdx::physics::turbulence_nu_t(0.0, 2.0e-5, 0.0, 1.0, c),
        std::invalid_argument, "SA eddy viscosity requires molecular viscosity");
}

void consistent_with_audited_closure_wrapper()
{
    const auto c = sa_controls(cfdx::physics::NegativeNuTildePolicy::CLAMP_ZERO);
    const double transport = cfdx::physics::turbulence_nu_t(0.0, 2.0e-5, 0.0, 1.0, c);
    const double wrapper = cfdx::physics::spalart_allmaras_nu_t(
        2.0e-5, 1.0, 1.0e-5, cfdx::physics::NegativeNuTildePolicy::CLAMP_ZERO);
    EXPECT_NEAR(transport, wrapper, 1e-15);
}

} // namespace

int main()
{
    run_case("SA transport: default policy is CLAMP_ZERO (historical)", default_policy_is_clamp_zero);
    run_case("SA transport: CLAMP_ZERO clamps negative nu_tilde to zero", clamp_zero_clamps_negative_nu_tilde);
    run_case("SA transport: CLAMP_ZERO keeps the positive formula unchanged", clamp_zero_keeps_positive_formula);
    run_case("SA transport: REJECT throws domain_error on negative nu_tilde", reject_throws_domain_error_on_negative);
    run_case("SA transport: REJECT throws invalid_argument on non-finite nu_tilde", reject_throws_invalid_argument_on_nonfinite);
    run_case("SA transport: REJECT keeps the positive path unchanged", reject_keeps_positive_path_unchanged);
    run_case("SA transport: molecular viscosity guard intact", molecular_viscosity_guard_intact);
    run_case("SA transport solver: REJECT refuses negative input transactionally", transport_solver_enforces_reject_policy_on_input_state);
    run_case("SA transport: consistent with the audited closure wrapper", consistent_with_audited_closure_wrapper);
    return run_all();
}
