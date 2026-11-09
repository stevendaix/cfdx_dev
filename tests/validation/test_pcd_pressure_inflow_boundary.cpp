// Pressure-space PCD inflow BC: auxiliary Kp/Fp only, default disabled.
// A two-cell FV mesh verifies face-flux sign (including local backflow)
// rather than patch label and ensures the legacy assembly is unchanged.
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh two_cell_mesh() {
    Mesh m;
    m.points().resize(12);
    const double p[12][3]={{0,0,0},{1,0,0},{2,0,0},{0,1,0},{1,1,0},{2,1,0},
                           {0,0,1},{1,0,1},{2,0,1},{0,1,1},{1,1,1},{2,1,1}};
    for (std::size_t i=0; i<12; ++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,6,9,3});
    m.faces().push_face({0,1,7,6});
    m.faces().push_face({3,9,10,4});
    m.faces().push_face({0,3,4,1});
    m.faces().push_face({6,7,10,9});
    m.faces().push_face({7,10,4,1});
    m.faces().push_face({2,5,11,8});
    m.faces().push_face({1,2,8,7});
    m.faces().push_face({4,10,11,5});
    m.faces().push_face({1,4,5,2});
    m.faces().push_face({7,8,11,10});
    m.ownership().resize(11);
    for(std::size_t f=0; f<11; ++f) {
        m.ownership().set_owner(f, f<=5 ? 0 : 1);
        m.ownership().set_neighbour(f, f==5 ? 1 : FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    m.cells().push_cell({5,6,7,8,9,10});
    Patch inlet; inlet.name="inlet"; inlet.type=PatchType::INLET;
    inlet.face_ids={0};
    Patch outlet; outlet.name="outlet"; outlet.type=PatchType::OUTLET;
    outlet.face_ids={6};
    Patch wall; wall.name="walls"; wall.type=PatchType::WALL;
    wall.face_ids={1,2,3,4,7,8,9,10};
    m.boundary().add_patch(inlet);
    m.boundary().add_patch(outlet);
    m.boundary().add_patch(wall);
    return m;
}

static double diagonal(const SparseMatrix& A, std::size_t row) {
    for (std::size_t k=A.row_offsets_data()[row]; k<A.row_offsets_data()[row+1]; ++k)
        if (A.columns_data()[k]==row) return A.values_data()[k];
    return 0.0;
}

int main() {
    run_case("pcd_inflow_dirichlet_is_opt_in_and_tracks_signed_face_flux", []{
        const Mesh mesh=two_cell_mesh();
        const FvGeometry geometry=build_fv_geometry(mesh);
        Field<double,Location::CELL> U(mesh.n_cells(),"U","m/s",3);
        ScalarBoundaryConditions pbc{
            {"inlet",{ScalarBoundaryType::ZERO_GRADIENT,0,0}},
            {"outlet",{ScalarBoundaryType::ZERO_GRADIENT,0,0}},
            {"walls",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}};
        for (double flow : {1.0,-1.0}) {
            for (std::size_t c=0;c<mesh.n_cells();++c) {
                U.component_data(0)[c]=flow;
                U.component_data(1)[c]=0.0;
                U.component_data(2)[c]=0.0;
            }
            VelocityBoundaryConditions ubc{
                {"inlet",{VelocityBoundaryCondition::Type::FIXED_VALUE,{flow,0,0}}},
                {"outlet",{VelocityBoundaryCondition::Type::FIXED_VALUE,{flow,0,0}}},
                {"walls",{VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}}}};
            // Gauge sentinel excludes reference-row pinning from this isolated
            // matrix-assembly contract; production retains the real gauge.
            const auto legacy=assemble_pcd_pressure_operators(
                mesh,geometry,U,1.0,0.01,ubc,pbc,mesh.n_cells(),ConvectionScheme::UPWIND);
            const auto explicit_off=assemble_pcd_pressure_operators(
                mesh,geometry,U,1.0,0.01,ubc,pbc,mesh.n_cells(),ConvectionScheme::UPWIND,false);
            const auto aux=assemble_pcd_pressure_operators(
                mesh,geometry,U,1.0,0.01,ubc,pbc,mesh.n_cells(),ConvectionScheme::UPWIND,true);
            for (std::size_t row=0;row<mesh.n_cells();++row) {
                for (std::size_t k=legacy.convection_diffusion.row_offsets_data()[row];
                     k<legacy.convection_diffusion.row_offsets_data()[row+1];++k) {
                    EXPECT_NEAR(legacy.convection_diffusion.values_data()[k],
                                explicit_off.convection_diffusion.values_data()[k],0.0);
                }
            }
            const std::size_t inflow=flow>0 ? 0 : 1;
            const std::size_t outflow=1-inflow;
            EXPECT_TRUE(diagonal(aux.convection_diffusion,inflow) >
                        diagonal(legacy.convection_diffusion,inflow));
            EXPECT_NEAR(diagonal(aux.convection_diffusion,outflow),
                        diagonal(legacy.convection_diffusion,outflow),1e-12);
            EXPECT_TRUE(diagonal(aux.laplacian,inflow) > diagonal(legacy.laplacian,inflow));
        }
    });
    return run_all();
}
