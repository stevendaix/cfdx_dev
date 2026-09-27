#include "cfdx/physics/boundary_constraint_fvm.h"
#include "common/test_harness.h"
#include <cmath>
#include <memory>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("face_resolved_dirichlet_uses_face_context", [] {
        Mesh mesh;
        mesh.faces().push_face({0, 1});
        mesh.faces().push_face({2, 3});
        Patch p; p.name = "inlet"; p.type = PatchType::INLET; p.face_ids = {0, 1};
        mesh.boundary().add_patch(p);
        std::vector<Vec3> centres{{0.0, 0.0, 0.0}, {0.0, 2.0, 0.0}};
        BoundaryConstraintMap constraints;
        constraints["inlet"].emplace_back("T", Dirichlet{
            std::make_shared<SpatialLinearValueProvider>(1.0, 0.0, 2.0, 0.0)});
        const auto resolved = resolve_scalar_boundary_constraints(
            mesh, centres, constraints, "T");
        EXPECT_TRUE(resolved.has(0));
        EXPECT_TRUE(resolved.has(1));
        EXPECT_TRUE(std::abs(resolved.conditions[0].value - 1.0) < 1e-14);
        EXPECT_TRUE(std::abs(resolved.conditions[1].value - 5.0) < 1e-14);
    });

    run_case("face_resolved_neumann_uses_time_context", [] {
        Mesh mesh;
        mesh.faces().push_face({0, 1});
        Patch p; p.name = "outlet"; p.type = PatchType::OUTLET; p.face_ids = {0};
        mesh.boundary().add_patch(p);
        std::vector<Vec3> centres{{1.0, 2.0, 3.0}};
        BoundaryConstraintMap constraints;
        constraints["outlet"].emplace_back("T", Neumann{
            std::make_shared<LinearTimeValueProvider>(2.0, 0.5)});
        const auto resolved = resolve_scalar_boundary_constraints(
            mesh, centres, constraints, "T", 4.0);
        EXPECT_TRUE(resolved.has(0));
        EXPECT_TRUE(resolved.conditions[0].type == ScalarBoundaryFaceCondition::Type::FIXED_GRADIENT);
        EXPECT_TRUE(std::abs(resolved.conditions[0].value - 4.0) < 1e-14);
    });

    run_case("missing_field_constraint_is_rejected", [] {
        Mesh mesh;
        mesh.faces().push_face({0, 1});
        Patch p; p.name = "wall"; p.type = PatchType::WALL; p.face_ids = {0};
        mesh.boundary().add_patch(p);
        std::vector<Vec3> centres{{0.0, 0.0, 0.0}};
        BoundaryConstraintMap constraints;
        constraints["wall"].emplace_back("U.x", Dirichlet{std::make_shared<ConstantValueProvider>(0.0)});
        EXPECT_THROW(resolve_scalar_boundary_constraints(mesh, centres, constraints, "p"), std::invalid_argument);
    });

    run_case("unsupported_robin_is_rejected_explicitly", [] {
        Mesh mesh;
        mesh.faces().push_face({0, 1});
        Patch p; p.name = "wall"; p.type = PatchType::WALL; p.face_ids = {0};
        mesh.boundary().add_patch(p);
        std::vector<Vec3> centres{{0.0, 0.0, 0.0}};
        BoundaryConstraintMap constraints;
        constraints["wall"].emplace_back("T", Robin{1.0, 1.0, std::make_shared<ConstantValueProvider>(0.0)});
        EXPECT_THROW(resolve_scalar_boundary_constraints(mesh, centres, constraints, "T"), std::invalid_argument);
    });

    run_case("flux_dependent_uses_inflow_provider_on_backflow", [] {
        Mesh mesh;
        mesh.faces().push_face({0, 1});
        Patch p; p.name = "outlet"; p.type = PatchType::OUTLET; p.face_ids = {0};
        mesh.boundary().add_patch(p);
        std::vector<Vec3> centres{{0.0, 0.0, 0.0}};
        BoundaryConstraintMap constraints;
        auto inflow = std::make_shared<ConstantValueProvider>(7.0);
        auto outflow = std::make_shared<ConstantValueProvider>(0.0);
        constraints["outlet"].emplace_back("T", FluxDependent{inflow, outflow});
        Field<double, Location::FACE> flux(mesh.n_faces());
        flux(0) = -2.0;
        const auto resolved = resolve_scalar_boundary_constraints(
            mesh, centres, constraints, "T", 0.0, &flux);
        EXPECT_TRUE(resolved.has(0));
        EXPECT_TRUE(resolved.conditions[0].type ==
                    ScalarBoundaryFaceCondition::Type::FIXED_VALUE);
        EXPECT_TRUE(std::abs(resolved.conditions[0].value - 7.0) < 1e-14);
    });

    return run_all();
}