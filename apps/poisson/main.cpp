#include "cfdx/core/solvers/scalar_diffusion.h"
#include "cfdx/io/gmsh/gmsh_importer.h"
#include "cfdx/core/memory/memory_planner.h"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::io::gmsh;

// ============================================
// Vertical Slice 1 — Poisson Solver Driver
// ============================================
// End-to-end demonstration of the canonical Phase 8 Dirichlet Poisson path:
//   Gmsh mesh import → memory planning → canonical PDE contract →
//   finite-volume assembly → CG solve.
//
// The source is a volumetric density S and the assembly integrates S_c V_c.
// Dirichlet values are supplied explicitly on boundary faces.

int main() {
    std::printf("=== CFDX Vertical Slice 1 — Poisson Solver ===\n");

    Mesh mesh;
    const bool mesh_ok =
        import_gmsh_mesh("tests/data/cavity_with_patches.msh", mesh);
    std::printf("Mesh import: %s (points=%zu, cells=%zu)\n",
                mesh_ok ? "PASS" : "FAIL",
                mesh.n_points(), mesh.numCells());
    if (!mesh_ok || mesh.numCells() == 0) {
        return 1;
    }

    std::vector<memory::BufferDescriptor> buffers;
    memory::MemoryPlanner planner;
    const auto budget_plan = planner.plan(
        buffers, 10, mesh.numCells(), mesh.numFaces(), 1, 5);
    std::printf("Budget: %.2f bytes/cell, feasible=%s\n",
                budget_plan.budget.bytes_per_cell,
                budget_plan.feasible ? "YES" : "NO");

    PoissonBoundaryCondition boundary =
        PoissonBoundaryCondition::dirichlet(mesh.n_faces());
    // Homogeneous Dirichlet data on every physical boundary face.
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) < 0) {
            boundary.face_values[f] = 0.0;
        }
    }

    const std::vector<double> source(mesh.numCells(), 1.0);
    const ScalarDiffusionConfig config{
        .diffusivity = 1.0,
        .max_iterations = 2000,
        .tolerance = 1e-10};

    const ScalarDiffusionResult result =
        solve_poisson_dirichlet(mesh, boundary, source, config);

    std::printf("Poisson system: %zu x %zu, nnz=%zu\n",
                result.matrix.numRows(),
                result.matrix.numCols(),
                result.matrix.nnz());
    std::printf("CG status=%s, iterations=%zu, residual=%.6e\n",
                to_string(result.linear_result.status),
                result.linear_result.iterations,
                result.linear_result.residual);

    return result.linear_result.converged() ? 0 : 1;
}
