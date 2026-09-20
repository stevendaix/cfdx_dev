#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/numerics/laplacian_assembly.h"
#include "cfdx/core/memory/memory_planner.h"
#include "cfdx/io/gmsh/gmsh_importer.h"
#include "cfdx/io/vtu/vtu_writer.h"
#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include <chrono>
constexpr double PI = 3.14159265358979323846;
inline double p_exact(double x, double y) {
    return std::sin(PI * x) * std::sin(PI * y);
}
inline double source_exact(double x, double y) {
    return 2.0 * PI * PI * std::sin(PI * x) * std::sin(PI * y);
}
int main(int argc, char** argv) {
    std::string mesh_path = "/tmp/mms_meshes/mms_square_32x32.msh";
    std::string output_path = "/tmp/mms_results/result.vtu";
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--mesh" && i + 1 < argc) mesh_path = argv[++i];
        else if (arg == "--output" && i + 1 < argc) output_path = argv[++i];
    }
    std::cout << "[DEBUG] Entrée dans main()\n";
    std::cout.flush();
    try {
        std::cout << "[DEBUG] Avant constructeur Mesh\n";
        std::cout.flush();
        cfdx::core::Mesh mesh;
        std::cout << "[DEBUG] Après constructeur Mesh, avant import\n";
        std::cout.flush();
        if (!cfdx::io::gmsh::import_gmsh_mesh(mesh_path, mesh)) {
            std::cerr << "Failed to load mesh: " << mesh_path << "\n";
            return 1;
        }
        if (!cfdx::io::gmsh::import_gmsh_mesh(mesh_path, mesh)) {
            std::cerr << "Failed to load mesh: " << mesh_path << "\n";
            return 1;
        }
        std::cout << "[DEBUG] Après import Gmsh, n_cells = " << mesh.n_cells() << "\n";
        std::cout.flush();
        const std::size_t n_cells = mesh.n_cells();
    if (n_cells == 0) {
        std::cerr << "ERREUR CRITIQUE : Le maillage est vide (n_cells == 0). Le parser Gmsh est peut-être encore un stub.\n";
        return 1;
    }
    std::cout << "[DEBUG] Après import. n_cells = " << mesh.n_cells() << "\n";
    std::cout << "  cells = " << mesh.n_cells() << "\n\n";
    std::cout << "[2/5] RCM & Memory Planner...\n";
    std::cout << "[DEBUG] Avant RCM\n";
    std::vector<uint32_t> owner, neighbour;
    cfdx::core::memory::RCMReorderer rcm;
    rcm.reorder(owner, neighbour, mesh.n_cells());
    cfdx::core::memory::MemoryPlanner planner;
    std::vector<cfdx::core::memory::BufferDescriptor> buffers;
    auto budget_plan = planner.plan(buffers, 10, mesh.n_cells(), mesh.n_faces(), 1, 5);
    std::cout << "  bytes/cell = " << budget_plan.budget.bytes_per_cell << "\n\n";
    std::cout << "[3/5] Identifying Dirichlet cells... (COMMENTÉ POUR TEST)\n";
    std::vector<bool> is_dirichlet(mesh.n_cells(), true);  // Tout Dirichlet = matrice identité
    std::cout << "  Dirichlet cells = SKIP (test)\n\n";
    std::cout << "  Dirichlet cells = SKIP (test simplifié)\n\n";
    std::cout << "[DEBUG] Avant assemble_laplacian_csr\n";
    std::cout.flush();
    std::cout << "[4/5] Laplacian assembly...\n";
    cfdx::core::SparseMatrix A;
    cfdx::core::Vector b;
    cfdx::core::assemble_laplacian_csr(mesh, A, b, is_dirichlet);
    std::cout << "  nnz = " << A.nnz() << "\n";
    std::cout << "  size = " << mesh.n_cells() << " x " << mesh.n_cells() << "\n\n";
    std::cout << "[5/5] CG solve...\n";
    cfdx::core::Vector x(mesh.n_cells(), 0.0);
    cfdx::core::SolverResult result = cfdx::core::solve_cg(A, b, x, 1000, 1e-10);
    std::cout << "  status = " << cfdx::core::to_string(result.status) << "\n";
    std::cout << "  iterations = " << result.iterations << "\n";
    std::cout << "  residual = " << result.residual << "\n\n";
    std::cout << "=== MMS Error Metrics ===\n";
    int n_side = static_cast<int>(std::round(std::sqrt(static_cast<double>(mesh.n_cells()))));
    const double V_cell = 1.0 / mesh.n_cells();
    double max_err = 0.0;
    double l2_err_sq = 0.0;
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        int ix = c % n_side;
        int iy = c / n_side;
        double xc = (ix + 0.5) / n_side;
        double yc = (iy + 0.5) / n_side;
        double p_num = x(c);
        double err = std::abs(p_num - p_exact(xc, yc));
        max_err = std::max(max_err, err);
        l2_err_sq += err * err * V_cell;
    }
    std::cout << "  L_inf error = " << max_err << "\n";
    std::cout << "  L_2 error   = " << std::sqrt(l2_err_sq) << "\n";
    } catch (const std::exception& e) {
        std::cerr << "[EXCEPTION] " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "[EXCEPTION INCONNUE]\n";
        return 1;
    }
    return 0;
}
