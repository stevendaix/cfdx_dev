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
    std::cout << "[3/5] Identifying Dirichlet cells...\n";
    std::vector<bool> is_dirichlet(mesh.n_cells(), false);

    // Méthode robuste : utiliser les PhysicalGroups du maillage (si disponibles)
    // Sinon : heuristique géométrique avec la taille caractéristique h
    const std::size_t n_pts = mesh.n_points();
    const std::size_t n_faces_now = mesh.n_faces();
    std::cout << "  Points: " << n_pts << ", Faces: " << n_faces_now << ", Cells: " << mesh.n_cells() << "\n";

    // Essayer d'utiliser les patches frontières depuis mesh.boundary()
    const auto& boundary_patches = mesh.boundary();
    bool has_boundary_patches = (boundary_patches.n_patches() > 0);
    std::cout << "  Boundary patches: " << (has_boundary_patches ? "YES" : "NO") << " (count=" << boundary_patches.n_patches() << ")\n";

    if (has_boundary_patches) {
        // Utiliser $PhysicalNames / $PhysicalGroups du maillage Gmsh
        // Les patches nommés "left", "right", "bottom", "top" sont Dirichlet (MMS)
        std::cout << "  PhysicalNames/PhysicalGroups détectés : " << boundary_patches.n_patches() << " patches\n";
        for (std::size_t p = 0; p < boundary_patches.n_patches(); ++p) {
            const auto& patch = boundary_patches.patch(p);
            const std::string& patch_name = patch.name;  // .name est un attribut public (string)
            std::cout << "  Patch: [" << p << "] name="" << patch_name << "" (faces=" << patch.size() << ")\n";
            // Marquer toutes les cellules adjacentes aux faces du patch comme Dirichlet
            for (auto f_id : patch.face_ids) {
                std::size_t f = static_cast<std::size_t>(f_id);
                if (f < mesh.n_faces()) {
                    std::size_t owner_cell = static_cast<std::size_t>(mesh.ownership().owner(f));
                    if (owner_cell < mesh.n_cells()) {
                        is_dirichlet[owner_cell] = true;
                    }
                }
            }
        }
    } else {
        std::cout << "  [WARN] Pas de PhysicalGroups/BoundaryPatches. Utilisation heuristique.\n";
        // Heuristique géométrique : centre de cellule proche du bord [0,1]²
        const auto& pts = mesh.points();
        const double* px = pts.x_data();
        const double* py = pts.y_data();
        const double* pz = pts.z_data();
        const auto* c_faces_data = mesh.cells().faces_data();
        const auto* c_offsets_data = mesh.cells().offsets_data();

        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            const auto off = c_offsets_data[c];
            const auto n = c_offsets_data[c + 1] - off;
            double cx = 0.0, cy = 0.0, count = 0.0;
            for (auto k = off; k < c_offsets_data[c + 1]; ++k) {
                std::size_t f_idx = c_faces_data[k];
                const auto f_off = mesh.faces().offsets_data()[f_idx];
                const auto f_n = mesh.faces().offsets_data()[f_idx + 1] - f_off;
                const auto* verts = mesh.faces().vertices_data();
                for (std::size_t v = 0; v < f_n; ++v) {
                    std::size_t node_idx = verts[f_off + v];
                    cx += px[node_idx]; cy += py[node_idx]; count += 1.0;
                }
            }
            if (count > 0) {
                cx /= count; cy /= count;
                double h = 1.0 / std::sqrt(static_cast<double>(mesh.n_cells()));
                const double eps = 0.5 * h;
                if (cx < eps || cx > 1.0 - eps || cy < eps || cy > 1.0 - eps) {
                    is_dirichlet[c] = true;
                }
            }
        }
    }

    std::size_t n_dirichlet = std::count(is_dirichlet.begin(), is_dirichlet.end(), true);
    std::cout << "  Dirichlet cells: " << n_dirichlet << " / " << mesh.n_cells()
              << " (" << (100.0 * n_dirichlet / mesh.n_cells()) << "%)\n\n";
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
