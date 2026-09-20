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
    std::cout << "=== Poisson MMS Validation ===\n";
    std::cout << "  Input : " << mesh_path << "\n";
    std::cout << "  Output: " << output_path << "\n\n";
    std::cout << "[1/5] Import Gmsh...\n";
    cfdx::core::Mesh mesh;
    if (!cfdx::io::import_gmsh_mesh(mesh_path, mesh)) {
        std::cerr << "Failed to load mesh: " << mesh_path << "\n";
        return 1;
    }
    std::cout << "  cells = " << mesh.n_cells() << "\n\n";
    std::cout << "[2/5] RCM & Memory Planner...\n";
    cfdx::core::RCMReorderer rcm;
    rcm.apply(mesh);
    cfdx::core::MemoryPlanner planner;
    auto budget = planner.estimate(mesh);
    std::cout << "  bytes/cell = " << budget.bytes_per_cell << "\n\n";
    std::cout << "[3/5] Identifying Dirichlet cells...\n";
    std::vector<bool> is_dirichlet(mesh.n_cells(), false);
    const auto& pts = mesh.points();
    const double* px = pts.x_data();
    const double* py = pts.y_data();
    const double* pz = pts.z_data();
    const auto* verts = mesh.faces().vertices_data();
    const auto* f_offsets = mesh.faces().offsets_data();
    const auto* c_faces = mesh.cells().faces_data();
    const auto* c_offsets = mesh.cells().offsets_data();
    std::vector<double> cx_arr(mesh.n_cells(), 0.0);
    std::vector<double> cy_arr(mesh.n_cells(), 0.0);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        const auto off = c_offsets[c];
        const auto n = c_offsets[c + 1] - off;
        double cx = 0.0, cy = 0.0, count = 0.0;
        for (auto k = off; k < c_offsets[c + 1]; ++k) {
            const auto f_idx = c_faces[k];
            const auto v_off = f_offsets[f_idx];
            const auto v_n = f_offsets[f_idx + 1] - v_off;
            for (auto vi = 0; vi < v_n; ++vi) {
                const auto v_idx = verts[v_off + vi];
                cx += px[v_idx]; cy += py[v_idx]; count += 1.0;
            }
        }
        cx_arr[c] = cx / std::max(count, 1.0);
        cy_arr[c] = cy / std::max(count, 1.0);
    }
    const double eps = 0.05;
    int dirichlet_count = 0;
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        bool is_edge = (cx_arr[c] < eps || cx_arr[c] > 1.0 - eps || cy_arr[c] < eps || cy_arr[c] > 1.0 - eps);
        is_dirichlet[c] = is_edge;
        if (is_edge) dirichlet_count++;
    }
    std::cout << "  Dirichlet cells = " << dirichlet_count << " / " << mesh.n_cells() << "\n\n";
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
    return 0;
}
