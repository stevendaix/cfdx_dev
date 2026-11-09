#include "cfdx/physics/pcd_fp_strong_component_anchor.h"
#include "cfdx/physics/pcd_fp_capture.h"
#include "cfdx/physics/pcd_fp_equilibration.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/advanced_preconditioners.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::core;
using namespace cfdx::physics;

static bool verify_synthetic() {
    SparseMatrix A(4, 4);
    A.push_back(0, 0, 2); A.push_back(0, 1, -1);
    A.push_back(1, 0, -1); A.push_back(1, 1, 2); A.push_back(1, 2, -1);
    A.push_back(2, 2, 2); A.push_back(2, 3, -2);
    A.push_back(3, 2, -1); A.push_back(3, 3, 1);
    A.finalize();
    auto gauge = anchor_pcd_fp_strong_components(A);
    if (gauge.anchored_rows != std::vector<std::size_t>{2}) return false;
    Vector b(4); b(0) = 1; b(1) = 2; b(2) = 3; b(3) = 4;
    const auto converted = gauge.rhs(b);
    if (converted(2) != 0 || converted(3) != 4) return false;
    if (A(0, 1) != -1 || A(2, 3) != -2) return false; // ORIGINAL unchanged
    if (gauge.matrix(2, 2) != 2 || gauge.matrix(2, 3) != 0 ||
        gauge.matrix(3, 2) != 0) return false;
    SparseMatrix B(2, 2);
    B.push_back(0, 0, 2); B.push_back(0, 1, -1);
    B.push_back(1, 0, -1); B.push_back(1, 1, 2);
    B.finalize();
    try { (void)anchor_pcd_fp_strong_components(B); return false; }
    catch (const std::runtime_error&) {}
    return true;
}

static bool replay_real_capture(const std::string& path) {
    const auto input = read_pcd_fp_capture(path);
    auto gauge = anchor_pcd_fp_strong_components(input.matrix);
    auto x = input.initial;
    gauge.pin_initial(x);
    const auto initial = x;
    auto rhs = gauge.rhs(input.rhs);
    ILU0Preconditioner pc;
    if (!pc.setup(gauge.matrix))
        throw std::runtime_error("PCD Fp anchored ILU0 setup rejected");
    const auto status = solve_fgmres(gauge.matrix, rhs, x,
        std::min<int>(512, static_cast<int>(x.size())), 100, 1e-8,
        &pc, {}, {}, false);
    const auto original = pcd_fp_true_residual(input.matrix, input.rhs, input.initial, x);
    const auto auxiliary = pcd_fp_true_residual(gauge.matrix, rhs, initial, x);
    std::cout << std::setprecision(16)
              << "PCD_FP_ANCHORED_REPLAY anchors=" << gauge.anchored_rows.size()
              << " rows=";
    for (const auto r : gauge.anchored_rows) std::cout << r << ",";
    std::cout << " status=" << static_cast<int>(status.status)
              << " iterations=" << status.iterations
              << " original_initial_l2=" << original.initial_l2
              << " original_true_relative_rhs=" << original.relative_rhs
              << " auxiliary_true_relative_rhs=" << auxiliary.relative_rhs
              << "\n";
    return status.status == SolverStatus::CONVERGED &&
        auxiliary.relative_rhs <= 1e-8 &&
        original.final_l2 < original.initial_l2 &&
        gauge.anchored_rows.size() >= 1 && x.is_valid();
}

int main(int argc, char** argv) {
    try {
        if (!verify_synthetic()) return 2;
        if (argc == 1) {
            std::cout << "PCD_FP_STRONG_COMPONENT_SYNTHETIC_PASS\n";
            return 0;
        }
        if (argc != 2) throw std::invalid_argument("usage: test_pcd_fp_strong_component_anchor [capture.cfp]");
        return replay_real_capture(argv[1]) ? 0 : 3;
    } catch (const std::exception& e) {
        std::cerr << "PCD_FP_STRONG_COMPONENT_TEST_ERROR: " << e.what() << "\n";
        return 4;
    }
}
