#include "cfdx/core/linalg/linear_system_capture.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "common/test_harness.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>

using namespace cfdx::core;
using namespace cfdx::testing;
namespace fs = std::filesystem;

int main() {
    const auto path = fs::temp_directory_path() / "cfdx_capture_qualification.cfdxlin";
    auto make_matrix = [] {
        SparseMatrix A(3, 3);
        A.push_back(0, 0, 5.0);
        A.push_back(0, 1, -2.0);
        A.push_back(1, 0, -2.0);
        A.push_back(1, 1, 5.0);
        A.push_back(1, 2, -1.0);
        A.push_back(2, 1, -1.0);
        A.push_back(2, 2, 3.0);
        A.finalize();
        return A;
    };
    run_case("capture_iteration_list_requires_unique_positive_numbers", [] {
        const auto selection = parse_linear_capture_iterations("20,21,600");
        EXPECT_TRUE(selection.size() == 3 && selection[2] == 600);
        EXPECT_THROW((void)parse_linear_capture_iterations("20,20"), std::invalid_argument);
        EXPECT_THROW((void)parse_linear_capture_iterations("20,abc"), std::invalid_argument);
        EXPECT_THROW((void)parse_linear_capture_iterations("0"), std::invalid_argument);
        EXPECT_THROW((void)parse_linear_capture_iterations("20,"), std::invalid_argument);
        EXPECT_THROW((void)parse_linear_capture_iterations(""), std::invalid_argument);
    });
    run_case("capture_exact_fp64_and_replay_initial_guess", [&] {
        const auto A = make_matrix();
        Vector b(3), x(3);
        b(0) = 1.25; b(1) = -2.0; b(2) = 3.25;
        x(0) = 0.125; x(1) = -0.125; x(2) = 0.75;
        write_linear_system_capture(path.string(), A, b, x,
            LinearProblemKind::Diffusion, NullSpaceModel::None,
            LinearScalingModel::None, 130, 1e-12);
        auto c = read_linear_system_capture(path.string());
        EXPECT_TRUE(c.matrix.nnz() == A.nnz());
        EXPECT_TRUE(c.max_iterations == 130);
        EXPECT_NEAR(c.tolerance, 1e-12, 0.0);
        for (std::size_t i = 0; i < b.size(); ++i) {
            EXPECT_NEAR(c.rhs(i), b(i), 0.0);
            EXPECT_NEAR(c.initial_solution(i), x(i), 0.0);
        }
        for (std::size_t i = 0; i < A.nnz(); ++i) {
            EXPECT_TRUE(c.matrix.columns_data()[i] == A.columns_data()[i]);
            EXPECT_NEAR(c.matrix.values_data()[i], A.values_data()[i], 0.0);
        }
        LinearSolverRequest request;
        request.krylov = KrylovModel::CG;
        request.preconditioner = PreconditionerModel::Jacobi;
        const auto report = solve_linear_system(c.matrix, c.rhs, c.initial_solution,
            c.problem, request, c.max_iterations, c.tolerance);
        EXPECT_TRUE(report.result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(report.physical_residual_relative < c.tolerance);
    });
    run_case("generic_dispatch_exports_only_when_opted_in", [&] {
        const auto A = make_matrix();
        Vector b(3, 1.0), x(3, 0.0);
        fs::remove(path);
        LinearSolverRequest request;
        request.krylov = KrylovModel::CG;
        request.preconditioner = PreconditionerModel::Jacobi;
        (void)solve_linear_system(A, b, x, LinearProblemKind::Diffusion,
            request, 50, 1e-12);
        EXPECT_TRUE(!fs::exists(path));
        request.linear_capture_path = path.string();
        x.fill(0.0);
        (void)solve_linear_system(A, b, x, LinearProblemKind::Diffusion,
            request, 50, 1e-12);
        const auto copy = read_linear_system_capture(path.string());
        EXPECT_NEAR(copy.initial_solution.norm2(), 0.0, 0.0);
        EXPECT_NEAR(copy.rhs(0), 1.0, 0.0);
    });
    run_case("checksum_detects_changed_matrix_without_accepting_it", [&] {
        std::ifstream source(path, std::ios::binary);
        const std::string data((std::istreambuf_iterator<char>(source)),
                               std::istreambuf_iterator<char>());
        auto modified = data;
        const auto marker = modified.find("4014000000000000"); // bit pattern for 5.0
        EXPECT_TRUE(marker != std::string::npos);
        modified.replace(marker, 16, "4018000000000000"); // 6.0
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        out << modified;
        out.close();
        EXPECT_THROW((void)read_linear_system_capture(path.string()), std::runtime_error);
    });
    run_case("reject_truncated_capture_and_retain_previous_data", [&] {
        std::ofstream out(path, std::ios::trunc);
        out << "CFDX_LINEAR_REPLAY_V1\n0 0 0\n3 3 7 50\n";
        out.close();
        EXPECT_THROW((void)read_linear_system_capture(path.string()), std::runtime_error);
    });
    run_case("refuse_nonfinite_or_bad_dimensions_before_writing", [&] {
        const auto A = make_matrix();
        Vector b(3, 1.0), x(3, 0.0);
        b(1) = std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW(write_linear_system_capture(path.string(), A, b, x,
            LinearProblemKind::General, NullSpaceModel::None,
            LinearScalingModel::None, 25, 1e-10), std::invalid_argument);
        EXPECT_THROW(write_linear_system_capture(path.string(), A,
            Vector(2, 1.0), x, LinearProblemKind::General,
            NullSpaceModel::None, LinearScalingModel::None, 25, 1e-10),
            std::invalid_argument);
        b.fill(1.0);
        EXPECT_THROW(write_linear_system_capture(path.string(), A, b, x,
            static_cast<LinearProblemKind>(999), NullSpaceModel::None,
            LinearScalingModel::None, 25, 1e-10), std::invalid_argument);
    });
    run_case("constant_nullspace_metadata_preserved", [&] {
        SparseMatrix A(2, 2);
        A.push_back(0, 0, 1.0); A.push_back(0, 1, -1.0);
        A.push_back(1, 0, -1.0); A.push_back(1, 1, 1.0);
        A.finalize();
        Vector b(2), x(2, 0.0);
        b(0) = -1.0; b(1) = 1.0;
        write_linear_system_capture(path.string(), A, b, x,
            LinearProblemKind::PressurePoisson, NullSpaceModel::Constant,
            LinearScalingModel::None, 70, 1e-11);
        const auto c = read_linear_system_capture(path.string());
        EXPECT_TRUE(c.null_space == NullSpaceModel::Constant);
        EXPECT_TRUE(c.problem == LinearProblemKind::PressurePoisson);
    });
    fs::remove(path);
    return run_all();
}
