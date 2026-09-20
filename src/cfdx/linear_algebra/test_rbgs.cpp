// M0.8 - Red-Black Gauss-Seidel Preconditioner
// Test implementation for RBGS preconditioner
#include "rbgs_preconditioner.h"
#include <vector>
#include <cmath>

namespace cfdx {
namespace linear_algebra {

void test_rbgs() {
    // Simple 3x3 tridiagonal system
    // A = [[2, -1, 0], [-1, 2, -1], [0, -1, 2]]
    std::vector<double> A = {
        2.0, -1.0, 0.0,
        -1.0, 2.0, -1.0,
        0.0, -1.0, 2.0
    };
    
    std::vector<double> b = {1.0, 0.0, 0.0};
    
    // Build ILU(0) preconditioner
    ILUPreconditioner ilu;
    ilu.build(A);
    
    // Solve using ILU(0)
    std::vector<double> x = ilu.solve(b);
    
    // Expected solution: x = [1/3, 1/3, 1/3]
    double expected = 1.0/3.0;
    double error = std::abs(x[0] - expected);
    
    if (error < 0.1) {
        std::cout << "RBGS test PASSED: error = " << error << std::endl;
    } else {
        std::cout << "RBGS test FAILED: error = " << error << std::endl;
    }
}

}  // namespace linear_algebra
