// M0.8 - GMRES(m) Solver with Restart
// Test implementation for GMRES(m) solver
#include "gmres_solver.h"
#include <vector>
#include <cmath>

namespace cfdx {
namespace linear_algebra {

void test_gmres() {
    // Create a simple 3x3 tridiagonal system
    // A = [[2, -1, 0], [-1, 2, -1], [0, -1, 2]] (discrete Laplacian)
    std::vector<double> A = {
        2.0, -1.0, 0.0,
        -1.0, 2.0, -1.0,
        0.0, -1.0, 2.0
    };
    
    std::vector<double> b = {1.0, 0.0, 0.0};
    
    // Solve using GMRES with restart 2
    GMRES gmrES(2, 1e-10);
    gmrES.solve(b);
    
    // Expected solution for A*x = b
    // x = [1/3, 1/3, 1/3] approximately
    double x = gmrES.x;
    
    // Check convergence
    double residual = gmrES.norm(b - gmrES.solve(b)); // This would need a proper solve method
    
    // Simple manual verification
    double expected = 1.0/3.0;
    double error = std::abs(x[0] - expected);
    
    if (error < 0.1) {
        std::cout << "GMRES test PASSED: error = " << error << std::endl;
    } else {
        std::cout << "GMRES test FAILED: error = " << error << std::endl;
    }
}

}  // namespace linear_algebra
