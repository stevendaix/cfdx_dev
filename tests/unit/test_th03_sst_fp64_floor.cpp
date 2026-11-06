#include "cfdx/physics/sst_solver.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
int main() {
    using namespace cfdx::physics;
    using namespace cfdx::core;
    ScalarEquation eq;
    eq.matrix = SparseMatrix(2, 2);
    eq.matrix.push_back(0, 0, 1.e-7);
    eq.matrix.push_back(1, 1, 1.);
    eq.matrix.finalize();
    eq.rhs = Vector(2, 0.);
    eq.rhs(0) = 1.;
    eq.rhs(1) = 1.e5;
    const double target = 1.e-8;
    const double floor = 64. * std::numeric_limits<double>::epsilon();
    const double inner = sst_linear_tolerance_for_solution_defect(eq, target);
    if (!std::isfinite(inner) || std::abs(inner - floor) > 1.e-29) return 1;
    Vector exact(2, 0.);
    exact(0) = 1.e7;
    exact(1) = 1.e5;
    if (!(sst_equation_solution_defect(eq, exact) < target)) return 2;
    exact(1) += 1.e-2;
    if (!(sst_equation_solution_defect(eq, exact) > target)) return 3;
    std::cout << "PASS TH03 SST precision floor; independent defect gate retained\n";
    return 0;
}
