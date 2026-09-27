#include "cfdx/core/linalg/mgr_preconditioner.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace cfdx::core;

static SparseMatrix make_system(double fine_diag = 4.0, double coarse_diag = 3.0) {
    SparseMatrix A(8, 8);
    for (std::size_t i = 0; i < 4; ++i) {
        A.push_back(i, i, fine_diag);
        A.push_back(i, 4 + i, -1.0);
        A.push_back(4 + i, i, -1.0);
        A.push_back(4 + i, 4 + i, coarse_diag);
    }
    A.finalize();
    return A;
}

int main() {
    const std::vector<std::size_t> fine{0, 1, 2, 3};
    const std::vector<std::size_t> coarse{4, 5, 6, 7};

    NativeMGRPreconditioner mgr(fine, coarse);
    SparseMatrix A = make_system();
    assert(mgr.setup(A));
    assert(mgr.setup_count() == 1);
    assert(mgr.fine_size() == 4);
    assert(mgr.coarse_size() == 4);

    Vector r(8, 0.0);
    for (std::size_t i = 0; i < 8; ++i) r(i) = 1.0 + static_cast<double>(i);
    Vector z(8, 0.0);
    assert(mgr.apply(r, z));
    assert(z.is_valid());
    assert(z.norm_inf() > 0.0);

    SparseMatrix A2 = make_system(4.5, 3.25);
    assert(mgr.update_values(A2));
    assert(mgr.numeric_updates() == 1);

    Vector z2(8, 0.0);
    assert(mgr.apply(r, z2));
    assert(z2.is_valid());
    assert(std::abs(z2.norm_inf() - z.norm_inf()) > 1e-14);

    SparseMatrix changed = A2;
    // SparseMatrix is intentionally immutable after finalize; a new
    // topology is therefore represented by a fresh assembled matrix.
    SparseMatrix changed_pattern(8, 8);
    for (std::size_t i = 0; i < 4; ++i) {
        changed_pattern.push_back(i, i, 4.5);
        changed_pattern.push_back(i, 4 + i, -1.0);
        changed_pattern.push_back(4 + i, i, -1.0);
        changed_pattern.push_back(4 + i, 4 + i, 3.25);
    }
    changed_pattern.push_back(0, 5, 0.125);
    changed_pattern.finalize();
    assert(!mgr.update_values(changed_pattern));
    assert(mgr.last_error().find("pattern changed") != std::string::npos);

    std::cout << "Native MGR reduction tests passed\n";
    return 0;
}
