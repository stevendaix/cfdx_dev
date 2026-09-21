#include "gmres_solver.h"
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

namespace cfdx::core::linalg {

bool gmres_solve(const SparseMatrix& A,
                 const std::vector<double>& b,
                 std::vector<double>& x,
                 int restart,
                 double tolerance,
                 int max_iter) {
    int n = static_cast<int>(b.size());
    x.assign(n, 0.0);

    std::vector<double> r = b;
    std::vector<std::vector<double>> V(restart + 1, std::vector<double>(n, 0.0));
    std::vector<double> h(restart + 1, 0.0);

    double beta = std::sqrt(std::inner_product(r.begin(), r.end(), r.begin(), 0.0));
    if (beta == 0.0) return true;

    V[0] = r;
    for (int i = 0; i < n; ++i) V[0][i] /= beta;

    for (int iter = 1; iter <= max_iter; ++iter) {
        std::vector<double> w = A.matvec(V[iter % restart]);
        for (int i = 0; i < iter && i < restart; ++i) {
            double hi = std::inner_product(V[i].begin(), V[i].end(), w.begin(), 0.0);
            h[i] = hi;
            for (int j = 0; j < n; ++j) w[j] -= hi * V[i][j];
        }
        double h_iter = std::sqrt(std::inner_product(w.begin(), w.end(), w.begin(), 0.0));
        h[iter] = h_iter;
        for (int j = 0; j < n; ++j) V[iter][j] = w[j] / h_iter;

        double res_norm = std::abs(h[static_cast<std::size_t>(std::min(iter, restart) - 1)]);
        if (res_norm < tolerance) return true;
    }
    return false;
}
} // namespace
