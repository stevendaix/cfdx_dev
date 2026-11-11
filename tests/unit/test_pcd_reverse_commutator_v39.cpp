#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/pcd_schur.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::core;

namespace {
void require(bool condition, const std::string& msg) {
    if (!condition) throw std::runtime_error(msg);
}

SparseMatrix dense2(double a, double b, double c, double d) {
    SparseMatrix m(2, 2);
    m.push_back(0, 0, a); m.push_back(0, 1, b);
    m.push_back(1, 0, c); m.push_back(1, 1, d);
    m.finalize();
    return m;
}

SparseMatrix diag2(double a, double d) {
    SparseMatrix m(2, 2);
    m.push_back(0, 0, a); m.push_back(1, 1, d);
    m.finalize();
    return m;
}

Vector exact_inverse2(const SparseMatrix& m, const Vector& b) {
    require(b.size() == 2, "test inverse incorrect vector shape");
    const double a = m(0, 0), c = m(1, 0);
    const double d = m(1, 1), t = m(0, 1);
    const double det = a*d - t*c;
    require(std::abs(det) > 1e-12, "test inverse singular");
    Vector out(2, 0.0);
    out(0) = (d*b(0) - t*b(1)) / det;
    out(1) = (-c*b(0) + a*b(1)) / det;
    return out;
}

bool close(const Vector& a, const Vector& b, double tol = 2e-12) {
    return a.size() == b.size() &&
           std::abs(a(0)-b(0)) < tol && std::abs(a(1)-b(1)) < tol;
}
}

int main() {
    try {
        // All three operators differ and Kp/Fp/Mp do not commute.
        const SparseMatrix mp = diag2(2.0, 5.0);
        const SparseMatrix kp = dense2(4.0, 1.0, -0.5, 3.0);
        const SparseMatrix fp = dense2(2.0, 1.5, -0.7, 1.0);
        const SparseMatrix uu = diag2(3.0, 4.0);
        const SparseMatrix g = diag2(1.0, 1.0);
        const SparseMatrix d = diag2(-1.0, -1.0);
        const SparseMatrix c = diag2(0.1, 0.2);
        const BlockOperator blocks(uu, g, d, c);
        blocks.validate();
        const auto solve_k = [kp](const Vector& rhs, Vector& x) {
            x = exact_inverse2(kp, rhs);
            return true;
        };
        const auto solve_f = [fp](const Vector& rhs, Vector& x) {
            x = exact_inverse2(fp, rhs);
            return true;
        };
        using Mode = PcdSchurActionMode;
        PcdSchurApproximation reverse(mp, kp, fp, solve_k, solve_f,
            std::nullopt, std::nullopt,
            Mode::ReverseLumpedMassCommutatorExperimental);
        PcdSchurApproximation original(mp, kp, fp, solve_k, solve_f,
            std::nullopt, std::nullopt, Mode::LumpedMassCommutator);
        require(reverse.setup(blocks) && original.setup(blocks),
            "valid reverse PCD setup rejected");
        require(std::string(reverse.name()) ==
            "pcd_reverse_lumped_mass_v39_experimental", "mode mislabeled");

        Vector rhs(2, 0.0); rhs(0) = 3.0; rhs(1) = -2.0;
        Vector z(2, 111.0), zold(2, 0.0);
        require(reverse.apply(rhs, z), "reverse apply rejected");
        require(original.apply(rhs, zold), "original apply rejected");
        Vector minv(2, 0.0);
        minv(0) = rhs(0) / 2.0;
        minv(1) = rhs(1) / 5.0;
        Vector product(2, 0.0);
        fp.matvec_into(minv, product);
        const Vector expected = exact_inverse2(kp, product) * (-1.0);
        require(close(z, expected), "reverse order -Kp^-1 Fp Mp^-1 incorrect");
        require((z-zold).norm2() > 0.01,
            "reverse accidentally implemented as existing order");
        Vector scaled(2, 0.0);
        require(reverse.apply(rhs * 3.0, scaled) && close(scaled, z * 3.0),
            "reverse PCD action is not linear");

        const SparseMatrix bad_mass = dense2(2.0, 0.1, 0.1, 5.0);
        PcdSchurApproximation bad(bad_mass, kp, fp, solve_k, solve_f,
            std::nullopt, std::nullopt,
            Mode::ReverseLumpedMassCommutatorExperimental);
        require(!bad.setup(blocks), "non-lumped mass silently accepted");
        PcdSchurApproximation negative(diag2(2.0, -5.0), kp, fp,
            solve_k, solve_f, std::nullopt, std::nullopt,
            Mode::ReverseLumpedMassCommutatorExperimental);
        require(!negative.setup(blocks), "negative pressure mass accepted");

        const auto new_k = dense2(4.1, 1.0, -0.5, 3.0);
        require(!reverse.update_pressure_values(blocks, mp, new_k, fp),
            "Kp mutated without updating nested solver callback");
        require(reverse.apply(rhs, scaled) && close(scaled, z),
            "failed transactional update corrupted existing PCD action");
        require(reverse.update_pressure_values_and_solves(
                    blocks, mp, new_k, fp, solve_k, solve_f),
            "explicit PCD numeric update rejected");
        Vector unchanged(2, 77.0);
        require(!reverse.apply(rhs, unchanged),
            "stale Kp callback passed independent true-residual check");
        require(unchanged(0) == 77.0 && unchanged(1) == 77.0,
            "failed PCD action altered caller result");

        std::cout << "PCD_V39_REVERSE_PASS noncommuting=1 transactional=1 "
                     "stale_callback_rejected=1 physics_qualified=0\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "PCD_V39_REVERSE_FAIL " << error.what() << '\n';
        return 2;
    }
}
