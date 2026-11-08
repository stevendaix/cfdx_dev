// CPU read-only verification of a two-word pressure correction on the exact
// CFDX original CSR capture. This never certifies Rhie-Chow/FV face fluxes.
#include "cfdx/core/linalg/linear_system_capture.h"
#include "cfdx/core/linalg/compensated_pressure_balance.h"
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    using namespace cfdx::core;
    if (argc != 3) {
        std::cerr << "usage: cfdx_compensated_pressure_replay capture.cfdxlin correction.split\n";
        return 2;
    }
    try {
        auto capture = read_linear_system_capture(argv[1]);
        std::ifstream input(argv[2]);
        if (!input) throw std::runtime_error("cannot open split file");
        std::string magic, count, checksum, line;
        if (!std::getline(input, magic) || magic != "CFDX_PRESSURE_SPLIT_V1")
            throw std::runtime_error("invalid split version");
        if (!std::getline(input, line)) throw std::runtime_error("missing split header");
        std::istringstream head(line);
        std::size_t rows = 0;
        if (!(head >> rows >> checksum) || (head >> line) ||
            rows != capture.matrix.n_rows())
            throw std::runtime_error("split CSR dimensions mismatch");
        std::ostringstream expected;
        expected << std::hex << std::setfill('0') << std::setw(16) << capture.checksum;
        if (checksum != expected.str())
            throw std::runtime_error("split does not belong to this exact CSR capture");
        Vector high(rows, 0.0), low(rows, 0.0);
        for (std::size_t i = 0; i < rows; ++i) {
            if (!std::getline(input, line))
                throw std::runtime_error("truncated pressure split");
            std::istringstream values(line);
            std::string h, l, extra;
            if (!(values >> h >> l) || (values >> extra))
                throw std::runtime_error("invalid split row");
            std::size_t used_h = 0, used_l = 0;
            high(i) = std::stod(h, &used_h);
            low(i) = std::stod(l, &used_l);
            if (used_h != h.size() || used_l != l.size() ||
                !std::isfinite(high(i)) || !std::isfinite(low(i)))
                throw std::runtime_error("nonfinite or invalid split value");
        }
        if (std::getline(input, line))
            throw std::runtime_error("extra rows in split file");
        const auto p = audit_compensated_pressure_balance(
            capture.matrix, capture.rhs, high, low);
        std::cout << std::scientific << std::setprecision(17)
            << "CFDX_COMPENSATED_CSR rows=" << rows
            << " symmetric_m_matrix=" << p.symmetric_m_matrix
            << " original_relative=" << p.original_relative_l2
            << " algebraic_face_relative=" << p.algebraic_face_relative_l2
            << " max_row_face_identity=" << p.max_row_face_identity
            << " internal_pair_net=" << p.global_internal_pair_defect
            << " fv_face_certified=" << p.physical_face_flux_certified
            << " tolerance=" << capture.tolerance << '\n';
        const bool accepted = p.valid && p.symmetric_m_matrix &&
            p.original_relative_l2 <= capture.tolerance &&
            p.algebraic_face_relative_l2 <= capture.tolerance;
        if (!accepted) std::cerr << "NOT_ALGEBRAICALLY_CERTIFIED\n";
        return accepted ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "CFDX_COMPENSATED_CSR_FAIL: " << e.what() << '\n';
        return 1;
    }
}
