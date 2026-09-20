#include "convection_assembly.h"
#include <cstdio>

using namespace cfdx::core::numerics;

bool assembleConvectionCSR(const Mesh&, const std::vector<ScalarCellField>&,
                           const std::string&, int, SparseMatrix&) {
    std::printf("[ConvectionAssembly] Stub — requires velocity interpolation\n");
    return false;  // Stub — real implementation needs interpolation + flux calculation
}

bool assembleMomentumCSR(const Mesh&, const std::vector<ScalarCellField>&,
                         const ScalarCellField&, const std::string&, SparseMatrix&, Vector&) {
    std::printf("[Momentum] Stub — combines convection + diffusion + pressure gradient\n");
    return false;
}
