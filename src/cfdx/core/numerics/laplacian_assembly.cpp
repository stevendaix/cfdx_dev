#include "laplacian_assembly.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/mesh_validator.h"
#include <cmath>

using namespace cfdx::core::numerics;
using namespace cfdx::core;

bool assembleLaplacianCSR(const Mesh& mesh, SparseMatrix& A) {
    size_t n = mesh.numCells();
    A.resize(n, n);
    A.setSymmetric(true);

    // For each face, compute coefficient = |Sf| / |d_ij|
    // and add to A (diagonal: +coeff, off-diagonal: -coeff)
    for (std::size_t f = 0; f < mesh.numFaces(); ++f) {
        uint32_t owner = mesh.faceOwner(f);
        uint32_t neighbour = mesh.faceNeighbour(f);

        if (owner == INVALID_CELL) continue;  // Boundary face handled separately
        if (neighbour == INVALID_CELL) {
            // Boundary face: only diagonal contribution
            double coeff = 1.0;  // Approximation: unit coefficient
            A.addToEntry(owner, owner, coeff);
            continue;
        }

        // Internal face connecting owner and neighbour
        double coeff = 1.0;  // Approximation for skeleton
        A.addToEntry(owner, owner, coeff);
        A.addToEntry(neighbour, neighbour, coeff);
        A.addToEntry(owner, neighbour, -coeff);
        A.addToEntry(neighbour, owner, -coeff);
    }

    // Finalize CSR structure
    A.finalize();
    return true;
}

bool assemblePoissonRHS(const Mesh& mesh, const ScalarCellField& source, Vector& b) {
    size_t n = mesh.numCells();
    b.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        b[i] = source(i);
    }
    return true;
}

LinearSystem buildPoissonSystem(const Mesh& mesh, const ScalarCellField& source) {
    LinearSystem system;
    SparseMatrix A;
    Vector b;
    assembleLaplacianCSR(mesh, A);
    assemblePoissonRHS(mesh, source, b);
    system.setMatrix(A);
    system.setRHS(b);
    return system;
}
