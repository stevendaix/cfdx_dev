#pragma once
// Mock meshio header for builds without meshio
// This allows CFDX to compile without the external meshio dependency
// In production, this would be provided by the actual meshio library

namespace meshio {
    // Stub functions - real implementation would use meshio::read_mesh
    inline void initialize() {}
}
