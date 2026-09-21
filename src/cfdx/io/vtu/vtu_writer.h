// M0.9-T05 — Lightweight VTU Writer (XML VTK Unstructured Grid)
//
// Export CFDX mesh + fields to .vtu for ParaView/PyVista visualization.
// No external VTK dependency — pure XML output.
//
// Format: VTK XML Unstructured Grid (.vtu)
// Supports: polyhedral cells (decomposed to tet/hex/pyramid), cell/face/point fields.

#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/index_types.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdint>

namespace cfdx {
namespace io {

// VTK cell types (subset needed for polyhedral decomposition)
enum class VtkCellType : std::uint8_t {
    VERTEX = 1,
    LINE = 3,
    TRIANGLE = 5,
    QUAD = 9,
    TETRA = 10,
    HEXAHEDRON = 12,
    WEDGE = 13,      // prism
    PYRAMID = 14,
    POLYHEDRON = 42  // VTK 9.0+ native polyhedron
};

// Convert CFDX face to VTK faces (triangulation for polyhedral cells)
struct VtuWriter {
    VtuWriter() = default;

    // Write mesh + fields to .vtu file
    // fields_cell: map of name -> ScalarCellField
    // fields_face: map of name -> ScalarFaceField (interpolated to cell centers for visualization)
    // fields_point: map of name -> ScalarPointField
    bool write(const std::string& filename,
               const cfdx::core::Mesh& mesh,
               const std::map<std::string, cfdx::core::ScalarCellField>& fields_cell = {},
               const std::map<std::string, cfdx::core::ScalarFaceField>& fields_face = {},
               const std::map<std::string, cfdx::core::ScalarPointField>& fields_point = {});

private:
    // Decompose polyhedral cells to VTK-supported types (tets)
    void decompose_polyhedra(const cfdx::core::Mesh& mesh,
                             std::vector<std::vector<cfdx::core::PointIndex>>& vtk_cells,
                             std::vector<VtkCellType>& vtk_cell_types,
                             std::vector<std::uint32_t>& cell_face_offsets,
                             std::vector<cfdx::core::FaceIndex>& cell_face_indices);

    // Write XML header + points
    void write_header(std::ofstream& os, const cfdx::core::Mesh& mesh,
                      const std::vector<std::vector<cfdx::core::PointIndex>>& vtk_cells,
                      const std::vector<VtkCellType>& vtk_cell_types);

    // Write cell connectivity
    void write_cells(std::ofstream& os,
                     const std::vector<std::vector<cfdx::core::PointIndex>>& vtk_cells,
                     const std::vector<VtkCellType>& vtk_cell_types);

    // Write fields
    void write_cell_fields(std::ofstream& os,
                           const cfdx::core::Mesh& mesh,
                           const std::map<std::string, cfdx::core::ScalarCellField>& fields,
                           const std::vector<std::size_t>& vtk_to_original_cell);

    void write_point_fields(std::ofstream& os,
                            const cfdx::core::Mesh& mesh,
                            const std::map<std::string, cfdx::core::ScalarPointField>& fields);

    // Helper: escape XML special chars
    static std::string xml_escape(const std::string& s);
};

}  // namespace io
}  // namespace cfdx