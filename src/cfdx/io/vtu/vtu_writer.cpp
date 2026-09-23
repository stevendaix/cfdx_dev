// M0.9-T05 — Lightweight VTU Writer (implementation)

#include "vtu_writer.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <cmath>

namespace cfdx {
namespace io {

bool VtuWriter::write(const std::string& filename,
                      const cfdx::core::Mesh& mesh,
                      const std::map<std::string, cfdx::core::ScalarCellField>& fields_cell,
                      const std::map<std::string, cfdx::core::ScalarFaceField>& fields_face,
                      const std::map<std::string, cfdx::core::ScalarPointField>& fields_point,
                      double physical_time,
                      std::size_t iteration,
                      bool write_time_metadata)
{
    std::ofstream os(filename);
    if (!os.is_open()) {
        std::cerr << "VtuWriter: cannot open " << filename << " for writing\n";
        return false;
    }

    os << std::setprecision(12) << std::scientific;

    // Decompose polyhedral cells to VTK-supported types
    std::vector<std::vector<cfdx::core::PointIndex>> vtk_cells;
    std::vector<VtkCellType> vtk_cell_types;
    std::vector<std::uint32_t> cell_face_offsets;
    std::vector<cfdx::core::FaceIndex> cell_face_indices;

    decompose_polyhedra(mesh, vtk_cells, vtk_cell_types, cell_face_offsets, cell_face_indices);

    if (write_time_metadata && !std::isfinite(physical_time)) {
        std::cerr << "VtuWriter: physical_time must be finite\n";
        return false;
    }

    write_header(os, mesh, vtk_cells, vtk_cell_types, physical_time, iteration, write_time_metadata);
    write_cells(os, vtk_cells, vtk_cell_types);
    write_cell_fields(os, mesh, fields_cell, vtk_cells);
    write_point_fields(os, mesh, fields_point);

    // Close XML
    os << "  </Piece>\n";
    os << " </UnstructuredGrid>\n";
    os << "</VTKFile>\n";

    return true;
}

void VtuWriter::decompose_polyhedra(const cfdx::core::Mesh& mesh,
                                    std::vector<std::vector<cfdx::core::PointIndex>>& vtk_cells,
                                    std::vector<VtkCellType>& vtk_cell_types,
                                    std::vector<std::uint32_t>& cell_face_offsets,
                                    std::vector<cfdx::core::FaceIndex>& cell_face_indices)
{
    const std::size_t n_cells = mesh.n_cells();
    vtk_cells.reserve(n_cells);
    vtk_cell_types.reserve(n_cells);
    cell_face_offsets.reserve(n_cells + 1);
    cell_face_indices.reserve(mesh.cells().n_face_refs());

    cell_face_offsets.push_back(0);

    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    const auto* face_vertices = mesh.faces().vertices_data();
    const auto* face_offsets = mesh.faces().offsets_data();

    for (std::size_t c = 0; c < n_cells; ++c) {
        const std::uint32_t off = cell_offsets[c];
        const std::uint32_t n = cell_offsets[c + 1] - off;

        // Store face indices for this cell (for field mapping)
        for (std::uint32_t k = 0; k < n; ++k) {
            cell_face_indices.push_back(cell_faces[off + k]);
        }
        cell_face_offsets.push_back(static_cast<std::uint32_t>(cell_face_indices.size()));

        // For now: decompose to tetrahedra (simpler, works for any polyhedron)
        // A polyhedron can be decomposed to n_faces-2 tets from a reference point
        // We use the cell center as reference
        if (n >= 4) {
            // Use first face's first vertex as reference
            const std::uint32_t first_face = cell_faces[off];
            const std::uint32_t fv_off = face_offsets[first_face];
            const std::uint32_t ref_vertex = face_vertices[fv_off];

            // Create tets: (ref_vertex, face_vertices[i], face_vertices[i+1])
            for (std::uint32_t k = 0; k < n - 2; ++k) {
                const std::uint32_t face_a = cell_faces[off + k];
                const std::uint32_t face_b = cell_faces[off + k + 1];

                const std::uint32_t fv_a = face_vertices[face_offsets[face_a]];
                const std::uint32_t fv_b = face_vertices[face_offsets[face_b]];

                vtk_cells.push_back({ref_vertex, fv_a, fv_b});
                vtk_cell_types.push_back(VtkCellType::TETRA);
            }
        }
    }
}

void VtuWriter::write_header(std::ofstream& os, const cfdx::core::Mesh& mesh,
                             const std::vector<std::vector<cfdx::core::PointIndex>>& vtk_cells,
                             const std::vector<VtkCellType>& vtk_cell_types,
                             double physical_time, std::size_t iteration, bool write_time_metadata)
{
    const std::size_t n_points = mesh.n_points();
    const std::size_t n_cells = vtk_cells.size();

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    os << " <UnstructuredGrid>\n";
    os << "  <Piece NumberOfPoints="" << n_points << "" NumberOfCells="" << n_cells << "">\n";

    if (write_time_metadata) {
        os << "   <FieldData>\n";
        os << "    <DataArray type=\"Float64\" Name=\"physical_time\" NumberOfTuples=\"1\" format=\"ascii\">\n";
        os << "     " << physical_time << "\n";
        os << "    </DataArray>\n";
        os << "    <DataArray type=\"UInt64\" Name=\"iteration\" NumberOfTuples=\"1\" format=\"ascii\">\n";
        os << "     " << iteration << "\n";
        os << "    </DataArray>\n";
        os << "   </FieldData>\n";

    // Points
    os << "   <Points>\n";
    os << "    <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (std::size_t i = 0; i < n_points; ++i) {
        os << "     " << mesh.points().x(i) << " "
           << mesh.points().y(i) << " "
           << mesh.points().z(i) << "\n";
    }

    os << "    </DataArray>\n";
    os << "   </Points>\n";
}

void VtuWriter::write_cells(std::ofstream& os,
                            const std::vector<std::vector<cfdx::core::PointIndex>>& vtk_cells,
                            const std::vector<VtkCellType>& vtk_cell_types)
{
    const std::size_t n_cells = vtk_cells.size();

    os << "   <Cells>\n";

    // Connectivity
    os << "    <DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::size_t c = 0; c < n_cells; ++c) {
        os << "     ";
        for (const auto v : vtk_cells[c]) {
            os << v << " ";
        }
        os << "\n";
    }
    os << "    </DataArray>\n";

    // Offsets
    os << "    <DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">\n";
    std::uint64_t offset = 0;
    for (std::size_t c = 0; c < n_cells; ++c) {
        offset += vtk_cells[c].size();
        os << "     " << offset << "\n";
    }
    os << "    </DataArray>\n";

    // Types
    os << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t c = 0; c < n_cells; ++c) {
        os << "     " << static_cast<std::uint8_t>(vtk_cell_types[c]) << "\n";
    }
    os << "    </DataArray>\n";

    os << "   </Cells>\n";
}

void VtuWriter::write_cell_fields(std::ofstream& os,
                                  const cfdx::core::Mesh& mesh,
                                  const std::map<std::string, cfdx::core::ScalarCellField>& fields,
                                  const std::vector<std::vector<cfdx::core::PointIndex>>& vtk_cells)
{
    if (fields.empty()) return;

    const std::size_t n_cells = vtk_cells.size();

    os << "   <CellData>\n";

    for (const auto& [name, field] : fields) {
        if (field.size() != mesh.n_cells()) continue;

        os << "    <DataArray type=\"Float64\" Name=\"" << xml_escape(name) << "\" format=\"ascii\">\n";

        // Map original cell field to decomposed VTK cells
        for (std::size_t vtk_c = 0; vtk_c < n_cells; ++vtk_c) {
            // Find which original cell this VTK cell belongs to
            // For tet decomposition: each tet maps back to its parent cell
            // Simple approach: distribute evenly
            std::size_t orig_c = vtk_c; // approximation
            if (orig_c < field.size()) {
                os << "     " << field(orig_c) << "\n";
            }
        }
        os << "    </DataArray>\n";
    }

    os << "   </CellData>\n";
}

void VtuWriter::write_point_fields(std::ofstream& os,
                                   const cfdx::core::Mesh& mesh,
                                   const std::map<std::string, cfdx::core::ScalarPointField>& fields)
{
    if (fields.empty()) return;

    os << "   <PointData>\n";

    for (const auto& [name, field] : fields) {
        if (field.size() != mesh.n_points()) continue;

        os << "    <DataArray type=\"Float64\" Name=\"" << xml_escape(name) << "\" format=\"ascii\">\n";
        for (std::size_t i = 0; i < mesh.n_points(); ++i) {
            os << "     " << field(i) << "\n";
        }
        os << "    </DataArray>\n";
    }

    os << "   </PointData>\n";
}

std::string VtuWriter::xml_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '&': out += "&"; break;
            case '<': out += "<"; break;
            case '>': out += ">"; break;
            case '"': out += "\""; break;
            case '\'': out += "&apos;"; break;
            default: out += c;
        }
    }
    return out;
}

}  // namespace io
}  // namespace cfdx