// M0.9-T05 — Lightweight VTU Writer (implementation)

#include "vtu_writer.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <stdexcept>

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

    std::vector<std::vector<cfdx::core::PointIndex>> vtk_cells;
    std::vector<VtkCellType> vtk_cell_types;
    std::vector<std::uint32_t> cell_face_offsets;
    std::vector<cfdx::core::FaceIndex> cell_face_indices;
    std::vector<std::size_t> vtk_original_cell_indices;

    decompose_polyhedra(mesh, vtk_cells, vtk_cell_types, cell_face_offsets, cell_face_indices, vtk_original_cell_indices);

    if (write_time_metadata && !std::isfinite(physical_time)) {
        std::cerr << "VtuWriter: physical_time must be finite\n";
        return false;
    }

    write_header(os, mesh, vtk_cells, vtk_cell_types, physical_time, iteration, write_time_metadata);
    write_cells(os, vtk_cells, vtk_cell_types, cell_face_offsets, cell_face_indices, mesh);
    write_cell_fields(os, mesh, fields_cell, vtk_cells, vtk_original_cell_indices);
    write_point_fields(os, mesh, fields_point);

    os << "  </Piece>\n";
    os << " </UnstructuredGrid>\n";
    os << "</VTKFile>\n";

    return true;
}

void VtuWriter::decompose_polyhedra(const cfdx::core::Mesh& mesh,
                                    std::vector<std::vector<cfdx::core::PointIndex>>& vtk_cells,
                                    std::vector<VtkCellType>& vtk_cell_types,
                                    std::vector<std::uint32_t>& cell_face_offsets,
                                    std::vector<cfdx::core::FaceIndex>& cell_face_indices,
                                    std::vector<std::size_t>& vtk_original_cell_indices)
{
    const std::size_t n_cells = mesh.n_cells();
    vtk_cells.reserve(n_cells);
    vtk_cell_types.reserve(n_cells);
    cell_face_offsets.reserve(n_cells + 1);
    cell_face_indices.reserve(mesh.cells().n_face_refs());
    vtk_original_cell_indices.reserve(n_cells);

    cell_face_offsets.push_back(0);

    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    const auto* face_vertices = mesh.faces().vertices_data();
    const auto* face_offsets = mesh.faces().offsets_data();

    for (std::size_t c = 0; c < n_cells; ++c) {
        const std::uint32_t off = cell_offsets[c];
        const std::uint32_t n_faces = cell_offsets[c + 1] - off;

        std::vector<cfdx::core::PointIndex> connectivity;
        connectivity.reserve(1 + n_faces * 5);
                for (std::uint32_t k = 0; k < n_faces; ++k) {
            const std::uint32_t face = cell_faces[off + k];
            const std::uint32_t fv_begin = face_offsets[face];
            const std::uint32_t fv_end = face_offsets[face + 1];
            const std::uint32_t n_vertices = fv_end - fv_begin;

            if (n_vertices < 3) {
                throw std::runtime_error("VTU export: cell face has fewer than three vertices");
            }

            for (std::uint32_t v = fv_begin; v < fv_end; ++v) {
                const auto point = face_vertices[v];
                if (std::find(connectivity.begin(), connectivity.end(), point) == connectivity.end()) {
                    connectivity.push_back(point);
                }
            }
            cell_face_indices.push_back(face);
        }

        cell_face_offsets.push_back(static_cast<std::uint32_t>(cell_face_indices.size()));
        vtk_cells.push_back(std::move(connectivity));
        vtk_cell_types.push_back(VtkCellType::POLYHEDRON);
        vtk_original_cell_indices.push_back(c);
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

    if (write_time_metadata) {
        os << "  <FieldData>\n";
        os << "   <DataArray type=\"Float64\" Name=\"physical_time\" NumberOfTuples=\"1\" format=\"ascii\">\n";
        os << "    " << physical_time << "\n";
        os << "   </DataArray>\n";
        os << "   <DataArray type=\"UInt64\" Name=\"iteration\" NumberOfTuples=\"1\" format=\"ascii\">\n";
        os << "    " << iteration << "\n";
        os << "   </DataArray>\n";
        os << "  </FieldData>\n";
    }

    os << "  <Piece NumberOfPoints=\"" << n_points << "\" NumberOfCells=\"" << n_cells << "\">\n";
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
                            const std::vector<VtkCellType>& vtk_cell_types,
                            const std::vector<std::uint32_t>& cell_face_offsets,
                            const std::vector<cfdx::core::FaceIndex>& cell_face_indices,
                            const cfdx::core::Mesh& mesh)
{
    const std::size_t n_cells = vtk_cells.size();
    if (vtk_cell_types.size() != n_cells || cell_face_offsets.size() != n_cells + 1) {
        throw std::runtime_error("VTU export: cell topology arrays are inconsistent");
    }

    os << "   <Cells>\n";
    os << "    <DataArray type=\"UInt64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& cell : vtk_cells) {
        os << "     ";
        for (const auto point : cell) os << point << " ";
        os << "\n";
    }
    os << "    </DataArray>\n";

    os << "    <DataArray type=\"UInt64\" Name=\"offsets\" format=\"ascii\">\n";
    std::uint64_t offset = 0;
    for (const auto& cell : vtk_cells) {
        offset += cell.size();
        os << "     " << offset << "\n";
    }
    os << "    </DataArray>\n";

    os << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (const auto type : vtk_cell_types) {
        os << "     " << static_cast<std::uint8_t>(type) << "\n";
    }
    os << "    </DataArray>\n";

    bool has_polyhedron = false;
    for (const auto type : vtk_cell_types) has_polyhedron |= type == VtkCellType::POLYHEDRON;
    if (has_polyhedron) {
        os << "    <DataArray type=\"UInt64\" Name=\"faces\" format=\"ascii\">\n";
        std::uint64_t face_offset = 0;
        for (std::size_t c = 0; c < n_cells; ++c) {
            const auto begin = cell_face_offsets[c];
            const auto end = cell_face_offsets[c + 1];
            os << "     " << (end - begin) << " ";
            face_offset += 1;
            for (std::uint32_t i = begin; i < end; ++i) {
                const auto face = cell_face_indices[i];
                const auto fv_begin = mesh.faces().offsets_data()[face];
                const auto fv_end = mesh.faces().offsets_data()[face + 1];
                os << (fv_end - fv_begin) << " ";
                face_offset += 1;
                for (std::uint32_t v = fv_begin; v < fv_end; ++v) {
                    os << mesh.faces().vertices_data()[v] << " ";
                    ++face_offset;
                }
            }
            os << "\n";
        }
        os << "    </DataArray>\n";
        os << "    <DataArray type=\"UInt64\" Name=\"faceoffsets\" format=\"ascii\">\n";
        face_offset = 0;
        for (std::size_t c = 0; c < n_cells; ++c) {
            const auto begin = cell_face_offsets[c];
            const auto end = cell_face_offsets[c + 1];
            face_offset += 1;
            for (std::uint32_t i = begin; i < end; ++i) {
                const auto face = cell_face_indices[i];
                const auto fv_begin = mesh.faces().offsets_data()[face];
                const auto fv_end = mesh.faces().offsets_data()[face + 1];
                face_offset += 1 + (fv_end - fv_begin);
            }
            os << "     " << face_offset << "\n";
        }
        os << "    </DataArray>\n";
    }

    os << "   </Cells>\n";
}

void VtuWriter::write_cell_fields(std::ofstream& os,
                                  const cfdx::core::Mesh& mesh,
                                  const std::map<std::string, cfdx::core::ScalarCellField>& fields,
                                  const std::vector<std::vector<cfdx::core::PointIndex>>& vtk_cells,
                                  const std::vector<std::size_t>& vtk_original_cell_indices)
{
    if (fields.empty()) return;

    const std::size_t n_cells = vtk_cells.size();
    if (vtk_original_cell_indices.size() != n_cells) {
        throw std::runtime_error("VTU export: cell provenance size does not match VTK cell count");
    }

    os << "   <CellData>\n";

    for (const auto& [name, field] : fields) {
        if (field.size() != mesh.n_cells()) continue;

        os << "    <DataArray type=\"Float64\" Name=\"" << xml_escape(name) << "\" format=\"ascii\">\n";

        for (std::size_t vtk_c = 0; vtk_c < n_cells; ++vtk_c) {
            const std::size_t orig_c = vtk_original_cell_indices[vtk_c];
            if (orig_c >= field.size()) {
                throw std::runtime_error("VTU export: cell provenance index is out of range");
            }
            os << "     " << field(orig_c) << "\n";
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
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            case '\'': out += "&apos;"; break;
            default: out += c;
        }
    }
    return out;
}

}  // namespace io
}  // namespace cfdx
