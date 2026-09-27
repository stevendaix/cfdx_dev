// M0.10-T04: CFDX case-level HDF5 I/O
//
// Spécification CFDX v0.7 §21/22/23 (extension issue #425)
// ---------------------------------------------------------------------------
// Reads and writes the full ``case.cfdx.h5`` interchange format that bridges
// the Python conversion layer and the C++ numerical core.
//
// File layout (all at root level):
//   Datasets:
//     points            — (n_points x 3) float64
//     face_vertices     — (n_face_verts) uint64
//     face_offsets      — (n_faces+1) uint64  CSR offsets
//     owner             — (n_faces) uint64
//     neighbour         — (n_faces) int64
//     cell_faces        — (n_cell_face_refs) uint64
//     cell_offsets      — (n_cells+1) uint64  CSR offsets
//   Groups:
//     /fields/scalar/<name>   — 1-D float64 per cell
//     /fields/vector/<name>   — (n_cells x dim) float64, SoA layout
//   Attributes (root):
//     format_version, schema_version, cfdx_version   — version metadata
//     topology_hash, geometry_hash, mesh_hash         — FNV-1a integrity
//     creation_date, modification_date, dimension,
//     precision, endian
//     source_solver, source_format, source_version,
//     source_case_path, source_case_name            — SourceInfo
//     case_setup_json                                — CaseSetup as JSON
//     gap_report_json                                — GapAnalysis as JSON
//     mesh_topology                                  — "cfdx-csr-v1"
//     boundary_patches (optional)                    — semicolon-separated
#pragma once

#include "cfdx/io/cfdx_io/io_interface.h"
#include "cfdx/io/cfdx_io/case_schema.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include <string>
#include <vector>
#include <utility>

namespace cfdx {
namespace io {

// ---------------------------------------------------------------------------
// Case-level read
// ---------------------------------------------------------------------------
// Reads the full case file: mesh topology, source metadata, case setup,
// and gap-analysis report.  Integrity hashes are validated *only* when
// present in the file (Python-written files omit them).
//
// Returns true on success.  On failure, ``mesh`` is left empty and the
// error is logged to stderr.
bool read_case_cfdx_h5(const std::string& filename,
                       cfdx::core::Mesh& mesh,
                       SourceInfo& source,
                       CaseSetup& setup,
                       GapAnalysis& gap);

// ---------------------------------------------------------------------------
// Case-level write
// ---------------------------------------------------------------------------
// Writes the full case file: mesh topology (with integrity hashes), source
// metadata, case setup JSON, and gap-analysis JSON.
bool write_case_cfdx_h5(const std::string& filename,
                        const cfdx::core::Mesh& mesh,
                        const SourceInfo& source,
                        const CaseSetup& setup,
                        const GapAnalysis& gap);

// ---------------------------------------------------------------------------
// Field read/write helpers (operate on an existing case file)
// ---------------------------------------------------------------------------
// Reads all scalar fields from /fields/scalar/<name>.
bool read_scalar_fields_hdf5(const std::string& filename,
                             std::vector<std::pair<std::string, cfdx::core::ScalarCellField>>& fields);

// Reads all vector fields from /fields/vector/<name>.
bool read_vector_fields_hdf5(const std::string& filename,
                             std::vector<std::pair<std::string, cfdx::core::Vec3CellField>>& fields);

// Writes all scalar and vector fields from a ConversionResult to /fields.
bool write_fields_hdf5(const std::string& filename,
                       const std::vector<cfdx::core::ScalarCellField>& scalars,
                       const std::vector<cfdx::core::Vec3CellField>& vectors);

}  // namespace io
}  // namespace cfdx
