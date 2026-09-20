#pragma once

#include "cfdx/io/hdf5/hdf5_writer.h"
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/index_types.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/geometry/face_geometry.h"
#include <string>
#include <vector>
#include <cstdint>

// Meshio adapter: optional; if meshio/meshio.h unavailable, basic parser used instead
#ifndef HAS_MESHIO
#define HAS_MESHIO 0
#endif

namespace cfdx::io::openfoam {

using ScalarCellField = cfdx::core::ScalarCellField;
using Vec3CellField = cfdx::core::Vec3CellField;

struct PatchDef {
    std::string name;
    int type;
    int n_elements;
    std::vector<std::uint32_t> face_ids;
    PatchDef(const char* n, int t, int ne, std::vector<std::uint32_t> ids)
        : name(n), type(t), n_elements(ne), face_ids(std::move(ids)) {}
};

struct OpenFOAMFace {
    int n;
    std::vector<uint32_t> indices;
};

bool import_openfoam_case(const std::string& ofCasePath, cfdx::core::Mesh& mesh);

bool import_openfoam_field(const std::string& ofFieldName,
                           const std::string& ofCasePath,
                           ScalarCellField& cfdxField);

bool import_openfoam_field(const std::string& ofFieldName,
                           const std::string& ofCasePath,
                           Vec3CellField& cfdxField);

bool read_openfoam_mesh_meshio(const std::string& ofCasePath,
                               std::vector<cfdx::core::Vec3>& points,
                               std::vector<OpenFOAMFace>& faces,
                               std::vector<PatchDef>& patches);

} // namespace cfdx::io::openfoam