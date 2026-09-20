// M0.9-T01 — HDF5 writer

#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include <string>

namespace cfdx {
namespace io {

bool write_mesh_hdf5(const std::string& filename, const cfdx::core::Mesh& mesh);
bool write_field_hdf5(const std::string& filename,
                      const cfdx::core::Field<double, cfdx::core::Location::CELL>& field);

}  // namespace io
}  // namespace cfdx