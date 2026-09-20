// M0.9-T02 — HDF5 reader

#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include <string>

namespace cfdx {
namespace io {

bool read_mesh_hdf5(const std::string& filename, cfdx::core::Mesh& mesh);
bool read_field_hdf5(const std::string& filename,
                     cfdx::core::Field<double, cfdx::core::Location::CELL>& field);

}  // namespace io
}  // namespace cfdx