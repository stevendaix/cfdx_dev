#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/io/hdf5/hdf5_writer.h"

namespace cfdx::io {

bool read_mesh_hdf5(const std::string&, cfdx::core::Mesh&) {
    return false;
}

bool read_field_hdf5(
    const std::string&,
    cfdx::core::Field<double, cfdx::core::Location::CELL>&) {
    return false;
}

bool write_mesh_hdf5(const std::string&, const cfdx::core::Mesh&) {
    return false;
}

bool write_field_hdf5(
    const std::string&,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>&) {
    return false;
}

} // namespace cfdx::io
