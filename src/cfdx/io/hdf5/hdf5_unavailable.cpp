#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/io/hdf5/hdf5_writer.h"
#include "cfdx/io/hdf5/case_hdf5_io.h"

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


// Linkable unsupported-feature stubs. Without native HDF5 the production
// executable must still link and reject .cfdx.h5 inputs cleanly at runtime.
// These functions never claim to read or write an unsupported HDF5 case.
bool read_case_cfdx_h5(const std::string&, cfdx::core::Mesh&,
                       SourceInfo&, CaseSetup&, GapAnalysis&) {
    return false;
}

bool write_case_cfdx_h5(const std::string&, const cfdx::core::Mesh&,
                        const SourceInfo&, const CaseSetup&, const GapAnalysis&) {
    return false;
}

bool read_scalar_fields_hdf5(const std::string&,
    std::vector<std::pair<std::string, cfdx::core::ScalarCellField>>&) {
    return false;
}

bool read_vector_fields_hdf5(const std::string&,
    std::vector<std::pair<std::string, cfdx::core::Vec3CellField>>&) {
    return false;
}

bool write_fields_hdf5(const std::string&,
    const std::vector<cfdx::core::ScalarCellField>&,
    const std::vector<cfdx::core::Vec3CellField>&) {
    return false;
}

} // namespace cfdx::io
