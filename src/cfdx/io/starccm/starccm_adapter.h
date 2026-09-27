// M0.10-T03: STAR-CCM+ .sim adapter — stub
//
// Format parsing has moved to the Python conversion layer.
// Architecture: "Python knows the formats, C++ only knows CFDX canonical HDF5."
// This C++ stub only identifies the source by extension; all conversion
// delegates to the Python converter.
#pragma once

#include "cfdx/io/cfdx_io/io_interface.h"
#include <string>

namespace cfdx {
namespace io {
namespace starccm {

class StarCCMAdapter : public SolverAdapter {
public:
    StarCCMAdapter() = default;

    bool detect_source(const std::string& case_path,
                       SourceInfo& info) override;

    bool convert(const std::string& case_path,
                 ConversionResult& result) override;
};

}  // namespace starccm
}  // namespace io
}  // namespace cfdx
