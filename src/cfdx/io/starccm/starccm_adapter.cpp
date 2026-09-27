// M0.10-T03: STAR-CCM+ adapter — stub (delegates to Python)
//
// Format parsing has moved to the Python conversion layer.
// Architecture: "Python knows the formats, C++ only knows CFDX canonical HDF5."
#include "starccm_adapter.h"

#include <algorithm>
#include <cctype>
#include <string>

namespace cfdx {
namespace io {
namespace starccm {

static bool ends_with_ci(const std::string& s, const std::string& ext) {
    if (s.size() < ext.size()) return false;
    std::string tail = s.substr(s.size() - ext.size());
    std::transform(tail.begin(), tail.end(), tail.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return tail == ext;
}

// Cheap extension-only detection. The Python layer performs content validation.
bool StarCCMAdapter::detect_source(const std::string& case_path,
                                   SourceInfo& info) {
    if (!ends_with_ci(case_path, ".sim")) {
        return false;
    }
    info.solver = "STAR-CCM+";
    info.format = "sim";
    info.version = "unknown";
    info.case_path = case_path;
    info.case_name = case_path;
    return true;
}

// Conversion is performed entirely by the Python layer.
bool StarCCMAdapter::convert(const std::string& case_path,
                             ConversionResult& result) {
    (void)case_path;
    result.gap_report.unsupported_blocking(
        "conversion", "python_converter",
        "STAR-CCM+ .sim conversion has moved to the Python layer",
        "Use the CFDX Python converter (cfdx_convert CLI) for STAR-CCM+ format parsing");
    return false;
}

}  // namespace starccm
}  // namespace io
}  // namespace cfdx
