#pragma once

#include <string>
#include <optional>

namespace cfdx {

struct PythonBridgeResult {
    bool success = false;
    std::string errorMessage;
    std::string stdoutOutput;
    std::string stderrOutput;
    std::string summaryJson;
};

class PythonBridge {
public:
    PythonBridge() = default;
    ~PythonBridge() = default;

    PythonBridgeResult convert(const std::string& sourcePath, const std::string& outputPath);

private:
    std::optional<std::string> findPythonExecutable() const;
    std::optional<std::string> findConverterScript() const;
};

} // namespace cfdx