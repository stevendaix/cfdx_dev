#include "PythonBridge.h"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace cfdx {

std::optional<std::string> PythonBridge::findPythonExecutable() const {
    const char* envPython = std::getenv("CFDX_PYTHON_EXECUTABLE");
    if (envPython && fs::exists(envPython)) {
        return envPython;
    }

    std::array<const char*, 5> candidates = {
        "python3",
        "python",
        "/usr/bin/python3",
        "/usr/local/bin/python3",
        "/opt/homebrew/bin/python3"
    };

    for (const char* cmd : candidates) {
        std::array<char, 128> buffer;
        std::string whichCmd = "which " + std::string(cmd);
        std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(whichCmd.c_str(), "r"), pclose);
        if (pipe && std::fgets(buffer.data(), buffer.size(), pipe.get())) {
            std::string path = buffer.data();
            path.erase(path.find_last_not_of("\n\r") + 1);
            if (fs::exists(path)) {
                return path;
            }
        }
    }

    return std::nullopt;
}

std::optional<std::string> PythonBridge::findConverterScript() const {
    const char* envScript = std::getenv("CFDX_CONVERT_SCRIPT");
    if (envScript && fs::exists(envScript)) {
        return envScript;
    }

    fs::path baseDir = fs::current_path();
    std::vector<fs::path> searchPaths = {
        baseDir / "src" / "cfdx" / "python" / "cfdx" / "cli" / "convert_cli.py",
        baseDir / "python" / "cfdx" / "cli" / "convert_cli.py",
        baseDir / "cfdx" / "cli" / "convert_cli.py",
        baseDir / "convert_cli.py"
    };

    for (const auto& path : searchPaths) {
        if (fs::exists(path)) {
            return path.string();
        }
    }

    return std::nullopt;
}

PythonBridgeResult PythonBridge::convert(const std::string& sourcePath, const std::string& outputPath) {
    PythonBridgeResult result;

    auto pythonExe = findPythonExecutable();
    if (!pythonExe) {
        result.errorMessage = "Python interpreter not found. Set CFDX_PYTHON_EXECUTABLE or ensure 'python3' is in PATH.";
        return result;
    }

    auto scriptPath = findConverterScript();
    if (!scriptPath) {
        result.errorMessage = "Converter script not found. Set CFDX_CONVERT_SCRIPT or ensure src/cfdx/python/cfdx/cli/convert_cli.py exists.";
        return result;
    }

    std::vector<std::string> cmd = {
        *pythonExe,
        *scriptPath,
        "convert",
        sourcePath,
        outputPath
    };

    std::string cmdStr;
    for (const auto& arg : cmd) {
        if (!cmdStr.empty()) cmdStr += " ";
        cmdStr += arg;
    }

    std::array<char, 4096> buffer;
    std::string fullCmd = cmdStr + " 2>&1";

    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(fullCmd.c_str(), "r"), pclose);
    if (!pipe) {
        result.errorMessage = "Failed to launch Python subprocess: " + cmdStr;
        return result;
    }

    std::ostringstream output;
    while (std::fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        output << buffer.data();
    }

    int exitCode = pclose(pipe.release());
    if (exitCode == -1) {
        result.errorMessage = "Failed to wait for Python subprocess";
        return result;
    }

    result.stdoutOutput = output.str();
    result.success = (WIFEXITED(exitCode) && WEXITSTATUS(exitCode) == 0);

    if (!result.success) {
        result.errorMessage = "Python converter exited with code " + std::to_string(WEXITSTATUS(exitCode));
    }

    size_t jsonStart = result.stdoutOutput.rfind("{");
    if (jsonStart != std::string::npos) {
        result.summaryJson = result.stdoutOutput.substr(jsonStart);
    }

    return result;
}

} // namespace cfdx