#include "PythonBridge.h"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

enum class SourceFormat {
    Unknown,
    SU2,
    Fluent,
    StarCCM,
    Saturne,
    OpenFOAM,
    CFDX
};

SourceFormat detectFormat(const fs::path& path) {
    std::string ext = path.extension().string();
    std::string stem = path.stem().string();
    std::string full = path.filename().string();

    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    std::transform(stem.begin(), stem.end(), stem.begin(), ::tolower);
    std::transform(full.begin(), full.end(), full.begin(), ::tolower);

    if (ext == ".su2" || ext == ".cfg")
        return SourceFormat::SU2;

    if (ext == ".cas" || ext == ".dat" ||
        stem == "cas" || stem == "dat" ||
        full == "cas.h5" || full == "dat.h5" ||
        ext == ".h5") {
        return SourceFormat::Fluent;
    }

    if (ext == ".sim")
        return SourceFormat::StarCCM;

    if (ext == ".xml" || ext == ".py" || ext == ".med" || ext == ".cgns" || ext == ".vtk" || ext == ".vtu")
        return SourceFormat::Saturne;

    if (ext == ".cfdx.h5" || full == ".cfdx.h5")
        return SourceFormat::CFDX;

    if (fs::is_directory(path))
        return SourceFormat::OpenFOAM;

    return SourceFormat::Unknown;
}

const char* formatToString(SourceFormat fmt) {
    switch (fmt) {
        case SourceFormat::SU2: return "SU2";
        case SourceFormat::Fluent: return "Fluent";
        case SourceFormat::StarCCM: return "STAR-CCM+";
        case SourceFormat::Saturne: return "Code_Saturne";
        case SourceFormat::OpenFOAM: return "OpenFOAM";
        case SourceFormat::CFDX: return "CFDX (native)";
        default: return "Unknown";
    }
}

void printUsage(const char* prog) {
    std::cout << "Usage: " << prog << " convert <source> <output>\n\n";
    std::cout << "Convert various CFD formats to CFDX HDF5 (.cfdx.h5)\n\n";
    std::cout << "Examples:\n";
    std::cout << "  " << prog << " convert case.su2 output.cfdx.h5\n";
    std::cout << "  " << prog << " convert case.cas.h5 output.cfdx.h5\n";
    std::cout << "  " << prog << " convert OpenFOAM/ output.cfdx.h5\n";
    std::cout << "  " << prog << " convert case.sim output.cfdx.h5\n";
    std::cout << "  " << prog << " convert case.cfdx.h5 output.cfdx.h5  # normalize/re-export\n";
}

bool validateOutput(const std::string& outputPath);

int main(int argc, char** argv) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    std::string cmd = argv[1];
    if (cmd != "convert") {
        std::cerr << "Unknown command: " << cmd << "\n";
        printUsage(argv[0]);
        return 1;
    }

    if (argc != 4) {
        std::cerr << "Expected 2 arguments for convert: <source> <output>\n";
        printUsage(argv[0]);
        return 1;
    }

    fs::path sourcePath = argv[2];
    fs::path outputPath = argv[3];

    if (!fs::exists(sourcePath)) {
        std::cerr << "Error: Source does not exist: " << sourcePath << "\n";
        return 1;
    }

    SourceFormat fmt = detectFormat(sourcePath);
    if (fmt == SourceFormat::Unknown) {
        std::cerr << "Error: Cannot auto-detect format for: " << sourcePath << "\n";
        return 1;
    }

    std::cout << "CFDX Convert\n";
    std::cout << "  Source: " << sourcePath << "\n";
    std::cout << "  Format: " << formatToString(fmt) << "\n";
    std::cout << "  Output: " << outputPath << "\n\n";

    cfdx::PythonBridge bridge;
    cfdx::PythonBridgeResult result = bridge.convert(sourcePath.string(), outputPath.string());

    if (!result.success) {
        std::cerr << "\nConversion FAILED:\n";
        std::cerr << result.errorMessage << "\n";
        if (!result.stdoutOutput.empty()) {
            std::cerr << "--- Python stdout ---\n" << result.stdoutOutput << "\n";
        }
        if (!result.stderrOutput.empty()) {
            std::cerr << "--- Python stderr ---\n" << result.stderrOutput << "\n";
        }
        return 1;
    }

    std::cout << "\nConversion completed successfully.\n";

    if (!result.summaryJson.empty()) {
        std::cout << "\n--- Conversion Summary ---\n";
        std::cout << result.summaryJson << "\n";
    }

    std::cout << "\nValidating output file...\n";
    if (validateOutput(outputPath.string())) {
        std::cout << "Validation PASSED\n";
    } else {
        std::cerr << "Validation FAILED\n";
        return 1;
    }

    return 0;
}

bool validateOutput(const std::string& outputPath) {
    // Validation is done by the Python converter which writes proper CFDX HDF5 files.
    // The C++ side just does a basic file existence check.
    std::filesystem::path p(outputPath);
    if (!std::filesystem::exists(p)) {
        std::cerr << "Output file does not exist\n";
        return false;
    }
    if (std::filesystem::file_size(p) == 0) {
        std::cerr << "Output file is empty\n";
        return false;
    }
    std::cout << "  Output file: " << outputPath << " (" 
              << std::filesystem::file_size(p) << " bytes)\n";
    return true;
}