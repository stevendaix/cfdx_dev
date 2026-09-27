#include "PythonBridge.h"
#include "cfdx/io/hdf5/case_hdf5_io.h"
#include "cfdx/io/cfdx_io/case_schema.h"
#include "cfdx/io/cfdx_io/gap_analysis.h"
#include "cfdx/core/mesh/mesh.h"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

void printUsage(const char* prog) {
    std::cout << "Usage: " << prog << " convert <source> <output>\\n\\n";
    std::cout << "Convert solver cases to the canonical CFDX HDF5 format.\\n";
    std::cout << "Format detection and conversion are delegated to the Python layer.\\n";
}


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
    if (std::string(argv[1]) != "convert" || argc != 4) {
        std::cerr << "Expected: " << argv[0] << " convert <source> <output>\\n";
        printUsage(argv[0]);
        return 1;
    }

    const fs::path sourcePath = argv[2];
    const fs::path outputPath = argv[3];
    if (!fs::exists(sourcePath)) {
        std::cerr << "Error: source does not exist: " << sourcePath << "\\n";
        return 1;
    }

    std::cout << "CFDX Convert\\n"
              << "  Source: " << sourcePath << "\\n"
              << "  Output: " << outputPath << "\\n\\n";

    cfdx::PythonBridge bridge;
    const auto result = bridge.convert(sourcePath.string(), outputPath.string());
    if (!result.success) {
        std::cerr << "Conversion FAILED: " << result.errorMessage << "\\n";
        if (!result.stdoutOutput.empty())
            std::cerr << "--- Python output ---\\n" << result.stdoutOutput << "\\n";
        return 1;
    }

    cfdx::core::Mesh mesh;
    cfdx::io::SourceInfo source;
    cfdx::io::CaseSetup setup;
    cfdx::io::GapAnalysis gap;
    if (!cfdx::io::read_case_cfdx_h5(outputPath.string(), mesh, source, setup, gap)) {
        std::cerr << "Validation FAILED: output is not a valid CFDX HDF5 case.\\n";
        return 1;
    }

    std::cout << "Validation PASSED: CFDX HDF5 schema and mesh are readable.\\n";
    return 0;
}

