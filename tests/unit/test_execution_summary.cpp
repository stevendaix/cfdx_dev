#include "cfdx/io/runtime/execution_summary.h"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

int main()
{
    using cfdx::io::ExecutionSummary;
    using cfdx::io::write_execution_summary;

    const auto path = std::filesystem::temp_directory_path() /
                      "cfdx_execution_summary_test.json";
    std::error_code ec;
    std::filesystem::remove(path, ec);

    write_execution_summary(path.string(), ExecutionSummary{
        0, true, 7, 1, "converged \"normally\"", true, true});

    std::ifstream in(path);
    assert(in.good());
    std::stringstream buffer;
    buffer << in.rdbuf();
    const std::string json = buffer.str();

    assert(json.find("\"format\": \"CFDX-EXECUTION\"") != std::string::npos);
    assert(json.find("\"schema_version\": 1") != std::string::npos);
    assert(json.find("\"process_exit_code\": 0") != std::string::npos);
    assert(json.find("\"converged\": true") != std::string::npos);
    assert(json.find("\"iterations\": 7") != std::string::npos);
    assert(json.find("\\\"normally\\\"") != std::string::npos);
    assert(json.find("\"restart_dat\": true") != std::string::npos);
    assert(json.find("\"convergence_json\": true") != std::string::npos);

    write_execution_summary(path.string(), ExecutionSummary{
        1, false, 9, 0, "iteration limit", true, true});
    in.close();
    in.open(path);
    buffer.str({});
    buffer.clear();
    buffer << in.rdbuf();
    assert(buffer.str().find("\"process_exit_code\": 1") != std::string::npos);
    assert(buffer.str().find("\"converged\": false") != std::string::npos);

    std::filesystem::remove(path, ec);
    std::filesystem::remove(path.string() + ".tmp", ec);
    return 0;
}
