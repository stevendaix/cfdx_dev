#pragma once

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace cfdx::io {

struct ExecutionSummary {
    int process_exit_code = 0;
    bool converged = false;
    std::size_t iterations = 0;
    int convergence_status = 0;
    std::string convergence_reason;
    bool restart_artifact = false;
    bool convergence_artifact = false;
};

inline void execution_json_string(std::ostream& out, const std::string& value)
{
    out << '"';
    for (const char raw : value) {
        const auto c = static_cast<unsigned char>(raw);
        switch (c) {
        case '"': out << "\\\""; break;
        case '\\': out << "\\\\"; break;
        case '\b': out << "\\b"; break;
        case '\f': out << "\\f"; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default:
            if (c < 0x20U)
                out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                    << static_cast<unsigned int>(c) << std::dec << std::setfill(' ');
            else
                out << static_cast<char>(c);
        }
    }
    out << '"';
}

inline void write_execution_summary(
    const std::string& path,
    const ExecutionSummary& summary)
{
    if (summary.process_exit_code < 0)
        throw std::invalid_argument("execution process_exit_code must be non-negative");

    const std::filesystem::path target(path);
    if (target.empty())
        throw std::invalid_argument("execution summary path must not be empty");

    const auto parent = target.parent_path();
    if (!parent.empty())
        std::filesystem::create_directories(parent);
    const auto temporary = target.string() + ".tmp";

    try {
        std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
        if (!out)
            throw std::runtime_error("cannot open temporary execution summary");

        out << "{\n"
            << "  \"format\": \"CFDX-EXECUTION\",\n"
            << "  \"schema_version\": 1,\n"
            << "  \"process_exit_code\": " << summary.process_exit_code << ",\n"
            << "  \"converged\": " << (summary.converged ? "true" : "false") << ",\n"
            << "  \"iterations\": " << summary.iterations << ",\n"
            << "  \"convergence_status\": " << summary.convergence_status << ",\n"
            << "  \"convergence_reason\": ";
        execution_json_string(out, summary.convergence_reason);
        out << ",\n"
            << "  \"artifacts\": {\n"
            << "    \"restart_dat\": " << (summary.restart_artifact ? "true" : "false") << ",\n"
            << "    \"convergence_json\": " << (summary.convergence_artifact ? "true" : "false") << "\n"
            << "  }\n"
            << "}\n";
        out.close();
        if (!out)
            throw std::runtime_error("failed writing execution summary");

        std::error_code ec;
        std::filesystem::rename(temporary, target, ec);
        if (ec) {
            std::filesystem::remove(target, ec);
            ec.clear();
            std::filesystem::rename(temporary, target, ec);
            if (ec)
                throw std::runtime_error("failed to install execution summary: " + ec.message());
        }
    } catch (...) {
        std::error_code ec;
        std::filesystem::remove(temporary, ec);
        throw;
    }
}

} // namespace cfdx::io
