// Unit tests for the native probe CSV export (cfdx::io::write_probe_csv).
//
// Pins the on-disk schema shared with the Python reader
// (``cfdx.probe_validation.read_probe_csv``): the magic/version header, the
// column contract, deterministic row ordering, exact float round-trip, and the
// input validations that keep a malformed export from being written.
#include "cfdx/io/probe/probe_csv.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace cfdx::io;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {

std::string read_all(const std::filesystem::path& path)
{
    std::ifstream in(path);
    EXPECT_TRUE(in.is_open());
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
}

std::filesystem::path unique_temp_path(const char* suffix)
{
    static std::uint64_t counter = 0;
    auto base = std::filesystem::temp_directory_path() /
        ("cfdx_probe_csv_test_" + std::to_string(counter++) + suffix);
    std::filesystem::remove(base);
    return base;
}

std::vector<std::string> split_lines(const std::string& text)
{
    std::vector<std::string> lines;
    std::string current;
    for (char c : text) {
        if (c == '\n') {
            lines.push_back(current);
            current.clear();
        } else {
            current.push_back(c);
        }
    }
    if (!current.empty())
        lines.push_back(current);
    return lines;
}

} // namespace

int main()
{
    run_case("writes_native_schema_header", []() {
        const auto path = unique_temp_path(".csv");
        write_probe_csv(path.string(), {
            {IncompressibleProbeSample{"pressure", 1, 0.0, 1.0}},
        });
        const std::string text = read_all(path);
        std::vector<std::string> lines = split_lines(text);
        EXPECT_TRUE(lines.size() >= 3);
        EXPECT_TRUE(lines[0] == "# cfdx-probe-csv v1");
        EXPECT_TRUE(lines[1] == "# columns: probe,iteration,time,value");
        EXPECT_TRUE(lines[2] == "pressure,1,0.0,1.0");
    });

    run_case("orders_rows_by_name_then_iteration", []() {
        const auto path = unique_temp_path(".csv");
        write_probe_csv(path.string(), {
            {IncompressibleProbeSample{"zeta", 2, 0.2, 3.0}},
            {IncompressibleProbeSample{"alpha", 2, 0.2, 2.0}},
            {IncompressibleProbeSample{"alpha", 1, 0.1, 1.0}},
            {IncompressibleProbeSample{"zeta", 1, 0.1, 4.0}},
        });
        const std::string text = read_all(path);
        std::vector<std::string> lines = split_lines(text);
        EXPECT_TRUE(lines[0] == "# cfdx-probe-csv v1");
        EXPECT_TRUE(lines[1] == "# columns: probe,iteration,time,value");
        EXPECT_TRUE(lines[2] == "alpha,1,0.1,1.0");
        EXPECT_TRUE(lines[3] == "alpha,2,0.2,2.0");
        EXPECT_TRUE(lines[4] == "zeta,1,0.1,4.0");
        EXPECT_TRUE(lines[5] == "zeta,2,0.2,3.0");
    });

    run_case("round_trips_full_precision_values", []() {
        const auto path = unique_temp_path(".csv");
        const double third = 1.0 / 3.0;
        const double tiny = 1e-9;
        write_probe_csv(path.string(), {
            {IncompressibleProbeSample{"p", 1, 0.1, third}},
            {IncompressibleProbeSample{"p", 2, tiny, third * 10.0}},
        });
        const std::string text = read_all(path);
        std::vector<std::string> lines = split_lines(text);
        // C++ to_chars shortest must match the canonical float spelling that
        // Python's repr() would emit: a decimal point is always present.
        EXPECT_TRUE(lines[2].find(",0.1,") != std::string::npos);
        EXPECT_TRUE(lines[2].back() != ',');
        // The value column must round-trip to the identical double.
        std::string value_text = lines[2].substr(lines[2].rfind(',') + 1);
        EXPECT_TRUE(std::stod(value_text) == third);
    });

    run_case("atomic_publish_leaves_no_tmp", []() {
        const auto path = unique_temp_path(".csv");
        write_probe_csv(path.string(), {
            {IncompressibleProbeSample{"u_x", 1, 0.0, 0.0}},
        });
        EXPECT_TRUE(std::filesystem::exists(path));
        EXPECT_TRUE(!std::filesystem::exists(path.string() + ".tmp"));
    });

    run_case("rejects_empty_samples", []() {
        const auto path = unique_temp_path(".csv");
        EXPECT_THROW_WITH(write_probe_csv(path.string(), {}),
                          std::invalid_argument, "no probe samples to write");
        EXPECT_TRUE(!std::filesystem::exists(path));
    });

    run_case("rejects_empty_probe_name", []() {
        const auto path = unique_temp_path(".csv");
        EXPECT_THROW_WITH(write_probe_csv(path.string(), {
            {IncompressibleProbeSample{"", 1, 0.0, 1.0}},
        }), std::invalid_argument, "probe name must not be empty");
    });

    run_case("rejects_nonfinite_value", []() {
        const auto path = unique_temp_path(".csv");
        EXPECT_THROW_WITH(write_probe_csv(path.string(), {
            {IncompressibleProbeSample{"p", 1, 0.0,
                                       std::numeric_limits<double>::quiet_NaN()}},
        }), std::invalid_argument, "non-finite");
    });

    run_case("rejects_nonfinite_time", []() {
        const auto path = unique_temp_path(".csv");
        EXPECT_THROW_WITH(write_probe_csv(path.string(), {
            {IncompressibleProbeSample{"p", 1,
                                       std::numeric_limits<double>::infinity(), 1.0}},
        }), std::invalid_argument, "non-finite");
    });

    run_case("rejects_name_with_comma", []() {
        const auto path = unique_temp_path(".csv");
        EXPECT_THROW_WITH(write_probe_csv(path.string(), {
            {IncompressibleProbeSample{"bad,name", 1, 0.0, 1.0}},
        }), std::invalid_argument, "comma or newline");
    });

    return run_all();
}
