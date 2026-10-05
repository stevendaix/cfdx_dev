// Probe time-series CSV export (#426).
//
// Native point probes leave the solver only through a caller-supplied callback,
// and their samples were never retained, so there was no artifact a user could
// export or a downstream tool could read. This writes the samples as CSV.
//
// Schema (version 1), long format — one row per probe per iteration:
//
//   # cfdx-probe-csv v1
//   # columns: probe,iteration,time,value
//   <probe>,<iteration>,<time>,<value>
//
// Long format is used rather than one column per probe because the set of probe
// names is a run-time property: a wide header cannot be written before the
// probes are known, and a sparse wide table is indistinguishable from a missing
// sample. Rows are ordered by probe name then iteration, so the output is
// deterministic and a probe's history is contiguous, which keeps it diffable.
//
// Float encoding: each numeric field is serialized with the shortest
// round-tripping representation using Python repr notation (decimal for
// 10^-4 <= |v| < 10^16, scientific otherwise, with a trailing ".0" for integer
// magnitudes). This is byte-for-byte identical to Python's repr(float), so a
// file written by cfdx::io::write_probe_csv is reproducible verbatim through
// cfdx.probe_validation.write_probe_csv and vice-versa.
//
// The Python layer owns a reader for this exact schema
// (`cfdx.probe_validation.read_probe_csv`), so a file written here round-trips
// through the application layer unchanged.
#pragma once

#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::io {

inline constexpr int ProbeCsvVersion = 1;
inline constexpr const char* ProbeCsvMagic = "cfdx-probe-csv";

// Canonical float serialization for the probe CSV. The emitted text is the
// shortest representation that round-trips exactly, using Python ``repr`` float
// notation: decimal for 10^-4 <= |v| < 10^16, scientific otherwise, with a
// trailing ``.0`` for integer magnitudes. This is byte-for-byte identical to
// ``cfdx.probe_validation.write_probe_csv`` (which uses ``repr(float)``) across
// the full IEEE-754 exponent range, so a CSV produced here is reproducible
// verbatim from the Python layer and vice-versa. The shortest fixed-point form
// is only taken for in-range exponents (safely bounded in length); extreme
// exponents use the shortest scientific form, which never overflows the buffer.
inline std::string format_probe_csv_value(double value)
{
    if (!std::isfinite(value))
        throw std::invalid_argument("probe CSV values must be finite");
    char buffer[64];
    const auto scientific = std::to_chars(buffer, buffer + sizeof(buffer), value,
                                          std::chars_format::scientific);
    std::string text(buffer, scientific.ptr);
    const std::size_t e_pos = text.find_first_of("eE");
    const int exponent = (e_pos != std::string::npos)
        ? std::atoi(text.c_str() + e_pos + 1) : 0;
    if (exponent >= -4 && exponent < 16) {
        const auto fixed = std::to_chars(buffer, buffer + sizeof(buffer), value,
                                         std::chars_format::fixed);
        text.assign(buffer, fixed.ptr);
    }
    if (text.find('.') == std::string::npos &&
        text.find('e') == std::string::npos &&
        text.find('E') == std::string::npos)
        text += ".0";
    return text;
}

// Writes collected probe samples to CSV. Publication is atomic: the file is
// written to a sibling temporary, flushed, and renamed, so a reader never sees a
// half-written export.
inline void write_probe_csv(
    const std::string& path,
    const std::vector<cfdx::physics::IncompressibleProbeSample>& samples)
{
    if (samples.empty())
        throw std::invalid_argument("write_probe_csv: no probe samples to write");

    for (const auto& sample : samples) {
        if (sample.name.empty())
            throw std::invalid_argument("write_probe_csv: probe name must not be empty");
        if (!std::isfinite(sample.time) || !std::isfinite(sample.value))
            throw std::invalid_argument(
                "write_probe_csv: probe '" + sample.name + "' has a non-finite sample");
    }

    std::vector<const cfdx::physics::IncompressibleProbeSample*> ordered;
    ordered.reserve(samples.size());
    for (const auto& sample : samples) ordered.push_back(&sample);
    std::stable_sort(ordered.begin(), ordered.end(),
                     [](const auto* a, const auto* b) {
                         if (a->name != b->name) return a->name < b->name;
                         return a->iteration < b->iteration;
                     });

    const std::string tmp_path = path + ".tmp";
    std::ofstream out(tmp_path);
    if (!out)
        throw std::runtime_error("write_probe_csv: cannot open " + path);
    out << "# " << ProbeCsvMagic << " v" << ProbeCsvVersion << "\n";
    out << "# columns: probe,iteration,time,value\n";
    for (const auto* sample : ordered) {
        // A comma or a newline in a probe name would corrupt the row structure,
        // so the name is validated rather than escaped silently.
        if (sample->name.find(',') != std::string::npos ||
            sample->name.find('\n') != std::string::npos ||
            sample->name.find('\r') != std::string::npos)
            throw std::invalid_argument(
                "write_probe_csv: probe name must not contain a comma or newline: '" +
                sample->name + "'");
        out << sample->name << ',' << sample->iteration << ','
            << format_probe_csv_value(sample->time) << ','
            << format_probe_csv_value(sample->value) << "\n";
    }
    if (!out) {
        out.close();
        std::remove(tmp_path.c_str());
        throw std::runtime_error("write_probe_csv: write failed for " + path);
    }
    out.close();
    if (std::rename(tmp_path.c_str(), path.c_str()) != 0) {
        std::remove(tmp_path.c_str());
        throw std::runtime_error("write_probe_csv: cannot replace " + path);
    }
}

} // namespace cfdx::io