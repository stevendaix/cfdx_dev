// M0.10-T03: GapAnalysis report generation (Markdown + JSON)
//
// Spécification CFDX v0.7 issue #425 §Gap Analysis Engine
// ---------------------------------------------------------------------------
// Generates human-readable Markdown and machine-readable JSON reports
// from the accumulated findings of a conversion.
#include "gap_analysis.h"
#include "io_interface.h"

#include <sstream>
#include <iomanip>

namespace cfdx {
namespace io {

// ---------------------------------------------------------------------------
// JSON serialisation (manual, no external dependency)
// ---------------------------------------------------------------------------
static std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s) {
        switch (c) {
            case '"':  out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n";  break;
            case '\r': out += "\\r";  break;
            case '\t': out += "\\t";  break;
            default:   out += c;      break;
        }
    }
    return out;
}

std::string GapAnalysis::to_json() const {
    std::ostringstream ss;
    ss << std::boolalpha;
    ss << "{\n";
    ss << "  \"has_blocking\": " << has_blocking() << ",\n";
    ss << "  \"summary\": {\n";
    ss << "    \"supported\": " << n_supported() << ",\n";
    ss << "    \"approximated\": " << n_approximated() << ",\n";
    ss << "    \"unsupported_nonblocking\": " << n_unsupported_nonblocking() << ",\n";
    ss << "    \"unsupported_blocking\": " << n_unsupported_blocking() << ",\n";
    ss << "    \"unavailable\": " << n_unavailable() << ",\n";
    ss << "    \"total\": " << findings_.size() << "\n";
    ss << "  },\n";
    ss << "  \"findings\": [\n";
    for (std::size_t i = 0; i < findings_.size(); ++i) {
        const auto& f = findings_[i];
        ss << "    {\n";
        ss << "      \"severity\": \"" << to_string(f.severity) << "\",\n";
        ss << "      \"category\": \"" << json_escape(f.category) << "\",\n";
        ss << "      \"feature\": \"" << json_escape(f.feature) << "\",\n";
        ss << "      \"detail\": \"" << json_escape(f.detail) << "\",\n";
        ss << "      \"suggestion\": \"" << json_escape(f.suggestion) << "\"\n";
        ss << "    }";
        if (i + 1 < findings_.size()) ss << ",";
        ss << "\n";
    }
    ss << "  ]\n";
    ss << "}\n";
    return ss.str();
}

// ---------------------------------------------------------------------------
// Markdown serialisation
// ---------------------------------------------------------------------------
std::string GapAnalysis::to_markdown() const {
    std::ostringstream ss;

    ss << "# Gap Analysis Report\n\n";
    ss << "## Summary\n\n";
    ss << "| Severity | Count |\n";
    ss << "|----------|-------|\n";
    ss << "| Supported / Mapped | " << n_supported() << " |\n";
    ss << "| Approximated | " << n_approximated() << " |\n";
    ss << "| Unsupported (non-blocking) | " << n_unsupported_nonblocking() << " |\n";
    ss << "| Unsupported (blocking) | " << n_unsupported_blocking() << " |\n";
    ss << "| Unavailable (source info missing) | " << n_unavailable() << " |\n";
    ss << "| **Total** | " << findings_.size() << " |\n";
    ss << "\n";
    ss << "**Blocking incompatibilities:** "
       << (has_blocking() ? "YES — conversion blocked" : "NONE") << "\n\n";

    // Group by severity for readability
    static const struct { Severity sev; const char* header; } categories[] = {
        {Severity::SUPPORTED,             "## 1. Supported / Mapped Features\n"},
        {Severity::APPROXIMATED,          "## 2. Approximated / Fallback Mappings\n"},
        {Severity::UNAVAILABLE,           "## 3. Unavailable Source Information\n"},
        {Severity::UNSUPPORTED_NONBLOCK,  "## 4. Unsupported Features (non-blocking)\n"},
        {Severity::UNSUPPORTED_BLOCK,     "## 5. Unsupported Features (blocking)\n"},
    };

    for (const auto& cat : categories) {
        bool any = false;
        for (const auto& f : findings_) {
            if (f.severity == cat.sev) {
                if (!any) {
                    ss << cat.header;
                    ss << "\n| Category | Feature | Detail | Suggestion |\n";
                    ss << "|----------|---------|--------|------------|\n";
                    any = true;
                }
                ss << "| " << f.category << " | " << f.feature
                   << " | " << f.detail << " | " << (f.suggestion.empty() ? "—" : f.suggestion) << " |\n";
            }
        }
        if (!any) {
            ss << cat.header;
            ss << "\n_None._\n";
        }
        ss << "\n";
    }

    return ss.str();
}

}  // namespace io
}  // namespace cfdx
