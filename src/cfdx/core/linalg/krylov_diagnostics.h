#pragma once

// Krylov diagnostic tooling for tracking gauge-row behavior and
// preconditioner failures in coupled Schur solves.
//
// This header is intentionally dependency-light (no .cpp file) so it can be
// included from header-only solver and preconditioner code. All diagnostics
// are controlled by the CFDX_KRYLOV_DIAGNOSTICS environment variable:
//   - unset or "0": diagnostics are compile-time disabled (no runtime cost)
//   - "1": concise per-solve summary lines are emitted
//   - "verbose": per-iteration gauge tracking is emitted

#include <cstdlib>
#include <string>
#include <iostream>

namespace cfdx::core {

enum class KrylovDiagnosticLevel {
    OFF = 0,
    ON = 1,
    VERBOSE = 2
};

inline KrylovDiagnosticLevel krylov_diagnostic_level() {
    static const KrylovDiagnosticLevel level = []() {
        const char* env = std::getenv("CFDX_KRYLOV_DIAGNOSTICS");
        if (!env || env[0] == '\0' || std::string(env) == "0")
            return KrylovDiagnosticLevel::OFF;
        if (std::string(env) == "verbose")
            return KrylovDiagnosticLevel::VERBOSE;
        return KrylovDiagnosticLevel::ON;
    }();
    return level;
}

inline bool krylov_diagnostics_enabled() {
    return krylov_diagnostic_level() > KrylovDiagnosticLevel::OFF;
}

inline bool krylov_diagnostics_verbose() {
    return krylov_diagnostic_level() >= KrylovDiagnosticLevel::VERBOSE;
}

// RAII guard that emits exactly one failure line per scope exit.
// Used to ensure a preconditioner failure is reported once, at the
// root stage, without duplicate generic messages from enclosing code.
class KrylovFailLog {
public:
    KrylovFailLog(const char* stage, const char* context)
        : stage_(stage), context_(context), failed_(false) {}

    ~KrylovFailLog() {
        if (failed_) {
            std::cerr << "KRYLOV_FAIL stage=" << stage_
                      << " context=" << context_ << "\n";
        }
    }

    void mark_failed() { failed_ = true; }

private:
    const char* stage_;
    const char* context_;
    bool failed_;
};

}  // namespace cfdx::core