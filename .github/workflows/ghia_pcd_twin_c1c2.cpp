// ghia_pcd_twin_c1c2.cpp — jumeau #914 étape 2 : variantes d'action PCD à flux non nul.
//
// Branch: ci/long-campaign-diagnostics (never merged, PR #903). Lives under
// .github/workflows/ so pushing it only wakes the ghia-pcd-twin-c1c2 workflow.
//
// Purpose (issue #914): the sandbox/TS twin measured the production action
// -Fp^-1 on the iteration-2 dump: slow geometric outer convergence
// (~0.988/iter, ~1600 projected external iterations). This C++ twin replays
// the EXACT production solve path (real headers: solve_fgmres,
// ILU0Preconditioner, PcdSchurApproximation, CoupledBlockSchurAMGPrecondi-
// tioner) on the bit-exact %a dumps of diag-evidence/ghia-pcd-audit/
// (run 38865832905), with the FULL production budget (restart 4096,
// max_iter 5000, tol 1e-10, warm x0), for four pressure actions:
//
//   Fp  : z_p = -Fp^-1 r_p           (production PROVISIONAL action; inner
//                                     FGMRES/ILU0 512/200/1e-8 — exact mirror)
//   KP1 : z_p = -Kp^-1 r_p           (alternate single inverse; inner
//                                     FGMRES/ILU0 512/200/1e-8 — documented
//                                     deviation from production GMRES/NativeAMG
//                                     40/100, kept so KP1/C1/C2 share one
//                                     inner solver configuration)
//   C1  : z_p = -Mp^-1 Kp Fp^-1 r_p  (exact inverse of the announced Schur
//                                     approximation S ~= -Fp Kp^-1 Mp; one
//                                     solve + one matvec + diagonal Mp)
//   C2  : z_p = -Mp^-1 Fp Kp^-1 r_p  (mirror ordering; one solve + one matvec
//                                     + diagonal Mp)
//
// At zero mass flux Fp ~= Kp and C1/C2 reduce to -Mp^-1 (mass-like, the
// structurally correct Stokes Schur inverse), while -Fp^-1 reduces to
// -Kp^-1 (Laplacian-like): the actions are NOT equal at flux nul — this is
// why the variants are compared on the non-zero-flux iteration 2 dump.
//
// Fail-closed: the parsed matrices are FNV-1a-verified against meta.txt
// before any solve; a checksum mismatch aborts (exit 3). No convergence
// tolerance is relaxed and no MAX_ITER inner solve is silently accepted as
// converged (nested_schur_pressure_solve_usable mirror).
//
// Compile with the production flags (g++ -std=c++17 -O2 -DNDEBUG): the
// acceptance gates are marginal (9.76e-11 vs 1e-10), any FP-rounding
// divergence from the Release build would invalidate the trajectory.

#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace {

using cfdx::core::SparseMatrix;
using cfdx::core::Vector;
using cfdx::core::SolverResult;
using cfdx::core::SolverStatus;
using cfdx::core::KrylovControls;
using cfdx::core::PrecisionPolicy;
using cfdx::core::ILU0Preconditioner;
using cfdx::core::PcdSchurApproximation;
using cfdx::core::CoupledBlockSchurAMGPreconditioner;
using cfdx::core::CoupledBlockSchurOptions;
using cfdx::core::CoupledSchurFactorization;
using cfdx::core::CoupledSchurVelocityApproximation;
using cfdx::core::CoupledSchurApproximationModel;
using cfdx::core::solve_fgmres;

std::uint64_t fnv1a_init() { return 1469598103934665603ULL; }

std::uint64_t fnv1a_bytes(std::uint64_t h, const void* data, std::size_t size)
{
    const auto* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

// %a round-trips every finite double exactly (C99 hex float); C++ stream
// extraction accepts hex floats since C++11 (num_get -> strtod).
bool parse_vec_file(const std::string& path, Vector& out, std::uint64_t& fnv)
{
    std::ifstream in(path);
    if (!in)
        return false;
    std::string line;
    if (!std::getline(in, line)) // CFDX_VEC header
        return false;
    std::vector<double> vals;
    while (std::getline(in, line)) {
        if (line.find_first_not_of(" \t\r\n") == std::string::npos)
            continue;
        // The dump writer uses printf("%a"), e.g. 0x1.8p-4.  The
        // standard iostream num_get implementation on the Actions runner
        // does not reliably parse C99 hexadecimal floating-point tokens;
        // strtod does, and lets us reject trailing junk explicitly.
        char* end = nullptr;
        const double v = std::strtod(line.c_str(), &end);
        if (end == line.c_str() || !std::isfinite(v))
            return false;
        while (*end == ' ' || *end == '\\t' || *end == '\\r' || *end == '\\n')
            ++end;
        if (*end != '\\0')
            return false;
        vals.push_back(v);
    }
    fnv = fnv1a_init();
    fnv = fnv1a_bytes(fnv, vals.data(), vals.size() * sizeof(double));
    if (out.size() != vals.size())
        out.resize(vals.size());
    for (std::size_t i = 0; i < vals.size(); ++i)
        out(i) = vals[i];
    return true;
}

// Parse a key=value token like "nnz=76991".
bool kv_size(const std::string& tok, std::size_t& out)
{
    const std::size_t eq = tok.find('=');
    if (eq == std::string::npos)
        return false;
    out = static_cast<std::size_t>(
        std::strtoull(tok.c_str() + eq + 1, nullptr, 10));
    return true;
}

// Chunked CSR: prefix.NNNN files (dump_csr_chunks format). Each chunk carries
// the GLOBAL row_offsets slice for its row range plus the columns and values
// of those rows; the concatenation of all chunks in order is the full CSR
// triple, duplicates preserved.
bool parse_csr_chunks(const std::string& prefix, SparseMatrix& M,
                      std::uint64_t& fnv)
{
    std::vector<std::size_t> off;
    std::vector<std::uint32_t> cols;
    std::vector<double> vals;
    std::size_t n_rows = 0, n_cols = 0, total_nnz = 0;
    std::size_t expected_chunks = 0, seen = 0;
    for (std::size_t c = 0;; ++c) {
        char idx[32];
        std::snprintf(idx, sizeof(idx), "%04llu",
                      static_cast<unsigned long long>(c));
        std::ifstream in(prefix + "." + idx);
        if (!in)
            break;
        ++seen;
        std::string line;
        std::size_t row0 = 0, row1 = 0;
        if (!std::getline(in, line))
            return false; // CFDX_CSR_CHUNK header
        {
            std::istringstream h(line);
            std::string tok;
            while (h >> tok) {
                if (tok.rfind("chunks=", 0) == 0) {
                    if (!kv_size(tok, expected_chunks))
                        return false;
                } else if (tok.rfind("row0=", 0) == 0) {
                    if (!kv_size(tok, row0))
                        return false;
                } else if (tok.rfind("row_end=", 0) == 0) {
                    if (!kv_size(tok, row1))
                        return false;
                } else if (tok.rfind("n_rows=", 0) == 0) {
                    if (!kv_size(tok, n_rows))
                        return false;
                } else if (tok.rfind("n_cols=", 0) == 0) {
                    if (!kv_size(tok, n_cols))
                        return false;
                } else if (tok.rfind("nnz=", 0) == 0) {
                    if (!kv_size(tok, total_nnz))
                        return false;
                }
                // "CFDX_CSR_CHUNK", "name=...", "chunk=..." : skipped.
            }
        }
        if (n_rows == 0 || row1 <= row0 || row1 > n_rows)
            return false;
        if (off.empty())
            off.assign(n_rows + 1, 0);
        if (!std::getline(in, line) || line.rfind("OFFSETS", 0) != 0)
            return false;
        {
            std::istringstream ls(line);
            std::string tag;
            ls >> tag; // OFFSETS
            for (std::size_t i = row0; i <= row1; ++i) {
                std::size_t v = 0;
                if (!(ls >> v))
                    return false;
                off[i] = v;
            }
        }
        const std::size_t chunk_nnz = off[row1] - off[row0];
        if (!std::getline(in, line) || line.rfind("COLUMNS", 0) != 0)
            return false;
        {
            std::istringstream ls(line);
            std::string tag;
            ls >> tag; // COLUMNS
            std::size_t count = 0;
            std::size_t v = 0;
            while (ls >> v) {
                cols.push_back(static_cast<std::uint32_t>(v));
                ++count;
            }
            if (count != chunk_nnz)
                return false;
        }
        if (!std::getline(in, line) || line.rfind("VALUES", 0) != 0)
            return false;
        {
            std::istringstream ls(line);
            std::string tag;
            std::size_t count = 0;
            ls >> tag >> count;
            if (count != chunk_nnz)
                return false;
            for (std::size_t k = 0; k < count; ++k) {
                if (!std::getline(in, line))
                    return false;
                std::istringstream vs(line);
                double x = 0.0;
                if (!(vs >> x))
                    return false;
                vals.push_back(x);
            }
        }
    }
    if (seen == 0 || (expected_chunks != 0 && seen != expected_chunks))
        return false;
    if (n_rows == 0 || n_cols == 0 || vals.size() != total_nnz ||
        cols.size() != total_nnz)
        return false;
    fnv = fnv1a_init();
    fnv = fnv1a_bytes(fnv, vals.data(), vals.size() * sizeof(double));
    // Rebuild through the same SparseMatrix insertion path as production
    // assembly, in the exact dumped order (duplicates preserved).
    M = SparseMatrix(n_rows, n_cols);
    for (std::size_t i = 0; i < n_rows; ++i)
        for (std::size_t k = off[i]; k < off[i + 1]; ++k)
            M.push_back(i, cols[k], vals[k]);
    M.finalize();
    return M.n_rows() == n_rows && M.n_cols() == n_cols;
}

void matvec(const SparseMatrix& M, const Vector& x, Vector& y)
{
    const auto* ro = M.row_offsets_data();
    const auto* co = M.columns_data();
    const double* va = M.values_data();
    for (std::size_t i = 0; i < M.n_rows(); ++i) {
        double s = 0.0;
        for (std::size_t k = ro[i]; k < ro[i + 1]; ++k)
            s += va[k] * x(co[k]);
        y(i) = s;
    }
}

// Verbatim mirror of nested_schur_pressure_solve_usable (solver_fix.h):
// a MAX_ITER inner solve is usable only with a finite iterate (documented
// best-effort); any other status is a hard failure.
bool pressure_solve_usable(const SolverResult& r, const Vector& x,
                           std::size_t expected_size)
{
    if (r.status == SolverStatus::CONVERGED)
        return true;
    if (r.status != SolverStatus::MAX_ITER_REACHED)
        return false;
    if (x.size() != expected_size || !x.is_valid())
        return false;
    for (std::size_t i = 0; i < x.size(); ++i)
        if (!std::isfinite(x(i)))
            return false;
    return true;
}

const char* status_name(SolverStatus s)
{
    switch (s) {
    case SolverStatus::CONVERGED: return "CONVERGED";
    case SolverStatus::MAX_ITER_REACHED: return "MAX_ITER_REACHED";
    case SolverStatus::DIVERGED: return "DIVERGED";
    case SolverStatus::NOT_APPLICABLE: return "NOT_APPLICABLE";
    default: return "UNKNOWN";
    }
}

struct InnerStats {
    std::size_t applies = 0, total_iters = 0, max_iters = 0;
    std::size_t converged = 0, best_effort = 0, fail = 0;
};

int run_variant(const std::string& name,
                const SparseMatrix& A, const Vector& b, const Vector& x0,
                const SparseMatrix& Mp, const SparseMatrix& Kp,
                const SparseMatrix& Fp, const Vector& mpdiag,
                std::size_t nc, std::size_t reference_cell)
{
    InnerStats st;
    auto ilu_fp = std::make_unique<ILU0Preconditioner>();
    if (!ilu_fp->setup(Fp)) {
        std::fprintf(stderr, "TWIN_SETUP_FAIL name=%s op=Fp ilu0\n",
                     name.c_str());
        return 4;
    }
    std::unique_ptr<ILU0Preconditioner> ilu_kp;
    if (name != "Fp") {
        ilu_kp = std::make_unique<ILU0Preconditioner>();
        if (!ilu_kp->setup(Kp)) {
            std::fprintf(stderr, "TWIN_SETUP_FAIL name=%s op=Kp ilu0\n",
                         name.c_str());
            return 4;
        }
    }

    // Production mirror of solve_fp (solver_fix.h): FGMRES(Fp, restart
    // min(512, n), budget 200, tol 1e-8, prepared ILU0, default controls).
    auto solve_fp = [&](const Vector& rhs, Vector& x) -> bool {
        const int restart = std::min<int>(
            512, static_cast<int>(Fp.n_rows()));
        SolverResult r = solve_fgmres(Fp, rhs, x, restart, 200, 1e-8,
                                      ilu_fp.get(), KrylovControls{},
                                      PrecisionPolicy{}, false);
        st.applies++;
        st.total_iters += r.iterations;
        if (r.iterations > st.max_iters)
            st.max_iters = r.iterations;
        if (r.status == SolverStatus::CONVERGED) {
            st.converged++;
        } else {
            if (r.status == SolverStatus::MAX_ITER_REACHED)
                st.best_effort++;
            std::fprintf(stderr,
                          "COUPLED_NESTED_PCD_BEST_EFFORT variant=%s"
                          " operator=Fp status=%d iterations=%llu"
                          " residual=%.17g relative=%.17g\n",
                          name.c_str(), static_cast<int>(r.status),
                          static_cast<unsigned long long>(r.iterations),
                          r.residual, r.residual_relative);
        }
        const bool ok = pressure_solve_usable(r, x, Fp.n_rows());
        if (!ok)
            st.fail++;
        return ok;
    };
    // Variant inner Kp solve: FGMRES/ILU0 512/200/1e-8. Documented deviation
    // from the production Kp path (GMRES/NativeAMG 40/100) so that KP1, C1
    // and C2 share one inner solver configuration.
    auto solve_kp = [&](const Vector& rhs, Vector& x) -> bool {
        const int restart = std::min<int>(
            512, static_cast<int>(Kp.n_rows()));
        SolverResult r = solve_fgmres(Kp, rhs, x, restart, 200, 1e-8,
                                      ilu_kp.get(), KrylovControls{},
                                      PrecisionPolicy{}, false);
        st.applies++;
        st.total_iters += r.iterations;
        if (r.iterations > st.max_iters)
            st.max_iters = r.iterations;
        if (r.status == SolverStatus::CONVERGED) {
            st.converged++;
        } else {
            if (r.status == SolverStatus::MAX_ITER_REACHED)
                st.best_effort++;
            std::fprintf(stderr,
                          "COUPLED_NESTED_PCD_BEST_EFFORT variant=%s"
                          " operator=Kp status=%d iterations=%llu"
                          " residual=%.17g relative=%.17g\n",
                          name.c_str(), static_cast<int>(r.status),
                          static_cast<unsigned long long>(r.iterations),
                          r.residual, r.residual_relative);
        }
        const bool ok = pressure_solve_usable(r, x, Kp.n_rows());
        if (!ok)
            st.fail++;
        return ok;
    };

    // The convection_diffusion_solve callback IS the variant action: after it
    // returns, PcdSchurApproximation::apply negates the output, so
    //     action = -(callback output).
    PcdSchurApproximation::PressureSolve inner;
    Vector tmp(nc, 0.0);
    if (name == "Fp") {
        inner = solve_fp; // action -Fp^-1
    } else if (name == "KP1") {
        inner = solve_kp; // action -Kp^-1
    } else if (name == "C1") {
        // callback = Mp^-1 Kp Fp^-1 -> action -Mp^-1 Kp Fp^-1
        inner = [&](const Vector& rhs, Vector& x) -> bool {
            for (std::size_t i = 0; i < nc; ++i)
                tmp(i) = 0.0;
            if (!solve_fp(rhs, tmp))
                return false;
            matvec(Kp, tmp, x);
            for (std::size_t i = 0; i < nc; ++i)
                x(i) = x(i) / mpdiag(i);
            return true;
        };
    } else if (name == "C2") {
        // callback = Mp^-1 Fp Kp^-1 -> action -Mp^-1 Fp Kp^-1
        inner = [&](const Vector& rhs, Vector& x) -> bool {
            for (std::size_t i = 0; i < nc; ++i)
                tmp(i) = 0.0;
            if (!solve_kp(rhs, tmp))
                return false;
            matvec(Fp, tmp, x);
            for (std::size_t i = 0; i < nc; ++i)
                x(i) = x(i) / mpdiag(i);
            return true;
        };
    } else {
        return 1;
    }

    // Production mirror of the coupled PCD Schur construction (solver_fix.h):
    // Full/Block/PCD options, reference cell 528 (no fixed pressure boundary
    // on this case, so the gauge row is the reference row).
    auto pcd = std::make_unique<PcdSchurApproximation>(
        Mp, Kp, Fp, solve_kp, inner, std::nullopt,
        std::optional<std::size_t>(reference_cell));
    CoupledBlockSchurOptions options;
    options.factorization = CoupledSchurFactorization::Full;
    options.velocity_approximation =
        CoupledSchurVelocityApproximation::Block;
    options.schur_approximation = CoupledSchurApproximationModel::PCD;
    auto schur = std::make_unique<CoupledBlockSchurAMGPreconditioner>(
        nc, options, reference_cell);
    schur->set_pcd_schur(std::move(pcd));
    if (!schur->setup(A)) {
        std::fprintf(stderr, "TWIN_SETUP_FAIL name=%s coupled error=%s\n",
                     name.c_str(), schur->last_error().c_str());
        return 4;
    }

    // Production mirror of the outer call: FGMRES(A, b, x, 4096, 5000,
    // 1e-10, precond, controls{restart_min=restart_max=4096, adaptive=0}),
    // warm start x0, default setup_preconditioner=true (the driver's second
    // setup is clean per the production comment; reset preserves pcd_schur_).
    Vector x(x0.size(), 0.0);
    for (std::size_t i = 0; i < x0.size(); ++i)
        x(i) = x0(i);
    KrylovControls controls;
    controls.restart_min = 4096;
    controls.restart_max = 4096;
    controls.adaptive_restart = false;
    const auto t0 = std::chrono::steady_clock::now();
    SolverResult r = solve_fgmres(A, b, x, 4096, 5000, 1e-10, schur.get(),
                                  controls);
    const double dt =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
            .count();
    std::fprintf(stderr,
                 "TWIN_VARIANT name=%s status=%s status_int=%d"
                 " iterations=%llu residual=%.17g relative=%.17g"
                 " x_gauge=%.17g elapsed_s=%.3f\n",
                 name.c_str(), status_name(r.status),
                 static_cast<int>(r.status),
                 static_cast<unsigned long long>(r.iterations), r.residual,
                 r.residual_relative, x(3 * nc + reference_cell), dt);
    std::fprintf(stderr,
                 "TWIN_INNER name=%s applies=%llu total_iters=%llu"
                 " avg=%.2f max=%llu converged=%llu best_effort=%llu"
                 " fail=%llu\n",
                 name.c_str(),
                 static_cast<unsigned long long>(st.applies),
                 static_cast<unsigned long long>(st.total_iters),
                 st.applies ? double(st.total_iters) / double(st.applies)
                            : 0.0,
                 static_cast<unsigned long long>(st.max_iters),
                 static_cast<unsigned long long>(st.converged),
                 static_cast<unsigned long long>(st.best_effort),
                 static_cast<unsigned long long>(st.fail));
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    const std::string dir =
        argc > 1 ? argv[1] : "diag-evidence/ghia-pcd-audit";
    const std::size_t nc = 1024;
    const std::size_t n = 4 * nc;
    const std::size_t reference_cell = 528;

    SparseMatrix A, Mp, Kp, Fp;
    std::uint64_t fnv = 0;
    if (!parse_csr_chunks(dir + "/A_iter2.csr", A, fnv)) {
        std::fprintf(stderr, "TWIN_PARSE_FAIL file=A_iter2.csr\n");
        return 3;
    }
    const std::uint64_t fnvA = fnv;
    if (!parse_csr_chunks(dir + "/Mp.csr", Mp, fnv)) {
        std::fprintf(stderr, "TWIN_PARSE_FAIL file=Mp.csr\n");
        return 3;
    }
    const std::uint64_t fnvMp = fnv;
    if (!parse_csr_chunks(dir + "/Kp.csr", Kp, fnv)) {
        std::fprintf(stderr, "TWIN_PARSE_FAIL file=Kp.csr\n");
        return 3;
    }
    const std::uint64_t fnvKp = fnv;
    if (!parse_csr_chunks(dir + "/Fp.csr", Fp, fnv)) {
        std::fprintf(stderr, "TWIN_PARSE_FAIL file=Fp.csr\n");
        return 3;
    }
    const std::uint64_t fnvFp = fnv;
    Vector b(n, 0.0), x0(n, 0.0);
    if (!parse_vec_file(dir + "/b_iter2.vec", b, fnv)) {
        std::fprintf(stderr, "TWIN_PARSE_FAIL file=b_iter2.vec\n");
        return 3;
    }
    const std::uint64_t fnvB = fnv;
    if (!parse_vec_file(dir + "/x0_iter2.vec", x0, fnv)) {
        std::fprintf(stderr, "TWIN_PARSE_FAIL file=x0_iter2.vec\n");
        return 3;
    }
    const std::uint64_t fnvX0 = fnv;

    // Fail-closed checksum verification against meta.txt (run 38865832905).
    struct Check {
        const char* name;
        std::uint64_t got;
        const char* expected;
    };
    const Check checks[] = {
        {"A_iter2", fnvA, "f1c712aa91af8a05"},
        {"Mp", fnvMp, "8e8450068b1f4383"},
        {"Kp", fnvKp, "a55ea9de9852207a"},
        {"Fp", fnvFp, "482578abba20a0ba"},
        {"b_iter2", fnvB, "ba89d2d797327475"},
        {"x0_iter2", fnvX0, "62e4545f548ca62c"},
    };
    char got[32];
    bool all_ok = true;
    for (const auto& c : checks) {
        std::snprintf(got, sizeof(got), "%016llx",
                      static_cast<unsigned long long>(c.got));
        const bool ok = std::string(got) == c.expected;
        all_ok = all_ok && ok;
        std::fprintf(stderr,
                     "TWIN_CHECKSUM name=%s got=%s expected=%s ok=%d\n",
                     c.name, got, c.expected, ok ? 1 : 0);
    }
    if (!all_ok) {
        std::fprintf(stderr, "TWIN_EXIT_CODE=3 checksum mismatch — abort\n");
        return 3;
    }

    // Structural checks.
    std::fprintf(stderr, "TWIN_DIMS n=%llu nnz_A=%llu gauge_row=%llu\n",
                 static_cast<unsigned long long>(A.n_rows()),
                 static_cast<unsigned long long>(A.nnz()),
                 static_cast<unsigned long long>(3 * nc + reference_cell));
    {
        const auto* ro = A.row_offsets_data();
        const auto* co = A.columns_data();
        const double* va = A.values_data();
        const std::size_t gr = 3 * nc + reference_cell;
        const bool gauge_ok = ro[gr + 1] - ro[gr] == 1 && co[ro[gr]] == gr &&
                              va[ro[gr]] == 1.0;
        std::fprintf(stderr, "TWIN_GAUGE row_identity=%d x0_gauge=%.17g\n",
                     gauge_ok ? 1 : 0, x0(gr));
        if (!gauge_ok) {
            std::fprintf(stderr, "TWIN_EXIT_CODE=3 gauge row invalid\n");
            return 3;
        }
    }
    // Mp diagonal extraction (uniform-mesh mass operator check).
    Vector mpdiag(nc, 0.0);
    {
        const auto* ro = Mp.row_offsets_data();
        const auto* co = Mp.columns_data();
        const double* va = Mp.values_data();
        double mn = 1e300, mx = -1e300;
        for (std::size_t i = 0; i < nc; ++i) {
            bool found = false;
            for (std::size_t k = ro[i]; k < ro[i + 1]; ++k) {
                if (co[k] == i) {
                    mpdiag(i) = va[k];
                    found = true;
                }
            }
            if (!found || mpdiag(i) == 0.0 || !std::isfinite(mpdiag(i))) {
                std::fprintf(stderr,
                             "TWIN_EXIT_CODE=3 Mp diagonal missing at %llu\n",
                             static_cast<unsigned long long>(i));
                return 3;
            }
            mn = std::min(mn, mpdiag(i));
            mx = std::max(mx, mpdiag(i));
        }
        std::fprintf(stderr, "TWIN_MP_DIAG min=%.17g max=%.17g uniform=%d\n",
                     mn, mx, (mn == mx) ? 1 : 0);
    }

    // Variants, in protocol order: production baseline first, then C1
    // (recorded before C2), then C2, then the KP1 single-inverse alternate.
    const char* order[] = {"Fp", "C1", "C2", "KP1"};
    for (const char* v : order) {
        const int rc = run_variant(v, A, b, x0, Mp, Kp, Fp, mpdiag, nc,
                                   reference_cell);
        if (rc != 0) {
            std::fprintf(stderr, "TWIN_EXIT_CODE=%d variant=%s\n", rc, v);
            return rc;
        }
    }
    std::fprintf(stderr, "TWIN_SUMMARY variants=4 all_completed=1\n");
    std::fprintf(stderr, "TWIN_EXIT_CODE=0\n");
    return 0;
}
