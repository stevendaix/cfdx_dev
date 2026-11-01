#!/usr/bin/env python3
"""Apply PR #1149 forensics-v2 evidence-only edits to a clean worktree.

Fail-closed: every anchor must match exactly once, and the run aborts
without writing anything if an anchor is missing, duplicated, or if a
previous application is detected. No tolerance, gate, oracle or numerical
path is modified by these splices.
"""

from __future__ import annotations

import sys
from pathlib import Path

SOLVER = Path("src/cfdx/physics/steady_incompressible_solver.h")
CMAKE = Path("CMakeLists.txt")
DOC = Path("docs/validation/NACA0012_RUNTIME_FORENSIC_AUDIT_2026-10-31.md")
APPENDIX = (
    Path(sys.argv[1])
    if len(sys.argv) > 1
    else Path("diag-tools/pr1149_naca_forensics_v2/docs_appendix.md")
)
APPENDIX_MARKER = "## Chronology, stencil quality and iteration-62 snapshot"
STENCIL_MARKER = "NACA_LS_STENCIL_QUALITY"

# Single splice: one-shot least-squares stencil audit (evidence only, gated by
# the existing forensic callback) + Matrix-Market dump extended to iterations
# listed in CFDX_NACA_MATRIX_DUMP_AT (tag .iter<K>). Both live inside the
# forensic gate of the pressure corrector; the iteration-1 dump behaviour is
# unchanged when CFDX_NACA_MATRIX_DUMP_AT is unset.
SOLVER_OLD = r'''                // Optional Matrix-Market snapshot of the exact matrix and
                // RHS passed to Krylov. No dump is written unless requested.
                const char* dump_prefix = std::getenv("CFDX_NACA_MATRIX_DUMP_PREFIX");
                if (dump_prefix && *dump_prefix && iter == 1 && corr == 0) {
                    const std::string dump_stem =
                        std::string(dump_prefix) + ".cells" + std::to_string(nc);
'''

SOLVER_NEW = r'''                // One-shot least-squares stencil audit (evidence only). The
                // momentum viscous terms reconstruct velocity gradients with
                // the one-ring internal-neighbour least-squares stencil that
                // compute_gradient_least_squares uses under
                // EXCLUDE_BOUNDARY. StencilQuality is a pure geometric
                // property of the 3x3 normal equations, so zero-valued
                // samples reproduce the production conditioning exactly.
                // The mesh is two-dimensional and coplanar, so rank 2 is the
                // expected stencil; the reported pivot ratio is the
                // in-plane conditioning of the pivoted normal equations.
                if (iter == 1 && corr == 0) {
                    const auto* ls_cell_faces = mesh.cells().faces_data();
                    const auto* ls_cell_offsets = mesh.cells().offsets_data();
                    std::size_t ls_rank3 = 0;
                    std::size_t ls_rank2 = 0;
                    std::size_t ls_rank1 = 0;
                    std::size_t ls_rank0 = 0;
                    std::size_t ls_samples_min =
                        std::numeric_limits<std::size_t>::max();
                    std::size_t ls_condition_max_cell = 0;
                    double ls_condition_max = 0.0;
                    std::vector<std::pair<double, std::size_t>> ls_worst;
                    for (std::size_t c = 0; c < nc; ++c) {
                        std::vector<Vec3> ls_centres;
                        std::vector<double> ls_zero_values;
                        for (Offset k = ls_cell_offsets[c];
                             k < ls_cell_offsets[c + 1]; ++k) {
                            const std::size_t f = ls_cell_faces[k];
                            const std::size_t ls_owner =
                                mesh.ownership().owner(f);
                            const std::int64_t ls_neighbour =
                                mesh.ownership().neighbour(f);
                            std::size_t ls_other = nc;
                            if (ls_owner == c) {
                                if (ls_neighbour >= 0)
                                    ls_other =
                                        static_cast<std::size_t>(ls_neighbour);
                            } else {
                                ls_other = ls_owner;
                            }
                            if (ls_other >= nc || ls_other == c) continue;
                            ls_centres.push_back(geometry.cell_centres[ls_other]);
                            ls_zero_values.push_back(0.0);
                        }
                        cfdx::core::StencilQuality ls_quality;
                        cfdx::core::least_squares_gradient(
                            geometry.cell_centres[c], 0.0, ls_centres,
                            ls_zero_values, &ls_quality);
                        if (ls_quality.rank == 3) {
                            ++ls_rank3;
                        } else if (ls_quality.rank == 2) {
                            ++ls_rank2;
                        } else if (ls_quality.rank == 1) {
                            ++ls_rank1;
                        } else {
                            ++ls_rank0;
                        }
                        ls_samples_min =
                            std::min(ls_samples_min, ls_quality.samples);
                        const double ls_condition =
                            ls_quality.min_pivot > 0.0
                                ? ls_quality.max_pivot / ls_quality.min_pivot
                                : std::numeric_limits<double>::infinity();
                        if (std::isfinite(ls_condition)) {
                            ls_worst.emplace_back(ls_condition, c);
                            if (ls_condition > ls_condition_max) {
                                ls_condition_max = ls_condition;
                                ls_condition_max_cell = c;
                            }
                        }
                    }
                    if (ls_worst.empty()) {
                        std::cerr << std::setprecision(17)
                                  << "NACA_LS_STENCIL_QUALITY"
                                  << " cells=" << nc
                                  << " rank3=" << ls_rank3
                                  << " rank2=" << ls_rank2
                                  << " rank1=" << ls_rank1
                                  << " rank0=" << ls_rank0
                                  << " samples_min="
                                  << (ls_samples_min ==
                                            std::numeric_limits<
                                                std::size_t>::max()
                                        ? 0
                                        : ls_samples_min)
                                  << " pivot_ratio_max=none"
                                  << " pivot_ratio_max_cell=none"
                                  << " pivot_ratio_max_xyz=none" << '\n';
                    } else {
                        const Vec3 ls_worst_point =
                            geometry.cell_centres[ls_condition_max_cell];
                        std::cerr << std::setprecision(17)
                                  << "NACA_LS_STENCIL_QUALITY"
                                  << " cells=" << nc
                                  << " rank3=" << ls_rank3
                                  << " rank2=" << ls_rank2
                                  << " rank1=" << ls_rank1
                                  << " rank0=" << ls_rank0
                                  << " samples_min="
                                  << (ls_samples_min ==
                                            std::numeric_limits<
                                                std::size_t>::max()
                                        ? 0
                                        : ls_samples_min)
                                  << " pivot_ratio_max=" << ls_condition_max
                                  << " pivot_ratio_max_cell="
                                  << ls_condition_max_cell
                                  << " pivot_ratio_max_xyz=" << ls_worst_point.x
                                  << ":" << ls_worst_point.y
                                  << ":" << ls_worst_point.z << '\n';
                        std::sort(ls_worst.begin(), ls_worst.end(),
                                  [](const std::pair<double, std::size_t>& a,
                                     const std::pair<double, std::size_t>& b) {
                                      return a.first > b.first;
                                  });
                        const std::size_t ls_reported =
                            std::min<std::size_t>(ls_worst.size(), 8);
                        for (std::size_t w = 0; w < ls_reported; ++w) {
                            const Vec3 ls_point =
                                geometry.cell_centres[ls_worst[w].second];
                            std::cerr << std::setprecision(17)
                                      << "NACA_LS_STENCIL_WORST"
                                      << " cell=" << ls_worst[w].second
                                      << " pivot_ratio=" << ls_worst[w].first
                                      << " xyz=" << ls_point.x << ":"
                                      << ls_point.y << ":" << ls_point.z
                                      << '\n';
                        }
                    }
                }

                // Optional Matrix-Market snapshot of the exact matrix and
                // RHS passed to Krylov. No dump is written unless requested.
                // CFDX_NACA_MATRIX_DUMP_AT lists additional nonlinear
                // iterations (comma-separated) whose first pressure corrector
                // is snapshotted with an .iter<K> tag, so the system at a
                // reproduced divergence can be compared with the first one.
                // Malformed tokens are ignored; unset restores the
                // iteration-1-only behaviour.
                const char* dump_prefix = std::getenv("CFDX_NACA_MATRIX_DUMP_PREFIX");
                const char* dump_at = std::getenv("CFDX_NACA_MATRIX_DUMP_AT");
                std::string dump_iteration_tag;
                if (dump_prefix && *dump_prefix && dump_at && *dump_at && corr == 0) {
                    std::istringstream dump_at_stream(dump_at);
                    std::string dump_at_token;
                    while (std::getline(dump_at_stream, dump_at_token, ',')) {
                        if (dump_at_token.empty()) continue;
                        try {
                            const unsigned long long parsed =
                                std::stoull(dump_at_token);
                            if (parsed == static_cast<unsigned long long>(iter)) {
                                dump_iteration_tag = std::string(".iter") +
                                    dump_at_token;
                                break;
                            }
                        } catch (const std::exception&) {
                            // Malformed evidence-only tokens are ignored.
                        }
                    }
                }
                if (dump_prefix && *dump_prefix && corr == 0 &&
                    (iter == 1 || !dump_iteration_tag.empty())) {
                    const std::string dump_stem =
                        std::string(dump_prefix) + ".cells" + std::to_string(nc)
                        + dump_iteration_tag;
'''

CMAKE_OLD_QUALIFICATION = '''        set_tests_properties(test_naca0012_qualification PROPERTIES
            ENVIRONMENT "CFDX_NACA_MATRIX_DUMP_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/naca-pressure-qualification;CFDX_NACA_AB_ALPHA_U=0.60;CFDX_NACA_AB_ALPHA_P=0.10")'''
CMAKE_NEW_QUALIFICATION = '''        set_tests_properties(test_naca0012_qualification PROPERTIES
            ENVIRONMENT "CFDX_NACA_MATRIX_DUMP_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/naca-pressure-qualification;CFDX_NACA_AB_ALPHA_U=0.60;CFDX_NACA_AB_ALPHA_P=0.10;CFDX_NACA_MATRIX_DUMP_AT=62")'''

CMAKE_OLD_WALL_SHEAR = '''        set_tests_properties(test_naca0012_wall_shear_audit PROPERTIES
            ENVIRONMENT "CFDX_NACA_MATRIX_DUMP_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/naca-pressure-wall-shear;CFDX_NACA_AB_ALPHA_U=0.60;CFDX_NACA_AB_ALPHA_P=0.10")'''
CMAKE_NEW_WALL_SHEAR = '''        set_tests_properties(test_naca0012_wall_shear_audit PROPERTIES
            ENVIRONMENT "CFDX_NACA_MATRIX_DUMP_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/naca-pressure-wall-shear;CFDX_NACA_AB_ALPHA_U=0.60;CFDX_NACA_AB_ALPHA_P=0.10;CFDX_NACA_MATRIX_DUMP_AT=62")'''

CMAKE_OLD_SPECTRUM = '''        add_test(NAME test_analyze_naca_pressure_spectrum
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/python/test_analyze_naca_pressure_spectrum.py)
        set_tests_properties(test_analyze_naca_pressure_spectrum PROPERTIES
            LABELS "qualification;validation;infrastructure")'''
CMAKE_NEW_SPECTRUM = '''        add_test(NAME test_analyze_naca_pressure_spectrum
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/python/test_analyze_naca_pressure_spectrum.py)
        set_tests_properties(test_analyze_naca_pressure_spectrum PROPERTIES
            LABELS "qualification;validation;infrastructure")
        add_test(NAME test_analyze_naca_pressure_chronology
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/python/test_analyze_naca_pressure_chronology.py)
        set_tests_properties(test_analyze_naca_pressure_chronology PROPERTIES
            LABELS "qualification;validation;infrastructure")'''


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise SystemExit(
            f"fail-closed: anchor '{label}' matched {count} times (expected 1)"
        )
    return text.replace(old, new)


def main() -> int:
    for path in (SOLVER, CMAKE, DOC, APPENDIX):
        if not path.is_file():
            raise SystemExit(f"fail-closed: missing input file {path}")

    solver = SOLVER.read_text(encoding="utf-8")
    if STENCIL_MARKER in solver:
        raise SystemExit("fail-closed: solver already carries the stencil audit")
    solver = replace_once(
        solver, SOLVER_OLD, SOLVER_NEW, "solver dump + stencil splice"
    )
    SOLVER.write_text(solver, encoding="utf-8")

    cmake = CMAKE.read_text(encoding="utf-8")
    cmake = replace_once(
        cmake,
        CMAKE_OLD_QUALIFICATION,
        CMAKE_NEW_QUALIFICATION,
        "cmake qualification environment",
    )
    cmake = replace_once(
        cmake,
        CMAKE_OLD_WALL_SHEAR,
        CMAKE_NEW_WALL_SHEAR,
        "cmake wall-shear environment",
    )
    cmake = replace_once(
        cmake, CMAKE_OLD_SPECTRUM, CMAKE_NEW_SPECTRUM, "cmake chronology test"
    )
    CMAKE.write_text(cmake, encoding="utf-8")

    doc = DOC.read_text(encoding="utf-8")
    if APPENDIX_MARKER in doc:
        raise SystemExit("fail-closed: doc already carries the v2 appendix")
    appendix = APPENDIX.read_text(encoding="utf-8")
    doc = doc.rstrip("\n") + "\n\n" + appendix.rstrip("\n") + "\n"
    DOC.write_text(doc, encoding="utf-8")

    print("applied: solver stencil audit + iter-K dump, cmake env + chronology test, doc appendix")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
