#!/usr/bin/env bash
# Reproducible, fail-closed local qualification of the supplied PR #1133 tree.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build="${CFDX_AUDIT_BUILD:-${root}/build-audit-pr1133}"
logs="${CFDX_AUDIT_LOGS:-${build}/audit-logs}"
mkdir -p "${logs}"
cmake -S "${root}" -B "${build}" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_ENABLE_MPI=OFF \
  -DCFDX_BUILD_TESTS=ON \
  -DCFDX_BUILD_RUNTIME_TESTS=OFF \
  -DCFDX_ENABLE_LONG_VALIDATION=ON \
  > "${logs}/configure.log" 2>&1
cmake --build "${build}" -j "${CFDX_AUDIT_JOBS:-4}" --target \
  test_n9_s9_rhie_chow_consistency \
  test_naca0012_qualification test_naca0012_wall_shear_audit \
  test_ghia_cavity test_wall_distance_methods test_wall_distance_turbulence \
  test_cg_solver test_amg_preconditioner_qualification \
  test_vof_curvature_contract test_phase10_turbulence_hardening \
  test_turbulence_qualification_equations test_production_turbulence_contract \
  test_thermal_vv test_radiation_vv \
  > "${logs}/build.log" 2>&1
quick='^(test_n9_s9_rhie_chow_consistency|test_ghia_cavity_quick|test_wall_distance_methods|test_wall_distance_turbulence|test_naca0012_quick|test_naca0012_wall_shear_audit_quick|test_cg_solver|test_amg_preconditioner_qualification|test_vof_curvature_contract|test_phase10_turbulence_hardening|test_turbulence_qualification_equations|test_production_turbulence_contract|test_thermal_vv|test_radiation_vv)$'
ctest --test-dir "${build}" -R "${quick}" -j 4 --output-on-failure \
  2>&1 | tee "${logs}/quick-ctest.log"
if [[ "${CFDX_AUDIT_LONG:-0}" == 1 ]]; then
  # Distinct campaigns: a failed NACA test must not suppress N9 evidence.
  for case in test_n9_physical_matrix test_naca0012_qualification test_naca0012_wall_shear_audit test_ghia_cavity_profile_study; do
    ctest --test-dir "${build}" -R "^${case}$" --output-on-failure -V \
      > "${logs}/${case}.log" 2>&1 || {
        echo "FAIL ${case}: see ${logs}/${case}.log" >&2
        exit 1
      }
  done
fi
echo "PASS selected qualification tests; logs: ${logs}"
