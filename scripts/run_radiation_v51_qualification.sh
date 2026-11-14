#!/bin/sh
# V51 full moment/wall acceleration, independent angular FV gates.
# No production solver or historical threshold is altered by this script.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v51"}
mode=${2:---quick}
case "$mode" in --quick|--stress|--extended) ;; *)
  echo "Usage: $0 [build-dir] [--quick|--stress|--extended]" >&2
  exit 2;; esac
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF \
  -DCFDX_BUILD_PYTHON=OFF -DCFDX_ENABLE_GPU=OFF \
  -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 3 --target \
  test_dom_holo_global_v51 test_dom_v51_anderson_math \
  test_dom_holo_variable_v50 \
  test_dom_smm_gray_wall_v49 test_dom_smm_second_moment_v48 \
  test_dom_energy_gate_v46
ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
  -R '^(test_dom_holo_global_v51|test_dom_v51_anderson_math|test_dom_holo_variable_v50|test_dom_smm_gray_wall_v49|test_dom_smm_second_moment_v48|test_dom_energy_gate_v46)$'
if [ "$mode" = "--stress" ] || [ "$mode" = "--extended" ]; then
  # Run separately: large angular FV tests must not compete for CPU/memory.
  ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
    -R '^test_dom_holo_v51_stress_(96|200|640)$'
fi
if [ "$mode" = "--extended" ]; then
  # Historical gray wall / angular cache regression, not a complete CFDX CI.
  cmake --build "$build" --parallel 3 --target \
    test_dom_moment_schur_v44 test_dom_oblique_scattering_v43 \
    test_dom_accepted_state_v42
  ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
    -R '^(test_dom_moment_schur_v44|test_dom_oblique_scattering_v43|test_dom_accepted_state_v42)$'
fi
