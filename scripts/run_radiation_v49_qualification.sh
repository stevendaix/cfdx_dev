#!/bin/sh
# V49 experimental gray SMM/HOLO qualification, no production flags changed.
# Quick mode builds a bounded subset. Full V47 long stays explicit opt-in.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v49"}
mode=${2:---quick}
case "$mode" in --quick|--extended) ;; *)
  echo "Usage: $0 [build-directory] [--quick|--extended]" >&2;exit 2;;esac
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF \
  -DCFDX_BUILD_PYTHON=OFF -DCFDX_ENABLE_GPU=OFF \
  -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 3 --target \
  test_dom_smm_gray_wall_v49 test_dom_smm_second_moment_v48 \
  test_dom_energy_gate_v46 test_dom_moment_schur_v44
ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
  -R '^(test_dom_smm_gray_wall_v49|test_dom_smm_second_moment_v48|test_dom_energy_gate_v46|test_dom_moment_schur_v44)$'
if [ "$mode" = "--extended" ]; then
  cmake --build "$build" --parallel 3 --target \
    test_dom_holo_face_closure_v47 test_dom_face_moment_v45 \
    test_dom_oblique_scattering_v43 test_dom_accepted_state_v42
  ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
    -R '^(test_dom_face_moment_v45|test_dom_oblique_scattering_v43|test_dom_accepted_state_v42)$'
  "$build/test_dom_holo_face_closure_v47" --sanitizer-small
  # Full V47 is a separate long campaign. Run manually if hardware permits:
  # ctest --test-dir "$build" -R '^test_dom_holo_face_closure_v47$'
fi
