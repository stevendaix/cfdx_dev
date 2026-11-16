#!/bin/sh
# Targeted V48 discrete second-moment qualification; no production switches.
# Run legacy full HOLO separately; its large anisotropic cases can take longer.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v48"}
mode=${2:---quick}
case "$mode" in --quick|--extended) ;; *)
  echo "Usage: $0 [build-directory] [--quick|--extended]" >&2; exit 2;; esac
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF \
  -DCFDX_BUILD_PYTHON=OFF -DCFDX_ENABLE_GPU=OFF \
  -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 3 --target \
  test_dom_smm_second_moment_v48 test_dom_holo_face_closure_v47 \
  test_dom_energy_gate_v46 test_dom_face_moment_v45 test_dom_moment_schur_v44
quick='^(test_dom_smm_second_moment_v48|test_dom_energy_gate_v46|test_dom_face_moment_v45|test_dom_moment_schur_v44)$'
ctest --test-dir "$build" --parallel 2 --output-on-failure --no-tests=error -R "$quick"
# Independent small V47 regression with same original acceptance criteria.
"$build/test_dom_holo_face_closure_v47" --sanitizer-small
if [ "$mode" = "--extended" ]; then
  cmake --build "$build" --parallel 3 --target \
    test_dom_oblique_scattering_v43 test_dom_accepted_state_v42 \
    test_participating_radiation_v37
  historic='^(test_dom_oblique_scattering_v43|test_dom_accepted_state_v42|test_participating_radiation_v37)$'
  ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error -R "$historic"
  # Full V47 long is deliberately opt-in and runs alone.
  ctest --test-dir "$build" --parallel 1 --output-on-failure \
    --no-tests=error -R '^test_dom_holo_face_closure_v47$'
fi
