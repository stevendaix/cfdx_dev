#!/bin/sh
# V50 opt-in HO/LO variable-optics qualification. Neither long physical
# campaigns nor production solver flags are changed by this helper.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v50"}
mode=${2:---quick}
case "$mode" in --quick|--stress|--extended) ;; *)
    echo "Usage: $0 [build-directory] [--quick|--stress|--extended]" >&2
    exit 2;; esac
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF \
    -DCFDX_BUILD_PYTHON=OFF -DCFDX_ENABLE_GPU=OFF \
    -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 3 --target \
    test_dom_holo_variable_v50 test_dom_smm_gray_wall_v49 \
    test_dom_smm_second_moment_v48 test_dom_energy_gate_v46
ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
    -R '^(test_dom_holo_variable_v50|test_dom_smm_gray_wall_v49|test_dom_smm_second_moment_v48|test_dom_energy_gate_v46)$'
if [ "$mode" = "--stress" ] || [ "$mode" = "--extended" ]; then
    ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
        -R '^test_dom_holo_v50_stress_(640|10000)$'
fi
if [ "$mode" = "--extended" ]; then
    # Long V44/V47 tests are deliberately separate and must be built/run
    # explicitly on adequate hardware. No automatic success credit.
    echo 'V44/V47 full campaigns and GPU/MPI/production CI remain separate.'
fi
