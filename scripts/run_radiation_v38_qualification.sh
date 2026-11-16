#!/bin/sh
# Targeted radiation verification on the selected V38 source tree.
# Functional tests only; the optional wall benchmark is NOT a whole-solver speed claim.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"${root}/build-radiation-v38"}
cmake -S "$root" -B "$build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCFDX_BUILD_TESTS=ON \
    -DCFDX_ENABLE_MPI=OFF \
    -DCFDX_BUILD_PYTHON=OFF \
    -DCFDX_ENABLE_GPU=OFF \
    -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 4 --target \
    test_radiation_dom_cache_v38 \
    test_participating_radiation_v37 \
    test_m1_m4_physics \
    test_radiation_vv \
    test_thermal_radiation_model_matrix
ctest --test-dir "$build" --output-on-failure \
    -R '^(test_radiation_dom_cache_v38|test_participating_radiation_v37|test_m1_m4_physics|test_radiation_vv|test_thermal_radiation_model_matrix)$'
if [ "${CFDX_V38_RUN_WALL_BENCHMARK:-0}" = 1 ]; then
    cmake --build "$build" --parallel 4 --target cfdx_benchmark_dom_cache_v38
    "$build/cfdx_benchmark_dom_cache_v38"
fi
