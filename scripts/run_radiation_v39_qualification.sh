#!/bin/sh
# Fast functional V39 radiation qualification; not a full industrial CFD V&V.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"${root}/build-radiation-v39"}
cmake -S "$root" -B "$build" \
    -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON \
    -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_PYTHON=OFF \
    -DCFDX_ENABLE_GPU=OFF -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 3 --target \
    test_p1_marshak_v39 test_dom_angular_refinement_v39 \
    test_radiation_dom_cache_v38 test_participating_radiation_v37 \
    test_radiation_vv test_m1_m4_physics test_thermal_radiation_model_matrix
ctest --test-dir "$build" --output-on-failure --no-tests=error \
    -R '^(test_p1_marshak_v39|test_dom_angular_refinement_v39|test_radiation_dom_cache_v38|test_participating_radiation_v37|test_radiation_vv|test_m1_m4_physics|test_thermal_radiation_model_matrix)$'
