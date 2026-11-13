#!/bin/sh
# CFDX V41 targeted thermal/radiation qualifications. No automatic promotion.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"${root}/build-radiation-v41"}
cmake -S "$root" -B "$build" \
    -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON \
    -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_PYTHON=OFF \
    -DCFDX_ENABLE_GPU=OFF -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 3 --target \
    test_cht_p1_face_coupled_v41 test_dom_separate_refinement_v41 \
    test_radiation_equal_error_cpu_v41 \
    test_p1_marshak_skew_v40 test_cht_p1_coupled_v40 \
    test_p1_dom_slab_oracle_v40 test_p1_marshak_v39 \
    test_dom_angular_refinement_v39 test_radiation_dom_cache_v38 \
    test_participating_radiation_v37 test_cht_nonconformal_transient_v37 \
    test_radiation_vv test_m1_m4_physics test_thermal_radiation_model_matrix
ctest --test-dir "$build" --output-on-failure --no-tests=error \
    -R '^(test_cht_p1_face_coupled_v41|test_dom_separate_refinement_v41|test_radiation_equal_error_cpu_v41|test_p1_marshak_skew_v40|test_cht_p1_coupled_v40|test_p1_dom_slab_oracle_v40|test_p1_marshak_v39|test_dom_angular_refinement_v39|test_radiation_dom_cache_v38|test_participating_radiation_v37|test_cht_nonconformal_transient_v37|test_radiation_vv|test_m1_m4_physics|test_thermal_radiation_model_matrix)$'
