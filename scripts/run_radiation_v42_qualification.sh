#!/bin/sh
# V42 radiation qualification. Quick = 11 actually replayed tests.
# --extended additionally compiles/runs 4 slower V37-V40 CHT/P1 suites.
# No physics thresholds are changed and no GPU/MPI qualification is implied.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"${root}/build-radiation-v42"}
mode=${2:-"--quick"}
case "$mode" in
    --quick|--extended) ;;
    *) echo "Usage: $0 [build_directory] [--quick|--extended]" >&2; exit 2;;
esac
cmake -S "$root" -B "$build" \
    -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON \
    -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_PYTHON=OFF \
    -DCFDX_ENABLE_GPU=OFF -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 3 --target \
    test_dom_accepted_state_v42 \
    test_cht_p1_face_coupled_v41 test_dom_separate_refinement_v41 \
    test_radiation_equal_error_cpu_v41 test_p1_marshak_skew_v40 \
    test_dom_angular_refinement_v39 test_radiation_dom_cache_v38 \
    test_participating_radiation_v37 test_radiation_vv \
    test_m1_m4_physics test_thermal_radiation_model_matrix
regex='^(test_dom_accepted_state_v42|test_cht_p1_face_coupled_v41|test_dom_separate_refinement_v41|test_radiation_equal_error_cpu_v41|test_p1_marshak_skew_v40|test_dom_angular_refinement_v39|test_radiation_dom_cache_v38|test_participating_radiation_v37|test_radiation_vv|test_m1_m4_physics|test_thermal_radiation_model_matrix)'
if [ "$mode" = "--extended" ]; then
    cmake --build "$build" --parallel 3 --target \
        test_cht_p1_coupled_v40 test_p1_dom_slab_oracle_v40 \
        test_p1_marshak_v39 test_cht_nonconformal_transient_v37
    regex="${regex%?}|test_cht_p1_coupled_v40|test_p1_dom_slab_oracle_v40|test_p1_marshak_v39|test_cht_nonconformal_transient_v37)"
fi
ctest --test-dir "$build" --output-on-failure --no-tests=error -R "$regex"
