#!/bin/sh
# CFDX V44 - exact angular moment Schur research qualification.
# --quick: V44 + V43 black-wall DSA + V42 true FV transport defect.
# --extended: add the existing V37-V43 radiation/CHT regression selection.
# The diagnostic fast compile switch affects large TEST drivers only, NOT
# the production solver binary or the cfdx_core optimization settings.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"${root}/build-radiation-v44"}
mode=${2:-"--quick"}
case "$mode" in
    --quick|--extended) ;;
    *) echo "Usage: $0 [build_directory] [--quick|--extended]" >&2; exit 2 ;;
esac
cmake -S "$root" -B "$build" \
    -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON \
    -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_PYTHON=OFF \
    -DCFDX_ENABLE_GPU=OFF -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 3 --target \
    test_dom_moment_schur_v44 test_dom_dsa_study_v43 test_dom_accepted_state_v42
regex='^(test_dom_moment_schur_v44|test_dom_dsa_study_v43|test_dom_accepted_state_v42)$'
if [ "$mode" = "--extended" ]; then
    cmake --build "$build" --parallel 3 --target \
        test_dom_oblique_scattering_v43 test_cht_p1_face_coupled_v41 \
        test_dom_separate_refinement_v41 test_radiation_equal_error_cpu_v41 \
        test_p1_marshak_skew_v40 test_dom_angular_refinement_v39 \
        test_radiation_dom_cache_v38 test_participating_radiation_v37 \
        test_radiation_vv test_m1_m4_physics test_thermal_radiation_model_matrix
    regex='^(test_dom_moment_schur_v44|test_dom_dsa_study_v43|test_dom_accepted_state_v42|test_dom_oblique_scattering_v43|test_cht_p1_face_coupled_v41|test_dom_separate_refinement_v41|test_radiation_equal_error_cpu_v41|test_p1_marshak_skew_v40|test_dom_angular_refinement_v39|test_radiation_dom_cache_v38|test_participating_radiation_v37|test_radiation_vv|test_m1_m4_physics|test_thermal_radiation_model_matrix)$'
fi
ctest --test-dir "$build" --output-on-failure --no-tests=error -R "$regex"
