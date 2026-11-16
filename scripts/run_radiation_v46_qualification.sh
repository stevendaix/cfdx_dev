#!/bin/sh
# V46: physically accepted DOM-Schur; extended mode includes historic gates.
# The optional O0 is TEST-ONLY and must not be used for performance claims.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v46"}
mode=${2:---quick}
case "$mode" in --quick|--extended) ;; *)
    echo "Usage: $0 [build-directory] [--quick|--extended]" >&2
    exit 2;; esac
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF \
    -DCFDX_BUILD_PYTHON=OFF -DCFDX_ENABLE_GPU=OFF \
    -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 3 --target \
    test_dom_energy_gate_v46 test_dom_face_moment_v45 \
    test_dom_moment_schur_v44 test_dom_accepted_state_v42
regex='^(test_dom_energy_gate_v46|test_dom_face_moment_v45|test_dom_moment_schur_v44|test_dom_accepted_state_v42)$'
if [ "$mode" = "--extended" ]; then
    cmake --build "$build" --parallel 3 --target \
        test_dom_oblique_scattering_v43 test_dom_dsa_study_v43 \
        test_participating_radiation_v37 test_radiation_dom_cache_v38 \
        test_dom_angular_refinement_v39 test_dom_separate_refinement_v41
    regex='^(test_dom_energy_gate_v46|test_dom_face_moment_v45|test_dom_moment_schur_v44|test_dom_accepted_state_v42|test_dom_oblique_scattering_v43|test_dom_dsa_study_v43|test_participating_radiation_v37|test_radiation_dom_cache_v38|test_dom_angular_refinement_v39|test_dom_separate_refinement_v41)$'
fi
ctest --test-dir "$build" --parallel 2 --output-on-failure --no-tests=error -R "$regex"
