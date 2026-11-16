#!/bin/sh
# CFDX V47: research HOLO on the actual DOM finite-volume operator.
# Diagnostic -O0 applies ONLY to header-heavy tests, not performance claims.
# Run the heavy HOLO and legacy DSA groups separately, and serialize the
# legacy group: some simultaneous oblique/DSA solves contend severely.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v47"}
mode=${2:---quick}
case "$mode" in --quick|--extended) ;; *)
    echo "Usage: $0 [build-directory] [--quick|--extended]" >&2
    exit 2;; esac
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF \
    -DCFDX_BUILD_PYTHON=OFF -DCFDX_ENABLE_GPU=OFF \
    -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 3 --target \
    test_dom_holo_face_closure_v47 test_dom_energy_gate_v46 \
    test_dom_face_moment_v45 test_dom_moment_schur_v44 \
    test_dom_accepted_state_v42
quick='^(test_dom_holo_face_closure_v47|test_dom_energy_gate_v46|test_dom_face_moment_v45|test_dom_moment_schur_v44|test_dom_accepted_state_v42)$'
ctest --test-dir "$build" --parallel 2 --output-on-failure --no-tests=error -R "$quick"
if [ "$mode" = "--extended" ]; then
    cmake --build "$build" --parallel 3 --target \
        test_dom_oblique_scattering_v43 test_dom_dsa_study_v43 \
        test_participating_radiation_v37 test_radiation_dom_cache_v38 \
        test_dom_angular_refinement_v39 test_dom_separate_refinement_v41
    historic='^(test_dom_oblique_scattering_v43|test_dom_dsa_study_v43|test_participating_radiation_v37|test_radiation_dom_cache_v38|test_dom_angular_refinement_v39|test_dom_separate_refinement_v41)$'
    ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error -R "$historic"
fi
