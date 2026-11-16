#!/bin/sh
# V45: face moment balance and Schur/Jacobi research qualification.
# No CUDA/MPI or production DOM code changes. Diagnostic O0 applies only
# to large validation drivers, NEVER use this mode for CPU speed claims.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v45"}
mode=${2:---quick}
case "$mode" in --quick|--extended|--spectrum) ;; *)
    echo "Usage: $0 [build-directory] [--quick|--extended|--spectrum]" >&2
    exit 2;; esac
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF \
    -DCFDX_BUILD_PYTHON=OFF -DCFDX_ENABLE_GPU=OFF \
    -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 3 --target \
    test_dom_face_moment_v45 test_dom_moment_schur_v44 \
    test_dom_oblique_scattering_v43 test_dom_accepted_state_v42
regex='^(test_dom_face_moment_v45|test_dom_moment_schur_v44|test_dom_oblique_scattering_v43|test_dom_accepted_state_v42)$'
if [ "$mode" = "--extended" ]; then
    cmake --build "$build" --parallel 3 --target \
        test_dom_dsa_study_v43 test_participating_radiation_v37 \
        test_radiation_dom_cache_v38 test_dom_angular_refinement_v39 \
        test_dom_separate_refinement_v41 test_cht_p1_face_coupled_v41
    regex='^(test_dom_face_moment_v45|test_dom_moment_schur_v44|test_dom_oblique_scattering_v43|test_dom_accepted_state_v42|test_dom_dsa_study_v43|test_participating_radiation_v37|test_radiation_dom_cache_v38|test_dom_angular_refinement_v39|test_dom_separate_refinement_v41|test_cht_p1_face_coupled_v41)$'
fi
ctest --test-dir "$build" --output-on-failure --no-tests=error -R "$regex"
if [ "$mode" = "--spectrum" ]; then
    evidence="$build/v45_spectral_evidence"
    mkdir -p "$evidence"
    "$build/test_dom_face_moment_v45" --spectrum "$evidence/V45" \
        > "$evidence/driver.log" 2>&1
    python3 "$root/scripts/analyze_dom_v45_spectrum.py" "$evidence" \
        --output "$evidence/spectral_report.json"
    echo "V45 spectral evidence: $evidence"
fi
