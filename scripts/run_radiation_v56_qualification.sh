#!/bin/sh
# Research only. Keep native GPU/MPI unqualified until real kernels/halos run.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v56"}
mode=${2:---quick}
case "$mode" in
    --quick|--extreme|--extended|--long) ;;
    --mpi) echo 'NOT_QUALIFIED: no native MPI HO hi+lo halo/reduction run' >&2; exit 3;;
    --gpu) echo 'NOT_QUALIFIED: no native GPU HO hi+lo arithmetic kernel' >&2; exit 3;;
    *) echo "Usage: $0 [build] [--quick|--extreme|--extended|--long|--mpi|--gpu]" >&2; exit 2;;
esac
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DCFDX_BUILD_TESTS=ON -DCFDX_BUILD_PYTHON=OFF \
    -DCFDX_ENABLE_MPI=OFF -DCFDX_ENABLE_GPU=OFF \
    -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 2 --target test_dom_holo_hilo_v56
ctest --test-dir "$build" --output-on-failure --parallel 1 --no-tests=error \
    -R '^(test_dom_holo_hilo_v56|test_dom_holo_hilo_v56_closure|test_dom_holo_hilo_v56_extreme_preflight)$'
if [ "$mode" = '--extreme' ] || [ "$mode" = '--extended' ] || [ "$mode" = '--long' ]; then
    ctest --test-dir "$build" --output-on-failure --parallel 1 --no-tests=error \
        -R '^test_dom_holo_hilo_v56_extreme_one$'
fi
if [ "$mode" = '--long' ]; then
    ctest --test-dir "$build" --output-on-failure --parallel 1 --no-tests=error \
        -R '^test_dom_holo_hilo_v56_extreme$'
fi
if [ "$mode" = '--extended' ]; then
    cmake --build "$build" --parallel 2 --target \
        test_dom_hilo_arithmetic_v55 test_dom_independent_angular_v54
    ctest --test-dir "$build" --output-on-failure --parallel 1 --no-tests=error \
        -R '^(test_dom_hilo_arithmetic_v55|test_dom_independent_angular_v54|test_dom_independent_angular_v54_analytic|test_dom_independent_angular_v54_stress)$'
fi
