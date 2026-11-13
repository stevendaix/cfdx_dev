#!/bin/sh
# Research-only V55: do not substitute any test for native GPU or MPI.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v55"}
mode=${2:---quick}
case "$mode" in
    --quick|--family|--extended) ;;
    --mpi)
        echo 'NOT_QUALIFIED: no native MPI rank/halo validation of V55 WLS.' >&2
        exit 3;;
    --gpu)
        echo 'NOT_QUALIFIED: no native CUDA/HIP/SYCL V55 arithmetic kernel.' >&2
        exit 3;;
    *) echo "Usage: $0 [build-directory] [--quick|--family|--extended|--gpu|--mpi]" >&2
       exit 2;;
esac
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DCFDX_BUILD_TESTS=ON -DCFDX_BUILD_PYTHON=OFF \
    -DCFDX_ENABLE_MPI=OFF -DCFDX_ENABLE_GPU=OFF \
    -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 2 --target \
    test_dom_hilo_arithmetic_v55 test_dom_wls_hilo_v55
ctest --test-dir "$build" --output-on-failure --parallel 1 --no-tests=error \
    -R '^(test_dom_hilo_arithmetic_v55|test_dom_wls_hilo_v55)$'
if [ "$mode" = '--family' ] || [ "$mode" = '--extended' ]; then
    ctest --test-dir "$build" --output-on-failure --parallel 1 --no-tests=error \
        -R '^(test_dom_wls_hilo_v55_family|test_dom_wls_hilo_v55_iterative_family)$'
fi
if [ "$mode" = '--extended' ]; then
    ctest --test-dir "$build" --output-on-failure --parallel 1 --no-tests=error \
        -R '^test_dom_wls_hilo_v55_extreme$'
    # Prior V54 independent angular + staged CHT/DOM regressions;
    # do not silently omit these targets or substitute the HO/LO identity.
    cmake --build "$build" --parallel 2 --target \
        test_dom_independent_angular_v54 test_cht_dom_staged_v54 \
        test_dom_stable_wls_v54
    ctest --test-dir "$build" --output-on-failure --parallel 1 --no-tests=error \
        -R '^(test_dom_independent_angular_v54|test_dom_independent_angular_v54_analytic|test_dom_independent_angular_v54_stress|test_cht_dom_staged_v54|test_cht_dom_staged_v54_long|test_dom_stable_wls_v54)$'
fi
