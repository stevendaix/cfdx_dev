#!/bin/sh
# V54 research-only qualification: independent ordinate/CHT/DOM WLS tests.
# A "gpu" or "mpi" request is NEVER silently replaced with CPU simulation.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v54"}
mode=${2:---quick}
case "$mode" in
  --quick|--stress|--extended) ;;
  --gpu)
    echo 'NOT_QUALIFIED: V54 has no compiled/executed GPU radiation backend or device run.' >&2
    exit 3;;
  --mpi)
    echo 'NOT_QUALIFIED: V54 has no executed MPI halo-exchange or rank conservation run.' >&2
    exit 3;;
  *) echo "Usage: $0 [build-directory] [--quick|--stress|--extended|--gpu|--mpi]" >&2
     exit 2;;
esac
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_BUILD_PYTHON=OFF \
  -DCFDX_ENABLE_MPI=OFF -DCFDX_ENABLE_GPU=OFF \
  -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 2 --target \
  test_dom_independent_angular_v54 test_cht_dom_staged_v54 test_dom_stable_wls_v54
ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
  -R '^(test_dom_independent_angular_v54|test_dom_independent_angular_v54_analytic|test_cht_dom_staged_v54|test_dom_stable_wls_v54)$'
if [ "$mode" = '--stress' ] || [ "$mode" = '--extended' ]; then
  ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
    -R '^(test_dom_independent_angular_v54_stress|test_cht_dom_staged_v54_long)$'
  cmake --build "$build" --parallel 2 --target \
    test_dom_holo_amg_v53 test_dom_smm_extreme_v53
  ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
    -R '^(test_dom_holo_amg_v53|test_dom_smm_extreme_v53)$'
fi
if [ "$mode" = '--extended' ]; then
  cmake --build "$build" --parallel 2 --target \
    test_dom_layered_real_transport_v52 test_cht_p1_face_coupled_v41 \
    test_dom_smm_tetrahedral_v52
  ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
    -R '^(test_dom_layered_real_transport_v52|test_cht_p1_face_coupled_v41|test_dom_smm_tetrahedral_v52)$'
fi
