#!/bin/sh
# V52 standalone research qualification, serial heavy transport.
# AMG must pass both strict algebraic and independent energy gates.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v52"}
mode=${2:---quick}
case "$mode" in --quick|--stress|--extended) ;; *)
  echo "Usage: $0 [build-dir] [--quick|--stress|--extended]" >&2;exit 2;;esac
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_PYTHON=OFF \
  -DCFDX_ENABLE_GPU=OFF -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 2 --target \
  test_dom_smm_tetrahedral_v52 test_dom_smm_anisotropic_v52 \
  test_dom_layered_oracle_v52 test_dom_layered_real_transport_v52 \
  test_dom_holo_refresh_640_v52
ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
  -R '^(test_dom_smm_tetrahedral_v52|test_dom_smm_anisotropic_v52|test_dom_layered_oracle_v52|test_dom_layered_real_transport_v52|test_dom_holo_refresh_640_v52)$'
if [ "$mode" = "--stress" ] || [ "$mode" = "--extended" ]; then
  ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
    -R '^test_dom_smm_anisotropic_v52_stress$'
fi
if [ "$mode" = "--extended" ]; then
  cmake --build "$build" --parallel 2 --target \
    test_dom_holo_global_v51 test_dom_smm_gray_wall_v49 \
    test_dom_smm_second_moment_v48
  ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
    -R '^(test_dom_holo_global_v51|test_dom_smm_gray_wall_v49|test_dom_smm_second_moment_v48)$'
fi
