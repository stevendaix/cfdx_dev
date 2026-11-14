#!/bin/sh
# V53 research-only: true HO+LO double-gate / directional AMG qualification.
# Heavy tests run serially; never reinterpret an interrupted build as a pass.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v53"}
mode=${2:---quick}
case "$mode" in --quick|--stress|--extended) ;; *)
  echo "usage: $0 [build-dir] [--quick|--stress|--extended]" >&2;exit 2;;esac
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_BUILD_PYTHON=OFF -DCFDX_ENABLE_MPI=OFF \
  -DCFDX_ENABLE_GPU=OFF -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON
cmake --build "$build" --parallel 2 --target \
  test_dom_holo_amg_v53 test_dom_smm_extreme_v53
ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
  -R '^(test_dom_holo_amg_v53|test_dom_smm_extreme_v53)$'
if [ "$mode" = "--stress" ] || [ "$mode" = "--extended" ]; then
  ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
    -R '^test_dom_holo_amg_v53_stress$'
fi
if [ "$mode" = "--extended" ]; then
  cmake --build "$build" --parallel 2 --target \
    test_dom_smm_anisotropic_v52 test_dom_smm_tetrahedral_v52 \
    test_dom_layered_oracle_v52
  ctest --test-dir "$build" --parallel 1 --output-on-failure --no-tests=error \
    -R '^(test_dom_smm_anisotropic_v52|test_dom_smm_tetrahedral_v52|test_dom_layered_oracle_v52)$'
fi
