#!/bin/sh
# CFDX V57: opt-in CPU research qualification, strict physical and algebraic
# thresholds. The actual native GPU and MPI tests have NOT been carried out.
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build=${1:-"$root/build-radiation-v57"}
mode=${2:---quick}
case "$mode" in
    --quick|--extended|--long) ;;
    --mpi) echo 'NOT_QUALIFIED: no native MPI HO hi+lo kernel, communication or reductions' >&2; exit 3 ;;
    --gpu) echo 'NOT_QUALIFIED: no native CUDA/HIP/SYCL HO hi+lo kernel' >&2; exit 3 ;;
    *) echo "Usage: $0 [build] [--quick|--extended|--long|--mpi|--gpu]" >&2;exit 2 ;;
esac
long=OFF
[ "$mode" = '--long' ] && long=ON
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_BUILD_PYTHON=OFF \
  -DCFDX_ENABLE_MPI=OFF -DCFDX_ENABLE_GPU=OFF \
  -DCFDX_ENABLE_LONG_VALIDATION="$long" \
  -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON \
  -DCFDX_ENABLE_RADIATION_RESEARCH_V57=ON
cmake --build "$build" --parallel 2 --target \
  test_dom_angular_hilo_v57 test_dom_holo_hilo_v57
ctest --test-dir "$build" --output-on-failure --parallel 1 --no-tests=error \
  -R '^test_dom_(angular_hilo_v57(_analytic|_extreme|_family)?|holo_hilo_v57(_extreme_one|_family_one)?)$'
if [ "$mode" = '--extended' ] || [ "$mode" = '--long' ]; then
    cmake --build "$build" --parallel 2 --target \
      test_dom_hilo_arithmetic_v55 test_dom_holo_hilo_v56 \
      test_dom_independent_angular_v54 test_dom_holo_amg_v53
    ctest --test-dir "$build" --output-on-failure --parallel 1 --no-tests=error \
      -R '^(test_dom_hilo_arithmetic_v55|test_dom_holo_hilo_v56|test_dom_independent_angular_v54(_analytic|_stress)?|test_dom_holo_amg_v53)$'
fi
if [ "$mode" = '--long' ];then
    # Full physically converged 640-cell campaign. Never replace with preflight.
    ctest --test-dir "$build" --output-on-failure --parallel 1 --no-tests=error \
      -R '^test_dom_holo_hilo_v57_(extreme|family)$'
fi
