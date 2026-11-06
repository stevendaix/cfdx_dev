#!/usr/bin/env bash
# Run actual production prerequisites followed by fail-closed TH-01/03/07/08 comparison.
# Usage: bash scripts/run_thermohydraulic_qualification.sh [solver-exports-dir]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
EVIDENCE_DIR="${1:-build/thermohydraulic/solver-exports}"
BUILD_DIR="${CFDX_TH_BUILD_DIR:-build/th-qualification}"
mkdir -p build/thermohydraulic
python scripts/validate_thermohydraulic_reference.py \
  --output build/thermohydraulic/reference-contracts.json
cmake -S . -B "$BUILD_DIR" -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_RUNTIME_TESTS=OFF
cmake --build "$BUILD_DIR" --parallel 2 --target \
  cfdx_production_solver test_ghia_cavity test_thermohydraulic_observables test_thermal_driver_vv \
  test_production_buoyancy_outer_loop_vv test_natural_convection_de_vahl_davis
ctest --test-dir "$BUILD_DIR" --output-on-failure --timeout 450 \
  -R '^(test_ghia_cavity_quick|test_thermohydraulic_observables|test_thermal_driver_vv|test_production_buoyancy_outer_loop_vv|test_natural_convection_de_vahl_davis)$'
# Mandatory physical stage: no complete solver evidence -> exit code 2.
python scripts/qualify_thermohydraulic_benchmarks.py \
  --evidence-dir "$EVIDENCE_DIR" \
  --output build/thermohydraulic/physical-prequalification.json \
  --require-pass
