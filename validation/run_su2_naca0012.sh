#!/bin/bash
# SU2 NACA0012 validation: convert -> solve -> compare with theory
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="$ROOT/tests/data/su2"
TMP_DIR="$(mktemp -d /tmp/cfdx_naca0012_XXXXXX)"
CLI="$ROOT/build/cfdx_convert"
SOLVER="$ROOT/build/cfdx_production_solver"

echo "=== SU2 NACA0012 Validation ==="
echo "Case: 2-D NACA0012, Euler, Mach 0.8, AoA 1.25 deg"
echo "Mesh: 5233 vertices, 10216 triangles"
echo ""

# 1. Convert SU2 -> CFDX
echo "[1/4] Converting SU2 -> CFDX..."
cp "$DATA_DIR/mesh_NACA0012_inv.su2" "$TMP_DIR/mesh.su2"
cp "$DATA_DIR/inv_NACA0012_basic.cfg" "$TMP_DIR/config.cfg"

"$CLI" convert "$TMP_DIR" "$TMP_DIR/case.cfdx.h5"
echo "  OK Conversion complete"

# 2. Run CFDX solver
echo "[2/4] Running CFDX solver..."
cd "$TMP_DIR"
"$SOLVER" --case case.cfdx.h5 --output solution.cfdx.h5 --max-iter 1000 2>&1 | tail -5
echo "  OK Solver complete"

# 3. Extract and compare results
echo "[3/4] Extracting results..."
python3 -c "
import h5py, numpy as np
with h5py.File('$TMP_DIR/solution.cfdx.h5', 'r') as f:
    p = f['fields/scalar/pressure'][:]
    U = f['fields/vector/velocity'][:]
    print(f'Pressure range: {p.min():.2f} - {p.max():.2f} Pa')
    print(f'Velocity range: {np.linalg.norm(U, axis=1).min():.1f} - {np.linalg.norm(U, axis=1).max():.1f} m/s')
    p0 = 101325.0; M = 0.8; g = 1.4
    p_theory = p0 * (1 + 0.2 * M**2)**(-3.5)
    print(f'Freestream pressure: {p0} Pa')
    print(f'Theory (isentropic): {p_theory:.2f} Pa')
    print(f'Computed min pressure: {p.min():.2f} Pa')
    print(f'Ratio computed/theory: {p.min()/p_theory:.4f}')
"

# 4. Summary
echo ""
echo "=== Validation Summary ==="
echo "  Case:        SU2 NACA0012 Euler"
echo "  Converter:   cfdx_convert (Python bridge)"
echo "  Solver:      cfdx_production_solver"
echo "  Output:      $TMP_DIR/solution.cfdx.h5"
echo ""
echo "To compare with SU2 results, run:"
echo "  cd $TMP_DIR && su2 config.cfg"
echo ""
echo "TMP_DIR=$TMP_DIR"
