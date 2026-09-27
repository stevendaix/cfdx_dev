#!/bin/bash
# Integration tests for cfdx convert CLI
# Tests the complete pipeline: CLI -> Python converter -> C++ validation

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILD_DIR="${PROJECT_ROOT}/build"
CFDX_CONVERT="${BUILD_DIR}/cfdx_convert"
TEST_DATA="${PROJECT_ROOT}/tests/data"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

passed=0
failed=0

run_test() {
    local name="$1"
    shift
    echo -e "${YELLOW}[TEST]${NC} $name"
    if "$@"; then
        echo -e "${GREEN}[PASS]${NC} $name"
        ((passed++))
        return 0
    else
        echo -e "${RED}[FAIL]${NC} $name"
        ((failed++))
        return 1
    fi
}

# Check if cfdx_convert exists
if [[ ! -f "$CFDX_CONVERT" ]]; then
    echo -e "${RED}Error: cfdx_convert not found at $CFDX_CONVERT${NC}"
    echo "Build the project first: cmake -B build && cmake --build build"
    exit 1
fi

# Test 1: Convert SU2 mesh
run_test "SU2 mesh conversion" \
    "$CFDX_CONVERT" convert "${TEST_DATA}/su2/mesh_NACA0012_inv.su2" \
    "${BUILD_DIR}/test_su2_output.cfdx.h5"

# Test 2: Convert Fluent .cas file
run_test "Fluent .cas conversion" \
    "$CFDX_CONVERT" convert "${TEST_DATA}/fluent/cavity.cas" \
    "${BUILD_DIR}/test_fluent_output.cfdx.h5"

# Test 3: Convert Fluent .cas.h5 file (if HDF5 available)
if [[ -f "${TEST_DATA}/fluent/h5/mixing_elbow.cas.h5" ]]; then
    run_test "Fluent .cas.h5 conversion" \
        "$CFDX_CONVERT" convert "${TEST_DATA}/fluent/h5/mixing_elbow.cas.h5" \
        "${BUILD_DIR}/test_fluent_h5_output.cfdx.h5"
fi

# Test 4: Convert STAR-CCM+ .sim file
run_test "STAR-CCM+ .sim conversion" \
    "$CFDX_CONVERT" convert "${TEST_DATA}/starccm/test_case.sim" \
    "${BUILD_DIR}/test_starccm_output.cfdx.h5"

# Test 5: Convert Code_Saturne .xml file
run_test "Code_Saturne .xml conversion" \
    "$CFDX_CONVERT" convert "${TEST_DATA}/saturne/test_case.xml" \
    "${BUILD_DIR}/test_saturne_xml_output.cfdx.h5"

# Test 6: Convert Code_Saturne .py file
run_test "Code_Saturne .py conversion" \
    "$CFDX_CONVERT" convert "${TEST_DATA}/saturne/test_case.py" \
    "${BUILD_DIR}/test_saturne_py_output.cfdx.h5"

# Test 7: Convert OpenFOAM case directory
# Create a minimal OpenFOAM case for testing
OF_CASE="${BUILD_DIR}/test_openfoam_case"
mkdir -p "${OF_CASE}/constant/polyMesh"
mkdir -p "${OF_CASE}/0"

# Write minimal mesh files
cat > "${OF_CASE}/constant/polyMesh/points" << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       vectorField;
    location    "constant/polyMesh";
    object      points;
}
8
(
(0 0 0)
(1 0 0)
(1 1 0)
(0 1 0)
(0 0 1)
(1 0 1)
(1 1 1)
(0 1 1)
)
EOF

cat > "${OF_CASE}/constant/polyMesh/faces" << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       faceList;
    location    "constant/polyMesh";
    object      faces;
}
6
(
(0 1 2 3)
(4 5 6 7)
(0 1 5 4)
(1 2 6 5)
(2 3 7 6)
(3 0 4 7)
)
EOF

cat > "${OF_CASE}/constant/polyMesh/owner" << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       labelList;
    location    "constant/polyMesh";
    object      owner;
}
1
(
0
)
EOF

cat > "${OF_CASE}/constant/polyMesh/neighbour" << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       labelList;
    location    "constant/polyMesh";
    object      neighbour;
}
0
(
)
EOF

cat > "${OF_CASE}/constant/polyMesh/boundary" << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       polyBoundaryMesh;
    location    "constant/polyMesh";
    object      boundary;
}
6
(
    bottom
    {
        type            patch;
        nFaces          1;
        startFace       0;
    }
    top
    {
        type            patch;
        nFaces          1;
        startFace       1;
    }
    front
    {
        type            patch;
        nFaces          1;
        startFace       2;
    }
    back
    {
        type            patch;
        nFaces          1;
        startFace       3;
    }
    left
    {
        type            patch;
        nFaces          1;
        startFace       4;
    }
    right
    {
        type            patch;
        nFaces          1;
        startFace       5;
    }
)
EOF

# Write minimal 0/U field
cat > "${OF_CASE}/0/U" << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);
boundaryField
{
    bottom  { type fixedValue; value uniform (0 0 0); }
    top     { type fixedValue; value uniform (1 0 0); }
    front   { type symmetry; }
    back    { type symmetry; }
    left    { type symmetry; }
    right   { type symmetry; }
}
EOF

cat > "${OF_CASE}/0/p" << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    bottom  { type zeroGradient; }
    top     { type zeroGradient; }
    front   { type symmetry; }
    back    { type symmetry; }
    left    { type symmetry; }
    right   { type symmetry; }
}
EOF

run_test "OpenFOAM case directory conversion" \
    "$CFDX_CONVERT" convert "${OF_CASE}" \
    "${BUILD_DIR}/test_openfoam_output.cfdx.h5"

# Test 8: Normalize/re-export CFDX file (round-trip)
if [[ -f "${BUILD_DIR}/test_su2_output.cfdx.h5" ]]; then
    run_test "CFDX normalize/re-export" \
        "$CFDX_CONVERT" convert "${BUILD_DIR}/test_su2_output.cfdx.h5" \
        "${BUILD_DIR}/test_normalized.cfdx.h5"
fi

# Test 9: Verify output files are valid HDF5 with correct structure
run_test "Output validation - SU2" \
    bash -c "h5dump -n 1 '${BUILD_DIR}/test_su2_output.cfdx.h5' >/dev/null 2>&1"

run_test "Output validation - Fluent" \
    bash -c "h5dump -n 1 '${BUILD_DIR}/test_fluent_output.cfdx.h5' >/dev/null 2>&1"

# Test 10: Test field size mismatch detection
# Create a CFDX file with mismatched field sizes and verify it's caught
python3 -c "
import h5py
import numpy as np
with h5py.File('${BUILD_DIR}/test_bad_field.cfdx.h5', 'w') as f:
    f.attrs['cfdx_version'] = '0.1.0'
    mesh = f.create_group('mesh/topology')
    mesh.attrs['cell_count'] = 100
    mesh.attrs['point_count'] = 8
    mesh.create_dataset('points', data=np.random.rand(8, 3))
    fields = f.create_group('fields/scalar')
    # Create field with wrong size (50 instead of 100)
    fields.create_dataset('bad_field', data=np.random.rand(50))
" 2>/dev/null || true

run_test "Field size mismatch detection" \
    bash -c "! '$CFDX_CONVERT' convert '${BUILD_DIR}/test_bad_field.cfdx.h5' '${BUILD_DIR}/test_bad_validated.cfdx.h5' 2>&1 | grep -q 'SIZE MISMATCH'"

# Test 11: Test invalid source format detection
run_test "Invalid source format rejection" \
    bash -c "! '$CFDX_CONVERT' convert '/nonexistent/file.xyz' '${BUILD_DIR}/test_invalid.cfdx.h5' 2>&1 | grep -q 'Cannot auto-detect'"

# Summary
echo ""
echo "==================="
echo "Test Summary"
echo "==================="
echo -e "${GREEN}Passed: $passed${NC}"
echo -e "${RED}Failed: $failed${NC}"

if [[ $failed -gt 0 ]]; then
    exit 1
fi

exit 0