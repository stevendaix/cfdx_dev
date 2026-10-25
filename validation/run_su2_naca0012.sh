#!/bin/bash
# Input-reader qualification only; the Python entry point records its scope.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
exec "${PYTHON:-python3}" "$ROOT/validation/run_su2_naca0012.py" "$@"
