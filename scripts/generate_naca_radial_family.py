"""Radial sensitivity family, not by itself a valid uniform refinement series."""

import argparse
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "generate_2d_validation_meshes.py"
SPEC = importlib.util.spec_from_file_location(
    "cfdx_generate_2d_validation_meshes", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MESHGEN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MESHGEN)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--surface-points", type=int, default=64)
    args = parser.parse_args()
    if args.surface_points < 16:
        raise ValueError("at least 16 half-surface points required")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for nr in (24, 36, 48, 72):
        p, c, b = MESHGEN.naca_o_grid(args.surface_points, nr, 20.0)
        dest = args.output_dir / f"naca0012_n{args.surface_points}_r{nr}.h5"
        MESHGEN.write_mesh(dest, p, c, b)
        print(f"WROTE {dest} cells={len(c)}", flush=True)


if __name__ == "__main__":
    main()
