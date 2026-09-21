#!/usr/bin/env python3
"""CFDX command-line entry point.

The CLI is intentionally thin: it validates/inspects native case files and
delegates numerical execution to the selected solver executable. No physics
implementation is duplicated in Python.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


REQUIRED_DATASETS = (
    "points",
    "face_vertices",
    "face_offsets",
    "owner",
    "neighbour",
    "cell_faces",
    "cell_offsets",
)


def _require_h5py():
    try:
        import h5py  # type: ignore
    except ImportError as exc:
        raise RuntimeError("h5py is required for CFDX HDF5 commands") from exc
    return h5py


def inspect_case(path: Path) -> dict:
    h5py = _require_h5py()
    with h5py.File(path, "r") as h5:
        datasets = {}
        def visitor(name, obj):
            if hasattr(obj, "shape"):
                datasets[name] = {
                    "shape": list(obj.shape),
                    "dtype": str(obj.dtype),
                    "bytes": int(obj.size * obj.dtype.itemsize),
                }
        h5.visititems(visitor)
        attrs = {str(k): v.item() if hasattr(v, "item") else v
                 for k, v in h5.attrs.items()}
        return {"file": str(path), "attributes": attrs, "datasets": datasets}


def check_case(path: Path) -> int:
    h5py = _require_h5py()
    errors = []
    with h5py.File(path, "r") as h5:
        missing = [name for name in REQUIRED_DATASETS if name not in h5]
        if missing:
            errors.append("missing required topology datasets: " + ", ".join(missing))
        if "points" in h5 and (len(h5["points"].shape) != 2 or h5["points"].shape[1] != 3):
            errors.append("points must have shape [nPoints, 3]")
        if "face_offsets" in h5 and "face_vertices" in h5:
            off = h5["face_offsets"][...]
            if len(off) == 0 or int(off[0]) != 0 or int(off[-1]) != h5["face_vertices"].size:
                errors.append("face CSR offsets are inconsistent")
        if "cell_offsets" in h5 and "cell_faces" in h5:
            off = h5["cell_offsets"][...]
            if len(off) == 0 or int(off[0]) != 0 or int(off[-1]) != h5["cell_faces"].size:
                errors.append("cell CSR offsets are inconsistent")
        schema = h5.attrs.get("schema_version")
        if schema is not None and int(schema) != 1:
            errors.append(f"unsupported schema_version={schema!r}")
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(f"CFDX CHECK: PASS {path}")
    return 0


def cmd_info(args) -> int:
    data = inspect_case(args.case)
    print(json.dumps({
        "file": data["file"],
        "attributes": data["attributes"],
        "dataset_count": len(data["datasets"]),
    }, indent=2, default=str))
    return 0


def cmd_inspect(args) -> int:
    print(json.dumps(inspect_case(args.case), indent=2, default=str))
    return 0


def cmd_check(args) -> int:
    return check_case(args.case)


def cmd_run(args) -> int:
    command = [args.solver, str(args.case), *args.solver_args]
    return subprocess.call(command)


def cmd_convert(args) -> int:
    script = Path(__file__).with_name("meshio_import.py")
    command = [sys.executable, str(script), str(args.input), str(args.output)]
    return subprocess.call(command)


def cmd_benchmark(args) -> int:
    if args.command:
        start = time.perf_counter()
        completed = subprocess.run(args.command, shell=True, check=False)
        elapsed = time.perf_counter() - start
        print(json.dumps({
            "command": args.command,
            "wall_time_s": elapsed,
            "return_code": completed.returncode,
        }, indent=2))
        return completed.returncode

    data = inspect_case(args.case)
    total_bytes = sum(item["bytes"] for item in data["datasets"].values())
    print(json.dumps({
        "case": str(args.case),
        "dataset_count": len(data["datasets"]),
        "stored_dataset_bytes": total_bytes,
    }, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cfdx")
    sub = parser.add_subparsers(dest="command", required=True)

    for name, func in (("check", cmd_check), ("info", cmd_info),
                       ("inspect", cmd_inspect)):
        p = sub.add_parser(name)
        p.add_argument("case", type=Path)
        p.set_defaults(func=func)

    p = sub.add_parser("run")
    p.add_argument("case", type=Path)
    p.add_argument("--solver", required=True)
    p.add_argument("solver_args", nargs=argparse.REMAINDER)
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("convert")
    p.add_argument("input", type=Path)
    p.add_argument("output", type=Path)
    p.set_defaults(func=cmd_convert)

    p = sub.add_parser("benchmark")
    p.add_argument("case", type=Path)
    p.add_argument("--command")
    p.set_defaults(func=cmd_benchmark)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        return int(args.func(args))
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"CFDX ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
