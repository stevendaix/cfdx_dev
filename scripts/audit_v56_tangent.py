"""Fail-closed evidence gate for a local convection tangent (NOT global JFNK)."""
from __future__ import annotations
import argparse
import json
import math
import pathlib
import re

FIELDS = ("schemes", "faces", "verified", "nonsmooth", "max_primal_error",
          "max_tangent_error", "max_cell_parity", "max_cell_tangent_error",
          "max_additivity", "kink_face_id", "kink_detected")


def audit(log: str) -> dict:
    lines = re.findall(r"^V56_BFS152_FACE_AUDIT: ([^\n]+)$", log, re.M)
    if len(lines) != 1:
        raise ValueError("V56 evidence must contain exactly one actual face audit")
    tokens = dict(s.split("=", 1) for s in lines[0].split() if "=" in s)
    if set(tokens) != set(FIELDS):
        raise ValueError("V56 incomplete or unexpected metric")
    data = {k: float(v) for k, v in tokens.items()}
    if any(not math.isfinite(v) for v in data.values()):
        raise ValueError("V56 nonfinite metric")
    for key in ("schemes", "faces", "verified", "nonsmooth", "kink_face_id", "kink_detected"):
        if data[key] != int(data[key]) or data[key] < 0:
            raise ValueError("V56 noninteger/count out of range")
    if not (data["schemes"] == 6 and data["faces"] == 634 and
            data["verified"] >= 1000 and data["kink_detected"] == 1 and
            0 <= data["kink_face_id"] < 634 and
            data["verified"] + data["nonsmooth"] == 6 * 634 and
            data["max_primal_error"] < 5e-12 and
            data["max_tangent_error"] < 2e-7 and
            data["max_cell_parity"] < 1e-10 and
            data["max_cell_tangent_error"] < 2e-6 and
            data["max_additivity"] < 1e-11):
        raise ValueError("V56 local face tangent proof fails strict audit limits")
    if re.search(r"(?im)^.*(?:qualified_global_jv|bfs_physically_converged|six_schur_qualified)\s*=\s*(?:true|1)\b", log):
        raise ValueError("V56 local face proof cannot qualify global Newton or BFS")
    return {
        "status": "LOCAL_CONVECTION_TANGENT_VERIFIED",
        "local_face_tangent_verified": True,
        "mesh_topology": "real_BFS152_V54_HDF5",
        "fields": "synthetic_controls_not_BFS_solution",
        "global_six_field_jv_qualified": False,
        "six_schur_global_newton_qualified": False,
        "bfs_physically_converged": False,
        "metrics": data,
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--log", type=pathlib.Path, required=True)
    p.add_argument("--output", type=pathlib.Path, required=True)
    args = p.parse_args()
    args.output.write_text(json.dumps(audit(args.log.read_text()), indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
