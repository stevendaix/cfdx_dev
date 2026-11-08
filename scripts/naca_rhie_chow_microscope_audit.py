"""Read the CFDX final-cell Rhie-Chow microscope; report, do not certify."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import re
from pathlib import Path

START = re.compile(r"=== CFDX MOMENTUM MICROSCOPE cell=(\d+) .* iteration=(\d+) ===")
FACE = re.compile(r"\bface=(\d+) owner=(\d+) neighbour=(-?\d+)\b")
VALUE = re.compile(r"\b(rhie_chow_expected|rhie_chow_actual|delta_rhie_chow)=([^\s]+)")
SUMMARY = re.compile(
    r"\b(rhie_chow_expected_linf|rhie_chow_actual_linf|delta_rhie_chow_linf)=([^\s]+)"
)


def parse_microscope(log: str) -> dict:
    starts = list(START.finditer(log))
    if len(starts) != 1:
        raise ValueError(
            "expected exactly one complete final-iteration momentum microscope"
        )
    selection = starts[0]
    end = log.find("=== END CFDX MOMENTUM MICROSCOPE ===", selection.end())
    if end < 0:
        raise ValueError("momentum microscope did not finish")
    fragment = log[selection.end() : end]
    face_rows = []
    summary = None
    for line in fragment.splitlines():
        matched = FACE.search(line)
        if matched:
            vals = dict(VALUE.findall(line))
            if set(vals) != {
                "rhie_chow_expected",
                "rhie_chow_actual",
                "delta_rhie_chow",
            }:
                raise ValueError("missing Rhie-Chow face microscope terms")
            parsed = {k: float(v) for k, v in vals.items()}
            if any(not math.isfinite(v) for v in parsed.values()):
                raise ValueError("nonfinite Rhie-Chow face term")
            err = parsed["rhie_chow_actual"] - parsed["rhie_chow_expected"]
            printed = parsed["delta_rhie_chow"]
            # The C++ microscope uses default output precision (six digits).
            if abs(printed - err) > 2e-5 * max(abs(printed), abs(err), 1e-12):
                raise ValueError("printed Rhie-Chow difference does not close")
            face_rows.append(
                {
                    "face": int(matched.group(1)),
                    "owner": int(matched.group(2)),
                    "neighbour": int(matched.group(3)),
                    **parsed,
                }
            )
        elif "rhie_chow_expected_linf=" in line:
            data = dict(SUMMARY.findall(line))
            if set(data) != set(
                (
                    "rhie_chow_expected_linf",
                    "rhie_chow_actual_linf",
                    "delta_rhie_chow_linf",
                )
            ):
                raise ValueError("incomplete Rhie-Chow final summary")
            summary = {k: float(v) for k, v in data.items()}
            if any(not math.isfinite(v) for v in summary.values()):
                raise ValueError("nonfinite Rhie-Chow summary")
    if not face_rows or summary is None:
        raise ValueError("missing Rhie-Chow face rows or final summary")
    internal = [r for r in face_rows if r["neighbour"] >= 0]
    if not internal:
        raise ValueError("no internal face for a Rhie-Chow operator check")
    localmax = max(abs(v["delta_rhie_chow"]) for v in internal)
    if localmax > 1.0001 * summary["delta_rhie_chow_linf"] + 1e-12:
        raise ValueError("Rhie-Chow summary smaller than internal face mismatch")
    return {
        "status": "DIAGNOSTIC_ONLY_NOT_PHYSICAL_QUALIFICATION",
        "cell": int(selection.group(1)),
        "iteration": int(selection.group(2)),
        "face_count": len(face_rows),
        "internal_face_count": len(internal),
        "summary": summary,
        "internal_max_absolute_difference": localmax,
        "internal_operator_discrepancy": [
            row
            for row in sorted(
                internal, key=lambda r: abs(r["delta_rhie_chow"]), reverse=True
            )
        ],
        "interpretation": (
            "rhie_chow_expected is the operator derived from pressure and "
            "the final reassembled momentum diagonal; rhie_chow_actual is "
            "HbyA predictor flux minus authoritative mass flux. Nonzero "
            "difference can indicate different states/pressure relaxation "
            "as well as a defect; it is not automatically an RC bug."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = parse_microscope(args.log.read_text(encoding="utf8"))
    data = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        if args.output.resolve() == args.log.resolve():
            parser.error("cannot overwrite input log")
        args.output.write_text(data, encoding="utf8")
    else:
        print(data, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
