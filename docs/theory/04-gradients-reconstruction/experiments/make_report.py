"""Generate a machine-readable convergence report."""

from __future__ import annotations

import json
from convergence import observed_order


def build_report(mesh_sizes: list[float], errors: list[float]) -> dict[str, object]:
    if len(mesh_sizes) != len(errors) or len(mesh_sizes) < 2:
        raise ValueError("at least two matching levels are required")
    orders = [
        observed_order(errors[i], errors[i + 1], mesh_sizes[i], mesh_sizes[i + 1])
        for i in range(len(errors) - 1)
    ]
    return {"mesh_sizes": mesh_sizes, "errors": errors, "observed_orders": orders}


if __name__ == "__main__":
    # Demonstration data only; not CFDX qualification evidence.
    report = build_report([0.25, 0.125, 0.0625], [4.0e-2, 1.0e-2, 2.5e-3])
    print(json.dumps(report, indent=2))
