"""Run the complete N2 reference gradient study and print JSON results."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from mesh import interior_indices, make_skewed_mesh
from quadratic_field import exact_gradient, field
from reconstruction import green_gauss, least_squares
from metrics import relative_l2_error
from convergence import observed_order


def study(method: str, levels: list[int]) -> dict[str, object]:
    errors = []
    hs = []
    condition_numbers = []

    for n in levels:
        mesh = make_skewed_mesh(n)
        values = field(mesh.centers[:, 0], mesh.centers[:, 1])
        exact = exact_gradient(mesh.centers[:, 0], mesh.centers[:, 1])

        if method == "green-gauss":
            numerical = green_gauss(mesh, values)
        elif method == "least-squares":
            numerical, cond = least_squares(mesh, values)
            condition_numbers.extend(cond[interior_indices(mesh)].tolist())
        elif method == "weighted-least-squares":
            numerical, cond = least_squares(mesh, values, weighted=True)
            condition_numbers.extend(cond[interior_indices(mesh)].tolist())
        else:
            raise ValueError(f"unknown method: {method}")

        idx = interior_indices(mesh)
        errors.append(relative_l2_error(numerical[idx], exact[idx], mesh.volumes[idx]))
        hs.append(1.0 / n)

    orders = [
        observed_order(errors[i], errors[i + 1], hs[i], hs[i + 1])
        for i in range(len(errors) - 1)
    ]
    result: dict[str, object] = {
        "method": method,
        "mesh_sizes": hs,
        "errors": errors,
        "observed_orders": orders,
    }
    if condition_numbers:
        result["max_condition_number"] = max(condition_numbers)
    return result


if __name__ == "__main__":
    levels = [8, 16, 32, 64]
    report = {
        "experiment": "N2 gradient reconstruction on a smooth skewed quadrilateral mesh",
        "qualification": False,
        "methods": [
            study("green-gauss", levels),
            study("least-squares", levels),
            study("weighted-least-squares", levels),
        ],
    }
    print(json.dumps(report, indent=2))
