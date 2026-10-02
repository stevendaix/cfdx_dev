"""Executable linear-field gradient experiment.

The analytic gradient of phi(x,y) = 2*x - 3*y + 1 is constant.
This provides a basic exactness check for reconstruction experiments.
"""

from __future__ import annotations

import numpy as np


def field(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return 2.0 * x - 3.0 * y + 1.0


def exact_gradient() -> np.ndarray:
    return np.array([2.0, -3.0])


if __name__ == "__main__":
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    values = field(points[:, 0], points[:, 1])
    print("values:", values)
    print("exact gradient:", exact_gradient())
