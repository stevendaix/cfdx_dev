"""Executable quadratic-field reference for gradient verification."""

from __future__ import annotations

import numpy as np


def field(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return x**2 + 2.0 * x * y + 3.0 * y**2


def exact_gradient(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.column_stack((2.0 * x + 2.0 * y, 2.0 * x + 6.0 * y))


if __name__ == "__main__":
    x = np.array([0.0, 0.5, 1.0])
    y = np.array([0.0, 0.5, 1.0])
    print(exact_gradient(x, y))
