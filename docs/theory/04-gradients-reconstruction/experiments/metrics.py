"""Reusable error metrics for gradient verification experiments."""

from __future__ import annotations

import numpy as np


def relative_l2_error(numerical: np.ndarray, exact: np.ndarray, weights: np.ndarray | None = None) -> float:
    numerical = np.asarray(numerical, dtype=float)
    exact = np.asarray(exact, dtype=float)
    if numerical.shape != exact.shape:
        raise ValueError("numerical and exact arrays must have the same shape")
    if weights is None:
        weights = np.ones(numerical.shape[0], dtype=float)
    weights = np.asarray(weights, dtype=float)
    if numerical.shape[0] != weights.shape[0]:
        raise ValueError("weights must match the number of points")
    numerator = np.sum(weights * np.sum((numerical - exact) ** 2, axis=1))
    denominator = np.sum(weights * np.sum(exact ** 2, axis=1))
    if denominator == 0.0:
        raise ValueError("exact field has zero norm")
    return float(np.sqrt(numerator / denominator))
