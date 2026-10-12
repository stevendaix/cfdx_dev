"""Small utility for observed convergence order."""

from __future__ import annotations

from math import log


def observed_order(error_coarse: float, error_fine: float, h_coarse: float, h_fine: float) -> float:
    if min(error_coarse, error_fine, h_coarse, h_fine) <= 0.0:
        raise ValueError("errors and mesh sizes must be positive")
    return log(error_coarse / error_fine) / log(h_coarse / h_fine)


if __name__ == "__main__":
    print(observed_order(0.01, 0.0025, 0.1, 0.05))
