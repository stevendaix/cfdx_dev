"""Quantitative point-probe validation contracts for GUI post-processing."""
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Sequence
import math

@dataclass(frozen=True)
class ProbeSample:
    iteration: int
    time: float
    value: float

@dataclass(frozen=True)
class ProbeSeries:
    name: str
    point: tuple[float, ...]
    samples: tuple[ProbeSample, ...]

    def validate_monotonic(self) -> None:
        for a, b in zip(self.samples, self.samples[1:]):
            if b.iteration < a.iteration or b.time < a.time:
                raise ValueError("probe samples must be monotonic")

    def max_relative_error(self, reference: Sequence[float]) -> float:
        if len(reference) != len(self.samples):
            raise ValueError("reference and probe sample lengths differ")
        self.validate_monotonic()
        maximum = 0.0
        for sample, expected in zip(self.samples, reference):
            if not math.isfinite(sample.value) or not math.isfinite(float(expected)):
                raise ValueError("probe validation requires finite values")
            scale = max(abs(float(expected)), 1.0)
            maximum = max(maximum, abs(sample.value - float(expected)) / scale)
        return maximum

    def validate_against(self, reference: Sequence[float], *, atol: float = 1e-12, rtol: float = 1e-8) -> float:
        error = self.max_relative_error(reference)
        for sample, expected in zip(self.samples, reference):
            if abs(sample.value - float(expected)) > atol + rtol * abs(float(expected)):
                raise ValueError(
                    f"probe {self.name!r} exceeds tolerance at iteration {sample.iteration}"
                )
        return error
