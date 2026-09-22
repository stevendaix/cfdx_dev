"""Headless monitor/probe contracts synchronized to solver iteration/time."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping


@dataclass(frozen=True)
class MonitorSample:
    iteration: int
    time: float
    values: Mapping[str, float]


@dataclass
class MonitorSeries:
    name: str
    samples: list[MonitorSample]

    def append(self, sample: MonitorSample) -> None:
        if self.samples and (
            sample.iteration < self.samples[-1].iteration
            or sample.time < self.samples[-1].time
        ):
            raise ValueError("monitor samples must be monotonic")
        self.samples.append(sample)

    def upsert(self, sample: MonitorSample) -> None:
        """Append a sample or merge values for the current iteration/time."""
        if not self.samples:
            self.samples.append(sample)
            return
        current = self.samples[-1]
        if sample.iteration < current.iteration or sample.time < current.time:
            raise ValueError("monitor samples must be monotonic")
        if sample.iteration == current.iteration and sample.time == current.time:
            values = dict(current.values)
            values.update(sample.values)
            self.samples[-1] = MonitorSample(current.iteration, current.time, values)
            return
        self.samples.append(sample)

    def at(self, iteration: int) -> MonitorSample:
        for sample in reversed(self.samples):
            if sample.iteration == iteration:
                return sample
        raise KeyError(iteration)
