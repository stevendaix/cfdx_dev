"""Headless monitor/probe contracts synchronized to solver iteration/time."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping
import json
from pathlib import Path


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

    def to_dict(self) -> dict[str, object]:
        """Return a stable, JSON-serializable monitor representation."""
        return {
            "name": self.name,
            "samples": [
                {
                    "iteration": sample.iteration,
                    "time": sample.time,
                    "values": dict(sample.values),
                }
                for sample in self.samples
            ],
        }

    def write_json(self, path: str | Path) -> Path:
        """Persist the monitor series without losing iteration/time provenance."""
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return destination

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "MonitorSeries":
        """Reconstruct a monitor series while applying the normal monotonic contract."""
        name = payload.get("name")
        raw_samples = payload.get("samples")
        if not isinstance(name, str) or not isinstance(raw_samples, list):
            raise ValueError("invalid monitor series payload")
        series = cls(name, [])
        for raw in raw_samples:
            if not isinstance(raw, Mapping):
                raise ValueError("invalid monitor sample payload")
            iteration = raw.get("iteration")
            time_value = raw.get("time")
            values = raw.get("values")
            if not isinstance(iteration, int) or isinstance(iteration, bool):
                raise ValueError("invalid monitor iteration")
            if not isinstance(time_value, (int, float)) or isinstance(time_value, bool):
                raise ValueError("invalid monitor time")
            if not isinstance(values, Mapping):
                raise ValueError("invalid monitor values")
            normalized = {str(key): float(value) for key, value in values.items()}
            series.append(MonitorSample(iteration, float(time_value), normalized))
        return series

    @classmethod
    def read_json(cls, path: str | Path) -> "MonitorSeries":
        """Load a persisted monitor series."""
        source = Path(path)
        return cls.from_dict(json.loads(source.read_text(encoding="utf-8")))

    def at(self, iteration: int) -> MonitorSample:
        for sample in reversed(self.samples):
            if sample.iteration == iteration:
                return sample
        raise KeyError(iteration)
