""""Headless monitor/probe contracts synchronized to solver iteration/time."""
from __future__ import annotations

import json
import math
from bisect import bisect_left
from dataclasses import dataclass
from pathlib import Path
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
        """Insert a sample in chronological order or merge an existing key.

        Solver stdout and stderr are consumed by independent threads, so their
        callbacks can arrive out of order even though each stream is ordered.
        ``upsert`` therefore keeps the public series ordered instead of treating
        callback arrival order as physical iteration order.
        """
        key = (sample.iteration, sample.time)
        keys = [(item.iteration, item.time) for item in self.samples]
        index = bisect_left(keys, key)
        if index < len(self.samples) and keys[index] == key:
            current = self.samples[index]
            values = dict(current.values)
            values.update(sample.values)
            self.samples[index] = MonitorSample(current.iteration, current.time, values)
            return
        self.samples.insert(index, sample)

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
            json.dumps(self.to_dict(), indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        return destination

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "MonitorSeries":
        """Reconstruct a monitor series while applying the normal monotonic contract."""
        name = payload.get("name")
        raw_samples = payload.get("samples")
        if not isinstance(name, str) or not name:
            raise ValueError("invalid monitor series name")
        if not isinstance(raw_samples, list):
            raise ValueError("invalid monitor series samples")
        series = cls(name, [])
        for raw in raw_samples:
            if not isinstance(raw, Mapping):
                raise ValueError("invalid monitor sample payload")
            iteration = raw.get("iteration")
            time_value = raw.get("time")
            values = raw.get("values")
            if not isinstance(iteration, int) or isinstance(iteration, bool) or iteration < 0:
                raise ValueError("invalid monitor iteration")
            if not isinstance(time_value, (int, float)) or isinstance(time_value, bool):
                raise ValueError("invalid monitor time")
            time_float = float(time_value)
            if not math.isfinite(time_float):
                raise ValueError("invalid monitor time")
            if not isinstance(values, Mapping):
                raise ValueError("invalid monitor values")
            normalized: dict[str, float] = {}
            for key, value in values.items():
                if not isinstance(key, str) or not key:
                    raise ValueError("invalid monitor value name")
                if not isinstance(value, (int, float)) or isinstance(value, bool):
                    raise ValueError("invalid monitor value")
                value_float = float(value)
                if not math.isfinite(value_float):
                    raise ValueError("invalid monitor value")
                normalized[key] = value_float
            series.append(MonitorSample(iteration, time_float, normalized))
        return series

    @classmethod
    def read_json(cls, path: str | Path) -> "MonitorSeries":
        """Load a persisted monitor series."""
        source = Path(path)
        return cls.from_dict(json.loads(source.read_text(encoding="utf-8"), parse_constant=lambda value: (_ for _ in ()).throw(ValueError(f"invalid JSON constant: {value}"))))

    def at(self, iteration: int) -> MonitorSample:
        for sample in reversed(self.samples):
            if sample.iteration == iteration:
                return sample
        raise KeyError(iteration)
