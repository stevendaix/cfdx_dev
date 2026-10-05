"""Point-probe model and CSV export contracts for the CFDX application layer.

The native solver already samples point probes and emits their time history to a
probe CSV (see ``cfdx::io::write_probe_csv`` and ``cfdx.probe_validation``). This
module owns the *definition* side of that contract — the probe catalogue the
Workbench editor populates and the solver argv it serializes to
(``--probe <name>:<x>,<y>,<z>:<field>`` / ``--probe-csv PATH``) — so the GUI and
the application layer never duplicate the probe schema that lives in C++.

It is a deliberately leaf module: only the standard library is imported, so it
can be used by :mod:`cfdx.case` and :mod:`cfdx.application.commands` without
introducing import cycles.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
import math


def nearest_probe(
    coordinates: Sequence[Sequence[float]],
    values: Sequence[float],
    point: Sequence[float],
) -> float:
    """Return the value at the nearest sample point."""
    if len(coordinates) != len(values) or not coordinates:
        raise ValueError("coordinates and values must be non-empty and aligned")
    if not point:
        raise ValueError("probe point must not be empty")
    dimensions = len(point)
    if any(len(coord) != dimensions for coord in coordinates):
        raise ValueError("coordinate dimensions must match the probe")
    index = min(
        range(len(coordinates)),
        key=lambda i: sum((coordinates[i][j] - point[j]) ** 2 for j in range(dimensions)),
    )
    return float(values[index])


class ProbeField(Enum):
    """A scalar field that may be sampled at a point probe.

    ``cli`` is the token accepted by the solver ``--probe`` option; ``unit`` is
    the SI unit reported in the probe CSV metadata. Members mirror
    ``cfdx::physics::IncompressibleProbeField`` so the Python and C++ schemas stay
    in lock-step.
    """

    PRESSURE = ("p", "Pa")
    U_X = ("u_x", "m/s")
    U_Y = ("u_y", "m/s")
    U_Z = ("u_z", "m/s")
    U_MAGNITUDE = ("u_mag", "m/s")

    def __init__(self, cli: str, unit: str) -> None:
        self.cli = cli
        self.unit = unit

    @classmethod
    def from_cli(cls, token: str) -> "ProbeField":
        token = (token or "").strip()
        for member in cls:
            if member.cli == token:
                return member
        raise ValueError(
            "--probe field must be one of " + ", ".join(member.cli for member in cls)
        )


@dataclass(frozen=True)
class Probe:
    """A single point probe: a labelled location and the field to sample."""

    name: str
    x: float
    y: float
    z: float
    field: ProbeField = ProbeField.PRESSURE

    def __post_init__(self) -> None:
        if not self.name or any(
            character in self.name for character in ",\r\n"
        ):
            raise ValueError(
                "probe name must be non-empty and free of commas/newlines: "
                f"{self.name!r}"
            )
        for axis_name, coordinate in (("x", self.x), ("y", self.y), ("z", self.z)):
            if not math.isfinite(coordinate):
                raise ValueError(
                    f"probe {self.name!r} has a non-finite {axis_name} coordinate"
                )

    def spec(self) -> str:
        """Serialize to the token accepted by the solver ``--probe`` option."""
        return f"{self.name}:{self.x},{self.y},{self.z}:{self.field.cli}"


@dataclass
class ProbeCatalog:
    """The mutable catalogue of probes attached to a case.

    The catalogue is the bridge between the Workbench probe editor and the
    solver argv: :meth:`specs` yields the ``--probe`` payloads and
    :meth:`csv_path` resolves the ``--probe-csv`` target. Probes are validated on
    insertion (finite location, non-empty name, no comma in the name) and their
    names must be unique, matching the constraints enforced by the native
    ``parse_probe``.
    """

    probes: list[Probe] = field(default_factory=list)

    def add(self, probe: Probe) -> None:
        self._validate_insertion(probe)
        self.probes.append(probe)

    def remove(self, name: str) -> bool:
        before = len(self.probes)
        self.probes = [probe for probe in self.probes if probe.name != name]
        return len(self.probes) != before

    def specs(self) -> list[str]:
        """Return the ``--probe`` payloads, one per probe, in declaration order."""
        return [probe.spec() for probe in self.probes]

    def names(self) -> list[str]:
        return [probe.name for probe in self.probes]

    def has(self, name: str) -> bool:
        return any(probe.name == name for probe in self.probes)

    def csv_path(self, project_path: str | Path) -> Path:
        """Resolve the ``--probe-csv`` target next to a saved case file."""
        case_path = Path(project_path)
        destination = case_path.parent / f"{case_path.stem}.probes.csv"
        destination.parent.mkdir(parents=True, exist_ok=True)
        return destination

    def _validate_insertion(self, probe: Probe) -> None:
        if self.has(probe.name):
            raise ValueError(f"probe name must be unique: {probe.name!r}")

    def validate(self) -> None:
        seen: set[str] = set()
        for probe in self.probes:
            probe.__post_init__()
            if probe.name in seen:
                raise ValueError(f"probe name must be unique: {probe.name!r}")
            seen.add(probe.name)
