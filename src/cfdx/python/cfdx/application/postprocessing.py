"""Headless post-processing capability registry shared by GUI, TUI and CLI.

The registry is metadata only: existing numerical/post-processing implementations remain
authoritative. GUI and other frontends enumerate this registry instead of maintaining
their own lists of derived fields, mesh-quality metrics or probe quantities.

The registry entries are *capability declarations*. A spec's ``availability`` is
"available" only when the quantity can actually be produced from a real production
dataset; "conditional" means it is computable but only from a richer artifact (a DAT
checkpoint) than the readable VTU result series. The registry never computes a value:
:func:`execute_derived` is the executable counterpart, and it is intentionally kept
separate so the registry stays a pure metadata contract that frontends can enumerate
without importing the computation layer.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Protocol

from ..derived import derived_vector_magnitude


class FieldDataset(Protocol):
    """Minimal read interface for the dataset a derived field is computed over."""

    def get(self, name: str) -> object | None: ...

    def names(self) -> tuple[str, ...]: ...


def _vector_field(dataset: object, name: str) -> tuple[tuple[float, ...], ...]:
    """Read a vector field as a tuple of per-cell component tuples.

    The production solver stores ``U`` as one cell field with three components
    per cell (``Field<double, Location::CELL>(n, "U", "m/s", 3)`` in
    ``apps/production_solver/main.cpp``), so a magnitude over the velocity must
    consume the vector array component-wise rather than three separately named
    scalars.
    """
    getter = getattr(dataset, "get", None)
    if not callable(getter):
        raise TypeError("dataset must expose a .get(name) accessor")
    raw = getter(name)
    if raw is None:
        raise KeyError(name)
    if not (
        isinstance(raw, (list, tuple)) and raw and isinstance(raw[0], (list, tuple))
    ):
        raise TypeError(f"field {name!r} is not a vector field")
    return tuple(tuple(float(value) for value in cell) for cell in raw)


@dataclass(frozen=True)
class _DerivedExecutor:
    stable_id: str
    source_field: str
    func: object


_DERIVED_EXECUTABLE: dict[str, _DerivedExecutor] = {
    "derived:u_magnitude": _DerivedExecutor(
        "derived:u_magnitude", "U", derived_vector_magnitude
    ),
}


def execute_derived(stable_id: str, dataset: FieldDataset) -> tuple[float, ...]:
    """Compute a declared derived field from a real dataset array.

    Only stable_ids declared in :data:`_DERIVED_EXECUTABLE` are accepted; every
    other id raises ``KeyError`` so a frontend cannot silently compute an
    unregistered quantity. The function performs no numerical tolerance checks:
    it is a pure arithmetic helper, and validation is the caller's responsibility.
    """
    spec = _DERIVED_EXECUTABLE.get(stable_id)
    if spec is None:
        raise KeyError(f"no executable derived field for: {stable_id!r}")
    return spec.func(_vector_field(dataset, spec.source_field))


_ALLOWED_CATEGORIES = frozenset(
    {
        "native_field",
        "derived_field",
        "mesh_quality",
        "probe",
        "integration",
        "engineering",
        "report",
    }
)
_ALLOWED_AVAILABILITY = frozenset({"available", "conditional", "unavailable"})
_ALLOWED_STATUS = frozenset(
    {"planned", "implemented", "verified", "validated", "qualified"}
)


@dataclass(frozen=True)
class PostProcessingSpec:
    """Stable metadata contract for one post-processing quantity/capability."""

    stable_id: str
    name: str
    category: str
    units: str | None
    inputs: tuple[str, ...]
    availability: str
    computation_owner: str
    validation_status: str

    def __post_init__(self) -> None:
        if not self.stable_id:
            raise ValueError("post-processing stable_id must not be empty")
        if not self.name:
            raise ValueError("post-processing name must not be empty")
        if self.category not in _ALLOWED_CATEGORIES:
            raise ValueError(f"unsupported post-processing category: {self.category!r}")
        if self.availability not in _ALLOWED_AVAILABILITY:
            raise ValueError(
                f"unsupported post-processing availability: {self.availability!r}"
            )
        if self.validation_status not in _ALLOWED_STATUS:
            raise ValueError(
                f"unsupported post-processing validation status: {self.validation_status!r}"
            )
        if not self.computation_owner:
            raise ValueError("post-processing computation_owner must not be empty")
        if len(set(self.inputs)) != len(self.inputs):
            raise ValueError("post-processing inputs must be unique")

    def to_dict(self) -> dict[str, object]:
        return {
            "stable_id": self.stable_id,
            "name": self.name,
            "category": self.category,
            "units": self.units,
            "inputs": list(self.inputs),
            "availability": self.availability,
            "computation_owner": self.computation_owner,
            "validation_status": self.validation_status,
        }


@dataclass
class PostProcessingRegistry:
    """Authoritative registry of post-processing capabilities."""

    _specs: dict[str, PostProcessingSpec] = field(default_factory=dict)

    def register(self, spec: PostProcessingSpec) -> None:
        if spec.stable_id in self._specs:
            raise ValueError(
                f"post-processing stable_id already registered: {spec.stable_id!r}"
            )
        self._specs[spec.stable_id] = spec

    def get(self, stable_id: str) -> PostProcessingSpec:
        try:
            return self._specs[stable_id]
        except KeyError as exc:
            raise KeyError(
                f"unknown post-processing capability: {stable_id!r}"
            ) from exc

    def all(self) -> tuple[PostProcessingSpec, ...]:
        return tuple(self._specs[key] for key in sorted(self._specs))

    def by_category(self, category: str) -> tuple[PostProcessingSpec, ...]:
        if category not in _ALLOWED_CATEGORIES:
            raise ValueError(f"unsupported post-processing category: {category!r}")
        return tuple(spec for spec in self.all() if spec.category == category)

    def available(self) -> tuple[PostProcessingSpec, ...]:
        return tuple(spec for spec in self.all() if spec.availability == "available")

    def extend(self, specs: Iterable[PostProcessingSpec]) -> None:
        for spec in specs:
            self.register(spec)


def build_post_processing_registry() -> PostProcessingRegistry:
    """Build the registry from already implemented CFDX post-processing owners."""

    registry = PostProcessingRegistry()

    registry.extend(
        (
            PostProcessingSpec(
                stable_id="field:p",
                name="Pressure",
                category="native_field",
                units="Pa",
                inputs=(),
                availability="available",
                computation_owner="cfdx::core::field",
                validation_status="verified",
            ),
            PostProcessingSpec(
                stable_id="field:U",
                name="Velocity",
                category="native_field",
                units="m/s",
                inputs=(),
                # The production solver's per-iteration VTU frames carry pressure
                # only (see ``apps/production_solver/main.cpp``: the
                # ``iteration_output_callback`` writes ``{{"p", pressure}}``).
                # Velocity is available through a loaded DAT checkpoint, not
                # through the readable result series, so the GUI must not
                # advertise it as a selectable VTU field.
                availability="conditional",
                computation_owner="cfdx::core::field",
                validation_status="verified",
            ),
            PostProcessingSpec(
                stable_id="derived:u_magnitude",
                name="Velocity magnitude",
                category="derived_field",
                units="m/s",
                inputs=("U",),
                availability="conditional",
                computation_owner="cfdx.python.application.postprocessing.execute_derived",
                validation_status="verified",
            ),
            PostProcessingSpec(
                stable_id="mesh_quality:skewness",
                name="Face skewness",
                category="mesh_quality",
                units=None,
                inputs=("face_centre", "owner_centre", "neighbour_centre", "Sf"),
                availability="available",
                computation_owner="cfdx::core::compute_face_quality",
                validation_status="verified",
            ),
            PostProcessingSpec(
                stable_id="mesh_quality:non_orthogonality",
                name="Face non-orthogonality",
                category="mesh_quality",
                units="rad",
                inputs=("face_centre", "owner_centre", "neighbour_centre", "Sf"),
                availability="available",
                computation_owner="cfdx::core::compute_face_quality",
                validation_status="verified",
            ),
            PostProcessingSpec(
                stable_id="probe:p",
                name="Point pressure probe",
                category="probe",
                units="Pa",
                inputs=("point", "p"),
                availability="available",
                computation_owner="cfdx.python.probe.ProbeField.PRESSURE",
                validation_status="validated",
            ),
            PostProcessingSpec(
                stable_id="probe:u_x",
                name="Point velocity X probe",
                category="probe",
                units="m/s",
                inputs=("point", "u_x"),
                availability="available",
                computation_owner="cfdx.python.probe.ProbeField.U_X",
                validation_status="verified",
            ),
            PostProcessingSpec(
                stable_id="probe:u_y",
                name="Point velocity Y probe",
                category="probe",
                units="m/s",
                inputs=("point", "u_y"),
                availability="available",
                computation_owner="cfdx.python.probe.ProbeField.U_Y",
                validation_status="verified",
            ),
            PostProcessingSpec(
                stable_id="probe:u_z",
                name="Point velocity Z probe",
                category="probe",
                units="m/s",
                inputs=("point", "u_z"),
                availability="available",
                computation_owner="cfdx.python.probe.ProbeField.U_Z",
                validation_status="verified",
            ),
            PostProcessingSpec(
                stable_id="probe:u_mag",
                name="Point velocity magnitude probe",
                category="probe",
                units="m/s",
                inputs=("point", "u_mag"),
                availability="available",
                computation_owner="cfdx.python.probe.ProbeField.U_MAGNITUDE",
                validation_status="verified",
            ),
        )
    )
    return registry
