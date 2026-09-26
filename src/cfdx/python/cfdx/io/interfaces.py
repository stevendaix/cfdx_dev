"""Common interfaces for solver I/O converters.

Every external adapter implements SolverAdapter and produces a ConversionResult.
Adapters translate external solver data into the stable CFDX intermediate
representation (CaseSetup + mesh + fields + GapAnalysis report).

Corresponds to the C++ header:
  src/cfdx/io/cfdx_io/io_interface.h
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

import numpy as np
from pydantic import BaseModel

from cfdx.io.schema import CaseSetup, SourceInfo, BoundarySpec, MaterialSpec
from cfdx.io.gap_analysis import GapAnalysis


@dataclass
class ConversionResult:
    """Everything an adapter produces after conversion."""
    source: SourceInfo
    mesh: Optional[dict] = None  # meshio.Mesh or dict with points/cells
    scalar_fields: list[tuple[str, np.ndarray]] = field(default_factory=list)
    vec_fields: list[tuple[str, np.ndarray]] = field(default_factory=list)
    setup: Optional[CaseSetup] = None
    gap_report: GapAnalysis = field(default_factory=GapAnalysis)
    available_fields: dict[str, list[str]] = field(default_factory=dict)
    field_data_source: str = ""

    @property
    def success(self) -> bool:
        return not self.gap_report.has_blocking()

    @property
    def has_mesh(self) -> bool:
        return self.mesh is not None and self.mesh.get("points") is not None


class SolverAdapter(ABC):
    """Abstract base class for all solver adapters.

    Each external solver provides a concrete adapter. Adapters never touch
    CFDX numerical-core internals; they only populate the common CFDX
    representation (mesh, CaseSetup, fields, GapAnalysis).
    """

    @property
    @abstractmethod
    def solver_name(self) -> str:
        """Return the solver name, e.g. 'Fluent'."""
        ...

    @property
    @abstractmethod
    def format_name(self) -> str:
        """Return the format name, e.g. 'cas_dat'."""
        ...

    @abstractmethod
    def convert(self, case_path: Union[str, Path], result: ConversionResult) -> bool:
        """Parse the source files and populate result.

        Returns True if the conversion completed without blocking
        incompatibilities.
        """
        ...

    def import_results(self, results_path: str, result: ConversionResult) -> bool:
        """Optional: results-only import (after mesh already loaded)."""
        return False

    def detect_source(self, case_path: str, info: SourceInfo) -> bool:
        """Detect solver version and source info. Default: not implemented."""
        return False
