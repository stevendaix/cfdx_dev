#!/usr/bin/env python3
"""CFDX native Python package — I/O adapters and Gap Analysis reporting.

This module implements the unified solver I/O converter architecture
described in issue #425 of https://github.com/stevendaix/cfdx_dev.
"""

__version__ = "0.1.0"

from cfdx.io.interfaces import (
    SolverAdapter,
    ConversionResult,
    SourceInfo,
)
from cfdx.io.gap_analysis import (
    GapAnalysis,
    GapAnalysisReport,
    Severity,
    Finding,
)
from cfdx.io.schema import (
    CaseSetup,
    BoundarySpec,
    MaterialSpec,
    InitialCondition,
    NumericalScheme,
    MeshMetadata,
    BCType,
    BCValueType,
    CFDX_SCHEMA_VERSION,
)
from cfdx.io.mapping_rules import MappingRules
from cfdx.io.gap_analysis import GapAnalysisReport

__all__ = [
    "SolverAdapter",
    "ConversionResult",
    "SourceInfo",
    "GapAnalysis",
    "Severity",
    "Finding",
    "CaseSetup",
    "BoundarySpec",
    "MaterialSpec",
    "InitialCondition",
    "NumericalScheme",
    "MeshMetadata",
    "BCType",
    "BCValueType",
    "MappingRules",
    "GapAnalysisReport",
    "CFDX_SCHEMA_VERSION",
]
