"""CFDX I/O subpackage — unified solver converters and Gap Analysis."""

from cfdx.io.converter import (
    ConversionResult as CanonicalConversionResult,
)
from cfdx.io.converter import (
    SolverConverter,
    available_converters,
    convert,
    detect_solver,
    get_converter,
)
from cfdx.io.gap_analysis import Finding, GapAnalysis, GapAnalysisReport, Severity
from cfdx.io.interfaces import ConversionResult, SolverAdapter, SourceInfo
from cfdx.io.mapping_rules import MappingRules
from cfdx.io.schema import CaseSetup

__all__ = [
    "CanonicalConversionResult",
    "CaseSetup",
    "ConversionResult",
    "Finding",
    "GapAnalysis",
    "GapAnalysisReport",
    "MappingRules",
    "Severity",
    "SolverAdapter",
    "SolverConverter",
    "SourceInfo",
    "available_converters",
    "convert",
    "detect_solver",
    "get_converter",
]
