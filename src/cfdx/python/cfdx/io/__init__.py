"""CFDX I/O subpackage — unified solver converters and Gap Analysis."""

from cfdx.io.interfaces import SolverAdapter, ConversionResult, SourceInfo
from cfdx.io.gap_analysis import GapAnalysis, GapAnalysisReport, Severity, Finding
from cfdx.io.schema import CaseSetup
from cfdx.io.mapping_rules import MappingRules
from cfdx.io.converter import (
    ConversionResult as CanonicalConversionResult,
    SolverConverter,
    get_converter,
    convert,
    detect_solver,
    available_converters,
)

__all__ = [
    "SolverAdapter",
    "ConversionResult",
    "SourceInfo",
    "GapAnalysis",
    "GapAnalysisReport",
    "Severity",
    "Finding",
    "CaseSetup",
    "MappingRules",
    "CanonicalConversionResult",
    "SolverConverter",
    "get_converter",
    "convert",
    "detect_solver",
    "available_converters",
]
