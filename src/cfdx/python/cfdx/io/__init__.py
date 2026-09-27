"""CFDX I/O subpackage — unified solver converters and Gap Analysis."""

from cfdx.io.interfaces import SolverAdapter, ConversionResult, SourceInfo
from cfdx.io.gap_analysis import GapAnalysis, GapAnalysisReport, Severity, Finding
from cfdx.io.schema import CaseSetup
from cfdx.io.mapping_rules import MappingRules

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
]
