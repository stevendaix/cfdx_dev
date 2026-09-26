"""Pytest configuration: add the CFDX Python package to sys.path."""
import sys
import os

import pytest

_pkg_root = os.path.join(os.path.dirname(__file__), "..", "src", "cfdx", "python")
_pkg_root = os.path.abspath(_pkg_root)
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))

pytest.DATA_DIR = DATA_DIR
