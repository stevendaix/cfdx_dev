"""Mapping rules loader and lookup for solver I/O converters.

Maps vendor-specific strings (Fluent zone types, OpenFOAM patch types,
SU2 markers, STAR-CCM+ region types) to CFDX internal taxonomies.

Rules are externalised in mapping_rules.yaml for easy updates without
recompilation.

Corresponds to the C++ header:
  src/cfdx/io/cfdx_io/mapping_rules.h
"""

from __future__ import annotations

import os
from pathlib import Path

from cfdx.io.schema import BCType, BCValueType


class MappingRules:
    """Maps vendor-specific strings to CFDX internal types.

    Load from a YAML file or use built-in defaults for Fluent, OpenFOAM,
    SU2, STAR-CCM+, and Code_Saturne.
    """

    # Default mappings (same as the C++ load_default)
    _DEFAULT_BC: dict[str, str] = {
        # Fluent
        "velocity-inlet": "INLET",
        "mass-flow-inlet": "INLET",
        "pressure-inlet": "INLET",
        "pressure-outlet": "PRESSURE_OUTLET",
        "outflow": "OUTLET",
        "wall": "WALL",
        "symmetry": "SYMMETRY",
        "periodic": "PERIODIC",
        "interface": "INTERFACE",
        "degenerate": "EMPTY",
        "internal": "INTERNAL",
        "slip-wall": "WALL",
        # OpenFOAM
        "fixedValue": "INLET",
        "fixedValueOutlet": "OUTLET",
        "zeroGradient": "INLET",
        "inlet": "INLET",
        "outlet": "OUTLET",
        "mappedWall": "WALL",
        "symmetryPlane": "SYMMETRY",
        "symmetry": "SYMMETRY",
        "empty": "EMPTY",
        "processor": "INTERFACE",
        "overset": "INTERFACE",
        # SU2
        " inlet": "INLET",
        "outlet": "PRESSURE_OUTLET",
        "wall": "WALL",
        "symmetry": "SYMMETRY",
        "periodic": "PERIODIC",
        "interface": "INTERFACE",
        "far": "OUTLET",
        "neumann": "OUTLET",
        # STAR-CCM+
        "Inlet": "INLET",
        "Pressure Outlet": "PRESSURE_OUTLET",
        "Outlet": "OUTLET",
        "Wall": "WALL",
        "Symmetry": "SYMMETRY",
        "Periodic": "PERIODIC",
        "Interface": "INTERFACE",
        "Free Stream": "INLET",
        "Opening": "PRESSURE_OUTLET",
        # Code_Saturne
        "far_field": "OUTLET",
        "opening": "PRESSURE_OUTLET",
    }

    _DEFAULT_BCVT: dict[str, str] = {
        "velocity-inlet": "FIXED",
        "mass-flow-inlet": "FIXED",
        "pressure-inlet": "FIXED",
        "pressure-outlet": "OUTLET_PRESSURE",
        "outflow": "ZERO_GRADIENT",
        "wall": "WALL_NO_SLIP",
        "symmetry": "ZERO_GRADIENT",
        "periodic": "MIXED",
        "interface": "MIXED",
        "fixedValue": "FIXED",
        "zeroGradient": "ZERO_GRADIENT",
        "mixed": "MIXED",
        "outletInlet": "OUTLET_PRESSURE",
        "inletOutlet": "FIXED",
        "inlet": "FIXED",
        "outlet": "ZERO_GRADIENT",
        "Inlet": "FIXED",
        "Pressure Outlet": "OUTLET_PRESSURE",
        "Outlet": "ZERO_GRADIENT",
        "Wall": "WALL_NO_SLIP",
        "Symmetry": "ZERO_GRADIENT",
        "Periodic": "MIXED",
        "Interface": "MIXED",
        "far_field": "FIXED",
        "opening": "OUTLET_PRESSURE",
    }

    _DEFAULT_TURBULENCE: dict[str, str] = {
        # Fluent
        "ke-realizable-viscous": "realizable_ke",
        "ke-standard-viscous": "k_epsilon",
        "ke-pseudo-viscous": "k_epsilon",
        "kw-viscous": "k_omega_sst",
        "mixing-length": "mixing_length",
        "transition-to-turbulence": "transition",
        "sst-ko-transport": "k_omega_sst",
        "sst-ko-Transport": "k_omega_sst",
        # OpenFOAM
        "kEpsilon": "k_epsilon",
        "kOmegaSST": "k_omega_sst",
        "kOmega": "k_omega",
        "realizableKE": "realizable_ke",
        "LRR": "reynolds_stress",
        "Boussinesq": "buoyant",
        # SU2
        "NONE": "laminar",
        "SA": "spalart_allmaras",
        "SST_SUTHERLAND": "k_omega_sst",
        "SA_NTS": "spalart_allmaras_nts",
        # STAR-CCM+
        "k-Epsilon": "k_epsilon",
        "k-Omega SST": "k_omega_sst",
        "Spalart-Allmaras": "spalart_allmaras",
        "Laminar": "laminar",
        "Realizable k-Epsilon": "realizable_ke",
        # Code_Saturne
        "KEPSILON": "k_epsilon",
        "KEPSIL": "k_epsilon",
        "KOMEGA": "k_omega",
        "K-OMEGA-SST": "k_omega_sst",
        "SPALART-ALLMARAS": "spalart_allmaras",
        "LES": "LES",
        "RANS": "RANS",
    }

    _DEFAULT_SCHEMES: dict[str, str] = {
        # Fluent
        "first-order": "first_order",
        "second-order-upwind": "second_order_upwind",
        "third-order-muscle": "third_order_muscle",
        "MUSCL": "MUSCL",
        "PRESTO!": "PRESTO",
        "STANDARD": "standard",
        "body-force-weighted": "body_force_weighted",
        "phase-coupled SIMPLE": "coupled",
        # OpenFOAM
        "Gauss linear": "linear",
        "Gauss linearUpwind": "linear_upwind",
        "Gauss upwind": "upwind",
        "Gauss limited": "limited",
        "Gauss LUST": "LUST",
        # SU2
        "GREEN_GAUSS_CELL_VOLUME": "green_gauss_cell",
        "LEAST_CONSERVATIVE": "least_squares",
    }

    _DEFAULT_PATCHES: dict[str, int] = {
        # WALL=0, INLET=1, OUTLET=2, SYMMETRY=3,
        # PERIODIC=4, INTERFACE=5, EMPTY=6, UNKNOWN=7
        "wall": 0,
        "mappedWall": 0,
        "velocity-inlet": 1,
        "pressure-inlet": 1,
        "mass-flow-inlet": 1,
        "inlet": 1,
        "pressure-outlet": 2,
        "outflow": 2,
        "outlet": 2,
        "open": 2,
        "symmetry": 3,
        "symmetryPlane": 3,
        "periodic": 4,
        "interface": 5,
        "processor": 5,
        "overset": 5,
        "internal": 5,
        "degenerate": 6,
        "empty": 6,
        "Wall": 0,
        "Inlet": 1,
        "Pressure Outlet": 2,
        "Outlet": 2,
        "Symmetry": 3,
        "Periodic": 4,
        "Interface": 5,
        "Far Field": 2,
        "far_field": 7,
        "opening": 2,
    }

    def __init__(self) -> None:
        self._bc_map: dict[str, BCType] = {}
        self._bcvt_map: dict[str, BCValueType] = {}
        self._turbulence_map: dict[str, str] = {}
        self._scheme_map: dict[str, str] = {}
        self._patch_map: dict[str, int] = {}

    @classmethod
    def load_default(cls) -> "MappingRules":
        rules = cls()
        for vendor, cfdx_name in cls._DEFAULT_BC.items():
            rules.add_bc_mapping(vendor, BCType[cfdx_name])
        for vendor, cfdx_name in cls._DEFAULT_BCVT.items():
            rules.add_bcvt_mapping(vendor, BCValueType[cfdx_name])
        for vendor, cfdx_name in cls._DEFAULT_TURBULENCE.items():
            rules.add_turbulence_mapping(vendor, cfdx_name)
        for vendor, cfdx_name in cls._DEFAULT_SCHEMES.items():
            rules.add_scheme_mapping(vendor, cfdx_name)
        for vendor, cfdx_type in cls._DEFAULT_PATCHES.items():
            rules.add_patch_mapping(vendor, cfdx_type)
        return rules

    @classmethod
    def load_from_file(cls, path: str) -> "MappingRules":
        rules = cls.load_default()

        if not os.path.exists(path):
            return rules

        try:
            import yaml  # type: ignore
        except ImportError:
            # Fallback: parse YAML manually (simplified)
            return cls._load_from_simple_yaml(path, rules)

        with open(path, "r") as f:
            data = yaml.safe_load(f)

        for vendor, cfdx_name in (data.get("bc") or {}).items():
            rules.add_bc_mapping(vendor, BCType[cfdx_name])
        for vendor, cfdx_name in (data.get("bcvt") or {}).items():
            rules.add_bcvt_mapping(vendor, BCValueType[cfdx_name])
        for vendor, cfdx_name in (data.get("turbulence") or {}).items():
            rules.add_turbulence_mapping(vendor, cfdx_name)
        for vendor, cfdx_name in (data.get("schemes") or {}).items():
            rules.add_scheme_mapping(vendor, cfdx_name)
        for vendor, cfdx_type in (data.get("patches") or {}).items():
            rules.add_patch_mapping(vendor, int(cfdx_type))

        return rules

    @classmethod
    def _load_from_simple_yaml(cls, path: str, rules: "MappingRules") -> "MappingRules":
        """Minimal YAML parser for the simple key:value format in mapping_rules.yaml."""
        section = ""
        with open(path, "r") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                if stripped.endswith(":") and not stripped.startswith('"'):
                    section = stripped[:-1].strip()
                    continue
                colon = stripped.find(":")
                if colon == -1:
                    continue
                key = stripped[:colon].strip()
                val = stripped[colon + 1 :].strip()
                # Strip quotes
                if len(key) >= 2 and key[0] in "\"'" and key[-1] == key[0]:
                    key = key[1:-1]

                if section == "bc" and val in BCType.__members__:
                    rules.add_bc_mapping(key, BCType[val])
                elif section == "bcvt" and val in BCValueType.__members__:
                    rules.add_bcvt_mapping(key, BCValueType[val])
                elif section == "turbulence":
                    rules.add_turbulence_mapping(key, val)
                elif section == "schemes":
                    rules.add_scheme_mapping(key, val)
                elif section == "patches":
                    try:
                        rules.add_patch_mapping(key, int(val))
                    except ValueError:
                        pass
        return rules

    def map_boundary_type(self, vendor_type: str) -> BCType:
        if vendor_type in self._bc_map:
            return self._bc_map[vendor_type]
        return BCType.UNKNOWN

    def map_bc_value_type(self, vendor_type: str) -> BCValueType:
        if vendor_type in self._bcvt_map:
            return self._bcvt_map[vendor_type]
        return BCValueType.UNKNOWN

    def map_turbulence(self, vendor_turbulence: str) -> str:
        return self._turbulence_map.get(vendor_turbulence, "")

    def map_scheme(self, vendor_scheme: str) -> str:
        return self._scheme_map.get(vendor_scheme, "")

    def map_patch_type(self, vendor_patch_type: str) -> int:
        return self._patch_map.get(vendor_patch_type, -1)

    def has_bc_mapping(self, vendor: str) -> bool:
        return vendor in self._bc_map

    def n_bc_mappings(self) -> int:
        return len(self._bc_map)

    def n_turbulence_mappings(self) -> int:
        return len(self._turbulence_map)

    def n_scheme_mappings(self) -> int:
        return len(self._scheme_map)

    # --- Modifiers ---
    def add_bc_mapping(self, vendor: str, cfdx_type: BCType) -> None:
        self._bc_map[vendor] = cfdx_type

    def add_bcvt_mapping(self, vendor: str, vtype: BCValueType) -> None:
        self._bcvt_map[vendor] = vtype

    def add_turbulence_mapping(self, vendor: str, cfdx_name: str) -> None:
        self._turbulence_map[vendor] = cfdx_name

    def add_scheme_mapping(self, vendor: str, cfdx_name: str) -> None:
        self._scheme_map[vendor] = cfdx_name

    def add_patch_mapping(self, vendor: str, cfdx_patch_type: int) -> None:
        self._patch_map[vendor] = cfdx_patch_type
