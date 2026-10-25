"""Canonical import data must preserve unknown properties and parseable gravity."""
from __future__ import annotations

import unittest

from cfdx.io.schema import CaseSetup, MaterialSpec
from cfdx.io.adapters.fluent import _extract_cas_h5_reference_values


class ImportedThermalContracts(unittest.TestCase):
    def test_unspecified_material_does_not_become_air(self):
        material = MaterialSpec(name="unknown")
        for name in ("density", "dynamic_viscosity", "thermal_conductivity",
                     "specific_heat", "molecular_weight"):
            self.assertEqual(getattr(material, name), 0.0)
        self.assertEqual(material.eos_model, "")
        self.assertEqual(MaterialSpec.model_validate_json(material.model_dump_json()), material)

    def test_explicit_material_values_survive_serialization(self):
        material = MaterialSpec(name="declared", density=1.5, dynamic_viscosity=0.002,
                                thermal_conductivity=2.0, specific_heat=900.0)
        loaded = MaterialSpec.model_validate_json(material.model_dump_json())
        self.assertEqual(loaded, material)

    def test_fluent_gravity_uses_native_three_scalar_contract(self):
        setup = CaseSetup()
        _extract_cas_h5_reference_values(
            {"gravity?": True, "gravity/x": 1.5, "gravity/y": -9.81, "gravity/z": 0.25}, setup)
        self.assertEqual([float(value) for value in setup.gravity_vector.split()],
                         [1.5, -9.81, 0.25])
        _extract_cas_h5_reference_values({"gravity?": False}, setup)
        self.assertEqual([float(value) for value in setup.gravity_vector.split()], [0.0, 0.0, 0.0])


if __name__ == "__main__":
    unittest.main(verbosity=2)
