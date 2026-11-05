"""Independent negative and positive contracts for NACA force+mesh evidence."""
import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import naca_force_mesh_qualification as m  # noqa: E402


def mesh(similar=True):
    return [{"n_surface": n, "n_radial": n // 2,
             "radius": 20.0, "radial_power": 1.8,
             "minimum_first_layer_spacing": (1 / n if similar else 1 / (n ** 1.6)),
             "inverted_or_degenerate_cells": 0} for n in m.LEVELS]


class NacaForceMeshTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.folder = Path(self.tmp.name)
        self.certs = {}
        self.histories = {}
        for n, cd in zip(m.LEVELS, [0.11, 0.13, 0.138]):
            c = self.folder / f'cert_{n}.csv'
            with c.open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=m.REQUIRED)
                writer.writeheader()
                writer.writerow(dict(case=f'N{n}', h=1/n, converged='yes',
                                     res_mom=1e-11, res_mass=2e-11,
                                     cd=cd, cdp=cd * 0.4, cdv=cd * 0.6,
                                     cl=0, cm=0))
            h = self.folder / f'history_{n}.csv'
            h.write_text('iteration,momentum_rhs_relative,continuity_normalized,'
                         'mass_local_linf,nonfinite_mass_faces,nonfinite_velocities\n'
                         '1,1e-11,1e-12,2e-11,0,0\n', encoding='utf-8')
            self.certs[n] = c
            self.histories[n] = h

    def test_valid_evidence_produces_only_richardson_candidate(self):
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'GCI_CANDIDATE')
        self.assertEqual(report['gci']['validation_reference'], 'NOT_PERFORMED')

    def test_incompatible_near_wall_refinement_refuses_gci(self):
        report = m.certify(mesh(False), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')
        self.assertEqual(report['gci']['status'], 'NOT_QUALIFIED')

    def test_missing_n256_certificate_cannot_be_promoted(self):
        self.certs[256].unlink()
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')
        self.assertEqual(report['cases'][2]['status'], 'NOT_QUALIFIED')

    def test_stale_incorrect_certificate_fails(self):
        self.certs[128].write_text(self.certs[128].read_text().replace('N128', 'N129'))
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')
        self.assertTrue(any('wrong mesh label' in item for item in report['reasons']))

    def test_mismatched_physical_history_fails(self):
        path = self.histories[64]
        path.write_text(path.read_text().replace('1e-11,1e-12', '8e-9,1e-12'))
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')
        self.assertTrue(any('stale' in item for item in report['reasons']))

    def test_non_monotone_drag_is_not_richardson(self):
        path = self.certs[256]
        path.write_text(path.read_text().replace('0.138', '0.115'))
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')


if __name__ == '__main__':
    unittest.main()
