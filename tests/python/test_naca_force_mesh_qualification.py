"""Independent negative and positive contracts for NACA force+mesh evidence."""
import csv
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
                         'mass_local_linf,nonfinite_mass_faces,nonfinite_velocities,'
                         'high_order_fraction,face_momentum_x_relative,'
                         'face_momentum_y_relative,face_momentum_z_relative,'
                         'face_momentum_x_linf,face_momentum_y_linf,face_momentum_z_linf,'
                         'convection_scheme\n'
                         '1,1e-11,1e-12,2e-11,0,0,1,1e-11,1e-11,0,1e-11,1e-11,0,TVD\n',
                         encoding='utf-8')
            self.certs[n] = c
            self.histories[n] = h

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

    def test_total_drag_gci_cannot_hide_nonmonotone_pressure_and_viscous_components(self):
        # Cd stays on an exact quadratic trend while Cd_pressure reverses.
        # A valid total Cd alone must not certify a pressure+skin-friction GCI.
        path = self.certs[128]
        with path.open(newline='') as stream:
            record = list(csv.DictReader(stream))[0]
        record['cdp'] = '0.06'
        record['cdv'] = str(float(record['cd']) - float(record['cdp']))
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=m.REQUIRED)
            writer.writeheader()
            writer.writerow(record)
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')
        self.assertTrue(any('cdp:' in reason for reason in report['reasons']))

    def test_independent_face_residual_blocks_false_matrix_convergence(self):
        path = self.histories[128]
        old = path.read_text()
        path.write_text(old.replace('1,1e-11,1e-12,2e-11,0,0,1,1e-11',
                                    '1,1e-11,1e-12,2e-11,0,0,1,0.002'))
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')
        self.assertTrue(any('independent local face momentum' in reason
                            for reason in report['reasons']))

    def test_missing_independent_residual_schema_blocks_gci(self):
        path = self.histories[64]
        path.write_text('iteration,momentum_rhs_relative,continuity_normalized,'
                        'mass_local_linf,nonfinite_mass_faces,nonfinite_velocities\n'
                        '1,1e-11,1e-12,2e-11,0,0\n')
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')
        self.assertTrue(any('incomplete physical history' in reason
                            for reason in report['reasons']))

    def test_first_order_upwind_cannot_claim_final_high_order_operator(self):
        path = self.histories[64]
        path.write_text(path.read_text().replace(',TVD\n', ',UPWIND\n'))
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')
        self.assertTrue(any('bounded second-order TVD' in reason
                            for reason in report['reasons']))

    def test_incomplete_final_convection_operator_blocks_gci(self):
        path = self.histories[256]
        path.write_text(path.read_text().replace(
            '1,1e-11,1e-12,2e-11,0,0,1,1e-11',
            '1,1e-11,1e-12,2e-11,0,0,0.5,1e-11'))
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')

    def test_negative_certificate_residual_is_not_a_small_residual(self):
        path = self.certs[64]
        path.write_text(path.read_text().replace('1e-11', '-1e-11'))
        history = self.histories[64]
        history.write_text(history.read_text().replace('1e-11', '-1e-11'))
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')

    def test_negative_continuity_in_history_cannot_be_certified(self):
        path = self.histories[128]
        path.write_text(path.read_text().replace('1e-12', '-1e-12'))
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')


if __name__ == '__main__':
    unittest.main()
        report = m.certify(mesh(), self.certs, self.histories)
        self.assertEqual(report['status'], 'BLOCKED')


if __name__ == '__main__':
    unittest.main()
