import unittest

from scripts.audit_bfs_v68_flux_swap import audit, read_records


def row(retry):
    fields = dict(
        iteration=23,
        retry=retry,
        face=3,
        phi_carrier=0.01,
        phi_physical=0.011,
        face_omega_owner=0.001,
        face_omega_neighbour=-0.001,
        omega_phys_linf=0.4,
        omega_carrier_linf=0.5,
        omega_effect_fraction=-0.25,
        omega_delta_linf=0.1,
        omega_worst_phys=0,
        omega_worst_carrier=1,
        omega_worst_delta=12,
        omega_phys_cell0=0.4,
        omega_carrier_cell0=0.5,
        omega_phys_cell12=0.2,
        omega_carrier_cell12=0.1,
        k_phys_linf=0.6,
        k_carrier_linf=0.65,
        k_delta_linf=0.05,
        face_identity_gap=1e-13,
        source_diffusion_gap=2e-14,
        same_six_state=1,
        sole_phi_changed=1,
        committed=0,
        qualified=0,
    )
    return "BFS_SST_FLUX_SWAP_V68 " + " ".join(f"{k}={v}" for k, v in fields.items())


class TestV68Evidence(unittest.TestCase):
    def setUp(self):
        self.data = "\n".join(row(i) for i in (0, 8, 10))

    def test_complete(self):
        self.assertEqual(len(read_records(self.data)), 3)

    def test_missing_probe(self):
        with self.assertRaises(ValueError):
            read_records(self.data.replace(row(8), ""))

    def test_reject_non_finite(self):
        with self.assertRaises(ValueError):
            read_records(
                self.data.replace("omega_delta_linf=0.1", "omega_delta_linf=nan")
            )

    def test_reject_source_shift(self):
        with self.assertRaises(ValueError):
            read_records(
                self.data.replace(
                    "source_diffusion_gap=2e-14", "source_diffusion_gap=0.1"
                )
            )

    def test_reject_committed_state(self):
        with self.assertRaises(ValueError):
            read_records(self.data.replace("committed=0", "committed=1"))

    def test_reject_wrong_fraction(self):
        with self.assertRaises(ValueError):
            read_records(
                self.data.replace(
                    "omega_effect_fraction=-0.25", "omega_effect_fraction=0.25"
                )
            )

    def test_reject_nonmatching_reference(self):
        with self.assertRaises(ValueError):
            audit(self.data, "")


if __name__ == "__main__":
    unittest.main()
