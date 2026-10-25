"""Independent acceptance examples for the T2 campaign's frozen criteria."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

path = Path(__file__).resolve().parents[2] / "scripts/t2_heated_channel_calibration.py"
spec = importlib.util.spec_from_file_location("t2_campaign", path)
campaign = importlib.util.module_from_spec(spec)
spec.loader.exec_module(campaign)


def observations():
    return [
        {"nx": nx, "ny": ny, "status": "PASS", "l2_poiseuille": 0.001, "err_rel": error}
        for (nx, ny), error in zip(campaign.MESHES, [0.16, 0.08, 0.04])
    ]


class CampaignVerdictTests(unittest.TestCase):
    def test_complete_campaign_meets_physical_targets(self):
        result = campaign.evaluate_campaign(observations(), True)
        self.assertTrue(result["qualification_passed"])
        self.assertEqual(result["verdict"], "physics_qualified")

    def test_execution_success_does_not_qualify_wrong_nusselt(self):
        data = observations()
        data[-1]["err_rel"] = 0.0993
        result = campaign.evaluate_campaign(data, True)
        self.assertTrue(result["execution_passed"])
        self.assertTrue(result["poiseuille_passed"])
        self.assertFalse(result["nusselt_passed"])
        self.assertFalse(result["qualification_passed"])
        self.assertEqual(result["verdict"], "execution_complete_physics_not_qualified")

    def test_poiseuille_gate_applies_to_every_mesh(self):
        data = observations()
        data[0]["l2_poiseuille"] = 0.006
        self.assertFalse(
            campaign.evaluate_campaign(data, True)["qualification_passed"]
        )

    def test_nan_infinity_negative_and_missing_metrics_cannot_pass(self):
        for metric in ("l2_poiseuille", "err_rel"):
            for invalid in (float("nan"), float("inf"), -1.0, None):
                with self.subTest(metric=metric, invalid=invalid):
                    for index in (0, -1):
                        data = observations()
                        data[index][metric] = invalid
                        self.assertFalse(
                            campaign.evaluate_campaign(data, True)[
                                "qualification_passed"
                            ]
                        )

    def test_missing_or_failed_mesh_does_not_qualify(self):
        for data in ([], observations()[:-1], observations()[1:]):
            self.assertFalse(
                campaign.evaluate_campaign(data, True)["execution_passed"]
            )
        data = observations()
        data[1]["status"] = "FAIL"
        self.assertFalse(
            campaign.evaluate_campaign(data, True)["qualification_passed"]
        )
        self.assertFalse(
            campaign.evaluate_campaign(observations(), False)["qualification_passed"]
        )

    def test_original_thresholds_are_inclusive_and_unchanged(self):
        self.assertEqual(campaign.L2_GATE, 0.005)
        self.assertEqual(campaign.ERR_TARGET, 0.05)
        data = observations()
        for entry in data:
            entry["l2_poiseuille"] = 0.005
        data[-1]["err_rel"] = 0.05
        self.assertTrue(
            campaign.evaluate_campaign(data, True)["qualification_passed"]
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
