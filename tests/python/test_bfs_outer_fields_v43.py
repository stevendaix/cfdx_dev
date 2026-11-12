import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from bfs_outer_fields_v43 import audit


def specimen(cycle=True, missing=False, fake=False):
    values=[]
    for k in range(1,17):
        values.append(f"GLOBAL_NONLINEAR_COMMITTED_V13 iteration={k} merit={1.0/k:.9e} qualified=0")
        if k<5 or (missing and k==12):
            continue
        four=1e-6 if cycle else 0.10
        count=k-4 if cycle else 0
        flagged = int(count>=8)
        row=[f"OUTER_FIELD_PERIOD_V43 iteration={k}"]
        for key in ["Ux_four","Uy_four","Uz_four","p_gauge_four","phi_four",
                    "k_four","omega_four","nut_four"]:
            row.append(f"{key}={four}")
        row.extend([f"max_four={four}","max_previous=0.5",
                    "worst_channel=0","worst_index=2",
                    f"same_phase_count={count}",
                    f"period_four_candidate={flagged}",
                    f"physical_qualified={int(fake)}"])
        values.append(" ".join(row))
    return "\n".join(values)


class TestV43Parser(unittest.TestCase):
    def test_valid_four_cycle_is_diagnostic_only(self):
        result=audit(specimen(cycle=True),True)
        self.assertTrue(result["field_period_four_candidate"])
        self.assertFalse(result["physical_converged"])
        self.assertEqual(result["committed_states"],16)

    def test_nonperiodic_fields_rejected_when_required(self):
        self.assertFalse(audit(specimen(cycle=False))["field_period_four_candidate"])
        with self.assertRaises(ValueError):
            audit(specimen(cycle=False),True)

    def test_incomplete_trace_fails_closed(self):
        with self.assertRaises(ValueError):
            audit(specimen(missing=True))

    def test_qualification_claim_fails_closed(self):
        with self.assertRaises(ValueError):
            audit(specimen(fake=True))

    def test_merit_alone_not_sufficient(self):
        merit_only="\n".join(line for line in specimen().splitlines()
                             if line.startswith("GLOBAL_NONLINEAR_COMMITTED"))
        with self.assertRaises(ValueError):
            audit(merit_only)


if __name__ == "__main__":
    unittest.main()
